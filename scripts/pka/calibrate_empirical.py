#!/usr/bin/env python3
"""Recalibrate one site kind of the empirical engine (zhyd mode), train split only.

Grid search over the kind's model pKa and desolvation scale factor, passed to
`zreduce pka --engine empirical --empirical-mode zhyd --parameters`. Every
grid point is evaluated on the whole training split (structures that do not
contain the residue are run once and reused), and the per-row predictions are
cached, so the full fit and the leave-one-cluster-out cross-validation are
pure computations over the cache.

Four variants are fitted: the default point, model pKa only, scale only, and
both. Each variant's leave-one-cluster-out CV MAE (refit without the cluster,
predict its rows) is reported; the variant with the lowest CV MAE is
selected, so a parameter is adopted only when adding it improves CV MAE.
Test and stress rows are never read. Only the Python standard library is
used.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import tempfile
from collections.abc import Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import calibrate
import pkadr

ZHYD_ARGS = ["--empirical-mode", "zhyd"]
# Parity-table model pKa values ("Model values and charges").
DEFAULT_MODEL_PKA = {
    "ASP": 3.80,
    "GLU": 4.50,
    "HIS": 6.50,
    "LYS": 10.50,
    "ARG": 12.50,
    "CYS": 9.00,
    "TYR": 10.00,
    "N_TERM": 8.00,
    "C_TERM": 3.20,
}
# Residue names whose presence in a structure lets the kind's parameters
# change its report.
RESIDUE_NAMES = {"CYS": "CYS", "TYR": "TYR", "HIS": "HIS", "LYS": "LYS"}

Point = tuple[float, float]  # (model pKa, desolvation scale)


def grid(lo: float, hi: float, step: float) -> list[float]:
    n = round((hi - lo) / step)
    return [round(lo + i * step, 6) for i in range(n + 1)]


def overrides(version: str, kind: str, point: Point, fixed: dict) -> dict:
    """Override file for `point`, on top of already fixed kinds."""
    model = {**fixed.get("model_pka", {}), kind: point[0]}
    scale = {**fixed.get("desolvation_scale", {}), kind: point[1]}
    return {"version": version, "model_pka": model, "desolvation_scale": scale}


def abs_errors(
    predictions: Sequence[float | None], rows: list[dict]
) -> list[float | None]:
    return [
        None if p is None else abs(p - float(r["experimental_pka"]))
        for p, r in zip(predictions, rows)
    ]


def mae(errors: Sequence[float | None], subset: list[int]) -> float:
    values = [errors[i] for i in subset if errors[i] is not None]
    return sum(values) / len(values) if values else float("inf")


def fit(
    errors: Mapping[Point, Sequence[float | None]],
    points: list[Point],
    subset: list[int],
    default: Point,
) -> Point:
    """Grid point with the lowest MAE over `subset`; ties go to the point
    closest to `default`."""

    def key(point: Point) -> tuple[float, float]:
        distance = abs(point[0] - default[0]) + abs(point[1] - default[1])
        return (round(mae(errors[point], subset), 9), distance)

    return min(points, key=key)


def loco_cv(
    errors: Mapping[Point, Sequence[float | None]],
    points: list[Point],
    clusters: list[str],
    default: Point,
) -> tuple[list[float | None], dict[str, Point]]:
    """Leave-one-cluster-out errors: each cluster's rows are predicted with
    the point fitted on the other clusters."""
    out: list[float | None] = [None] * len(clusters)
    chosen: dict[str, Point] = {}
    for cluster in sorted(set(clusters)):
        held = [i for i, c in enumerate(clusters) if c == cluster]
        rest = [i for i, c in enumerate(clusters) if c != cluster]
        point = fit(errors, points, rest, default)
        chosen[cluster] = point
        for i in held:
            out[i] = errors[point][i]
    return out, chosen


def variants(
    models: list[float], scales: list[float], default: Point
) -> dict[str, list[Point]]:
    return {
        "default": [default],
        "model_pka": [(m, default[1]) for m in models],
        "desolvation_scale": [(default[0], s) for s in scales],
        "both": [(m, s) for m in models for s in scales],
    }


def select(
    errors: Mapping[Point, Sequence[float | None]],
    rows: list[dict],
    kind: str,
    models: list[float],
    scales: list[float],
    default: Point,
) -> dict:
    """Fit every variant, cross-validate it and pick the lowest CV MAE."""
    all_rows = list(range(len(rows)))
    kind_rows = [i for i, r in enumerate(rows) if r["site"] == kind]
    clusters = [r["cluster_id"] for r in rows]
    summary: dict[str, dict] = {}
    for name, points in variants(models, scales, default).items():
        point = fit(errors, points, all_rows, default)
        cv, _ = loco_cv(errors, points, clusters, default)
        summary[name] = {
            "model_pka": point[0],
            "desolvation_scale": point[1],
            "train_mae": round(mae(errors[point], all_rows), 4),
            "train_mae_kind": round(mae(errors[point], kind_rows), 4),
            "cv_mae": round(mae(cv, all_rows), 4),
            "cv_mae_kind": round(mae(cv, kind_rows), 4),
        }
    best = min(summary, key=lambda n: (summary[n]["cv_mae"], n != "default"))
    return {"variants": summary, "selected": best}


def outer_fold(cluster: str, seed: str, folds: int) -> int:
    digest = hashlib.sha256(f"{seed}\0outer\0{cluster}".encode()).digest()
    return int.from_bytes(digest[:8], "big") % folds


def nested_cv(
    errors: Mapping[Point, Sequence[float | None]],
    rows: list[dict],
    kind: str,
    models: list[float],
    scales: list[float],
    default: Point,
    folds: int,
    seed: str,
) -> dict:
    """Estimate the whole selection procedure: for each outer cluster fold,
    run `select` (fit + inner leave-one-cluster-out choice of variant) on the
    other folds and predict the held-out fold with the chosen point."""
    out: list[float | None] = [None] * len(rows)
    chosen = []
    fold_of = [outer_fold(r["cluster_id"], seed, folds) for r in rows]
    for fold in range(folds):
        inner = [i for i, f in enumerate(fold_of) if f != fold]
        held = [i for i, f in enumerate(fold_of) if f == fold]
        if not held or not inner:
            continue
        sub_errors = {p: [e[i] for i in inner] for p, e in errors.items()}
        result = select(
            sub_errors, [rows[i] for i in inner], kind, models, scales, default
        )
        best = result["variants"][result["selected"]]
        point = (best["model_pka"], best["desolvation_scale"])
        chosen.append(
            {"fold": fold, "variant": result["selected"], "point": list(point)}
        )
        for i in held:
            out[i] = errors[point][i]
    every = list(range(len(rows)))
    kind_rows = [i for i, r in enumerate(rows) if r["site"] == kind]
    return {
        "folds": folds,
        "nested_cv_mae": round(mae(out, every), 4),
        "nested_cv_mae_kind": round(mae(out, kind_rows), 4),
        "default_mae": round(mae(errors[default], every), 4),
        "default_mae_kind": round(mae(errors[default], kind_rows), 4),
        "chosen": chosen,
    }


def fitting_rows(path: Path, rows: str, seed: str, train_fraction: float) -> list[dict]:
    """`train`: the cluster train split; `all`: every non-holdout row."""
    if rows == "train":
        return calibrate.train_rows(path, seed, train_fraction)
    with path.open(newline="") as handle:
        return [r for r in csv.DictReader(handle) if not r["holdout_family"]]


def contains_residue(cif: Path, name: str) -> bool:
    return f" {name} " in cif.read_text()


def evaluate_grid(
    rows: list[dict],
    kind: str,
    points: list[Point],
    fixed: dict,
    structures: Path,
    zreduce: Path,
    work: Path,
    jobs: int,
) -> dict[Point, list[float | None]]:
    pdbs = sorted(
        {r["pdb_id"] for r in rows if (structures / f"{r['pdb_id']}.cif").exists()}
    )
    name = RESIDUE_NAMES.get(kind)
    affected = [
        p
        for p in pdbs
        if name is None or contains_residue(structures / f"{p}.cif", name)
    ]
    static = [p for p in pdbs if p not in affected]

    def run(point: Point, targets: list[str], tag: str) -> dict[str, dict]:
        out_dir = work / tag
        out_dir.mkdir(parents=True, exist_ok=True)
        param_path = out_dir / "params.json"
        param_path.write_text(json.dumps(overrides("calibration", kind, point, fixed)))

        def one(pdb: str) -> tuple[str, dict]:
            out = out_dir / f"{pdb}.json"
            command = pkadr.pka_command(
                structures / f"{pdb}.cif",
                out,
                zreduce,
                param_path,
                "empirical",
                ZHYD_ARGS,
            )
            subprocess.run(command, check=True, stderr=subprocess.DEVNULL)
            return pdb, json.loads(out.read_text())

        with ThreadPoolExecutor(jobs) as pool:
            return dict(pool.map(one, targets))

    default = (DEFAULT_MODEL_PKA[kind], 1.0)
    base_reports = run(default, static, "static")
    table: dict[Point, list[float | None]] = {}
    for n, point in enumerate(points):
        reports = {**base_reports, **run(point, affected, f"p{n}")}
        predictions: list[float | None] = []
        for row in rows:
            report = reports.get(row["pdb_id"])
            site = pkadr.match_site(report, row) if report else None
            predictions.append(
                None if site is None else float(site["pka"]["predicted"])
            )
        table[point] = predictions
        if n % 50 == 0:
            print(f"{kind}: {n + 1}/{len(points)} grid points", flush=True)
    return table


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(__doc__ or "").strip().splitlines()[0]
    )
    parser.add_argument("canonical_csv", type=Path)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--zreduce", type=Path, default=Path("zig-out/bin/zreduce"))
    parser.add_argument("--kind", required=True, choices=sorted(DEFAULT_MODEL_PKA))
    parser.add_argument("--seed", default="zhyd-pka-split-v1")
    parser.add_argument("--train-fraction", type=float, default=0.8)
    parser.add_argument(
        "--rows",
        choices=["train", "all"],
        default="train",
        help="fit on the train split, or on every non-holdout row (merged sets)",
    )
    parser.add_argument(
        "--nested-folds",
        type=int,
        default=0,
        help="also estimate the procedure by nested cluster CV with this many folds",
    )
    parser.add_argument(
        "--model-range", type=float, default=1.5, help="± around the default"
    )
    parser.add_argument("--model-step", type=float, default=0.1)
    parser.add_argument("--scale-min", type=float, default=0.2)
    parser.add_argument("--scale-max", type=float, default=1.5)
    parser.add_argument("--scale-step", type=float, default=0.05)
    parser.add_argument(
        "--fixed", type=Path, help="override JSON of kinds already fitted"
    )
    parser.add_argument(
        "--cache", type=Path, help="prediction cache (read if present, else written)"
    )
    parser.add_argument("-j", "--jobs", type=int, default=8)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()

    rows = fitting_rows(args.canonical_csv, args.rows, args.seed, args.train_fraction)
    fixed = json.loads(args.fixed.read_text()) if args.fixed else {}
    default = (DEFAULT_MODEL_PKA[args.kind], 1.0)
    models = grid(
        default[0] - args.model_range, default[0] + args.model_range, args.model_step
    )
    scales = grid(args.scale_min, args.scale_max, args.scale_step)
    points = sorted(
        {p for ps in variants(models, scales, default).values() for p in ps}
    )

    if args.cache and args.cache.exists():
        cached = json.loads(args.cache.read_text())
        table = {(m, s): v for m, s, v in cached}
    else:
        with tempfile.TemporaryDirectory() as tmp:
            table = evaluate_grid(
                rows,
                args.kind,
                points,
                fixed,
                args.structures,
                args.zreduce,
                Path(tmp),
                args.jobs,
            )
        if args.cache:
            args.cache.write_text(
                json.dumps([[m, s, v] for (m, s), v in table.items()])
            )

    errors = {p: abs_errors(v, rows) for p, v in table.items()}
    result = select(errors, rows, args.kind, models, scales, default)
    chosen = result["variants"][result["selected"]]
    nested = (
        nested_cv(
            errors,
            rows,
            args.kind,
            models,
            scales,
            default,
            args.nested_folds,
            args.seed,
        )
        if args.nested_folds
        else None
    )
    record = {
        "kind": args.kind,
        "selected": result["selected"],
        "model_pka": chosen["model_pka"],
        "desolvation_scale": chosen["desolvation_scale"],
        "variants": result["variants"],
        "nested_cv": nested,
        "fit": {
            "objective": f"{args.rows}_rows_mae",
            "rows": args.rows,
            "selection": "lowest leave-one-cluster-out CV MAE over the fitting rows",
            "seed": args.seed,
            "train_fraction": args.train_fraction,
            "fitting_rows": len(rows),
            "fitting_clusters": len({r["cluster_id"] for r in rows}),
            "kind_rows": sum(r["site"] == args.kind for r in rows),
            "grid": {
                "model_pka": [models[0], models[-1], args.model_step],
                "desolvation_scale": [scales[0], scales[-1], args.scale_step],
            },
            "fixed": fixed,
            "canonical_sha256": sha256(args.canonical_csv),
        },
    }
    args.output.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record["variants"], indent=2))
    print(f"selected: {result['selected']}")


if __name__ == "__main__":
    main()
