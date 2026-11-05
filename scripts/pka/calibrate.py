#!/usr/bin/env python3
"""Fit `zreduce pka` parameters on the training split only.

Coordinate descent over a small, physically bounded parameter grid, minimizing
training-split MAE against experimental pKa values. Test and stress (holdout
family) rows are never read during fitting; evaluate them afterwards with
evaluate.py. Only the Python standard library is used.
"""

from __future__ import annotations

import argparse
import csv
import json
import tempfile
from pathlib import Path

import evaluate
import pkadr

# Search space: (name, candidate values). Bounds are kept physically plausible
# so the fit cannot, e.g., invert the sign of a term.
SEARCH_SPACE: list[tuple[str, list[float]]] = [
    ("burial_max_shift", [0.0, 0.5, 1.0, 1.5, 2.0, 3.0]),
    ("burial_surface_count", [6.0, 9.0, 12.0, 16.0, 20.0]),
    ("burial_core_count", [35.0, 45.0, 55.0, 70.0, 90.0]),
    ("hbond_strength", [0.0, 0.2, 0.35, 0.5, 0.7, 1.0]),
    ("hbond_max_shift", [0.75, 1.5, 2.0, 3.0]),
    ("hbond_cutoff", [3.5, 4.0, 4.5]),
    ("fixed_charge_scale", [0.0, 0.6, 1.2, 2.0, 3.0]),
    ("coupling_scale", [0.0, 0.5, 1.0, 2.0, 3.0, 4.0]),
]


def train_rows(canonical: Path, seed: str, train_fraction: float) -> list[dict]:
    with canonical.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    return [
        row
        for row in rows
        if evaluate.assign_split(row, seed, train_fraction) == "train"
    ]


def objective(
    rows: list[dict],
    params: dict,
    structures: Path,
    zreduce: Path,
    work: Path,
    kind: str = "mae",
) -> float:
    param_path = work / "params.json"
    param_path.write_text(json.dumps(params))
    joined, _ = pkadr.predict(rows, structures, zreduce, work / "reports", param_path)
    errors = [float(r["predicted_pka"]) - float(r["experimental_pka"]) for r in joined]
    if kind == "rmse":
        return (sum(e * e for e in errors) / len(errors)) ** 0.5
    return sum(abs(e) for e in errors) / len(errors)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(__doc__ or "").strip().splitlines()[0]
    )
    parser.add_argument("canonical_csv", type=Path)
    parser.add_argument("--structures", type=Path, required=True)
    parser.add_argument("--zreduce", type=Path, default=Path("zig-out/bin/zreduce"))
    parser.add_argument("--version", required=True, help="label for the fitted set")
    parser.add_argument("--seed", default="zhyd-pka-split-v1")
    parser.add_argument("--train-fraction", type=float, default=0.8)
    parser.add_argument("--sweeps", type=int, default=3)
    parser.add_argument("--objective", choices=("mae", "rmse"), default="mae")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()

    rows = train_rows(args.canonical_csv, args.seed, args.train_fraction)
    params: dict = {"version": args.version}
    history: list[dict] = []
    with tempfile.TemporaryDirectory() as tmp:
        work = Path(tmp)
        counter = iter(range(1_000_000))

        def score(candidate: dict) -> float:
            run_dir = work / str(next(counter))
            run_dir.mkdir()
            return objective(
                rows, candidate, args.structures, args.zreduce, run_dir, args.objective
            )

        best = score(params)
        history.append({"step": "initial", "train_objective": best})
        for sweep in range(args.sweeps):
            improved = False
            for name, values in SEARCH_SPACE:
                for value in values:
                    if params.get(name) == value:
                        continue
                    trial = {**params, name: value}
                    trial_score = score(trial)
                    if trial_score < best - 1e-4:
                        best, params, improved = trial_score, trial, True
                history.append(
                    {
                        "sweep": sweep,
                        "parameter": name,
                        "train_objective": round(best, 4),
                    }
                )
                print(
                    f"sweep {sweep} {name}: train {args.objective} {best:.4f}",
                    flush=True,
                )
            if not improved:
                break

    result = {
        "parameters": params,
        "fit": {
            "objective": f"train_split_{args.objective}",
            "seed": args.seed,
            "train_fraction": args.train_fraction,
            "train_rows": len(rows),
            "train_objective": round(best, 4),
            "history": history,
        },
    }
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(params, indent=2))


if __name__ == "__main__":
    main()
