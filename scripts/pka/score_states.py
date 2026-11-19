"""Score the protonation states `run --ph` would apply against experiment.

`run --ph` (docs/specs/zhyd-protonation-pka.md, Phase 4) protonates a site
when its predicted pKa is at or above the pH and keeps every other site in its
default state (HIS, ASP and GLU neutral/deprotonated as placed; LYS charged;
CYS thiol). This compares, per site kind, how often the default states and the
pKa-derived states match the experimental state (protonated when the
experimental pKa is at or above the pH):

    python3 scripts/pka/score_states.py merged.csv reports/ --rows fit --ph 7
    python3 scripts/pka/score_states.py merged.csv reports/ --rows holdout \
        --group pdb_id --ph 6 --ph 7 --ph 8

Only sites whose state differs from the default can gain or lose: "right" and
"wrong" count those, and the net gain is right - wrong. A bootstrap over
`--group` (clusters by default) gives a 95 % interval for the net gain; with a
single group (e.g. the SNase stress set, one cluster) use `--group pdb_id`,
which treats related structures as independent and so understates the spread.

`reports/` holds `<pdb_id>.json` from `zhyd pka`. Standard library only.
"""

from __future__ import annotations

import argparse
import csv
import json
import random
from collections.abc import Sequence
from pathlib import Path

from pkadr import match_site

# Default state of each kind `run --ph` can change: True = protonated.
# Mirrors src/ph_states.zig (TYR, ARG and termini are never changed).
DEFAULT_PROTONATED = {
    "HIS": False,
    "ASP": False,
    "GLU": False,
    "LYS": True,
    "CYS": True,
}
PHYSICAL = (-5.0, 20.0)  # pka.isPhysicalPka


def applied_state(kind: str, predicted: float | None, ph: float) -> bool:
    """The state `run --ph` applies (True = protonated)."""
    if predicted is None or not PHYSICAL[0] <= predicted <= PHYSICAL[1]:
        return DEFAULT_PROTONATED[kind]
    return predicted >= ph


def site_calls(
    rows: Sequence[dict[str, str]], reports: Path
) -> tuple[list[dict], list[str]]:
    """Join rows to report predictions; returns (calls, missing entry ids)."""
    cache: dict[str, dict] = {}
    calls: list[dict] = []
    missing: list[str] = []
    for row in rows:
        if row["site"] not in DEFAULT_PROTONATED:
            continue
        pdb = row["pdb_id"]
        if pdb not in cache:
            path = reports / f"{pdb}.json"
            cache[pdb] = json.loads(path.read_text()) if path.exists() else {}
        site = match_site(cache[pdb], row) if cache[pdb] else None
        if site is None:
            missing.append(row["entry_id"])
            continue
        calls.append(
            {
                "row": row,
                "kind": row["site"],
                "experimental": float(row["experimental_pka"]),
                "predicted": float(site["pka"]["predicted"]),
            }
        )
    return calls, missing


def tally(calls: Sequence[dict], ph: float) -> dict[str, int]:
    out = {"n": 0, "default": 0, "ph": 0, "right": 0, "wrong": 0}
    for call in calls:
        truth = call["experimental"] >= ph
        default = DEFAULT_PROTONATED[call["kind"]]
        applied = applied_state(call["kind"], call["predicted"], ph)
        out["n"] += 1
        out["default"] += default == truth
        out["ph"] += applied == truth
        if applied != default:
            out["right" if applied == truth else "wrong"] += 1
    return out


def bootstrap_net(
    calls: Sequence[dict], ph: float, group: str, reps: int, seed: int
) -> tuple[float, float] | None:
    """95 % percentile interval of the net gain, resampling whole groups."""
    groups: dict[str, list[dict]] = {}
    for call in calls:
        groups.setdefault(call["row"][group], []).append(call)
    if len(groups) < 2:
        return None
    nets: dict[str, int] = {}
    for key, members in groups.items():
        t = tally(members, ph)
        nets[key] = t["right"] - t["wrong"]
    keys = sorted(nets)
    rng = random.Random(seed)
    samples = sorted(sum(nets[rng.choice(keys)] for _ in keys) for _ in range(reps))
    return samples[int(0.025 * reps)], samples[int(0.975 * reps) - 1]


def select_rows(path: Path, which: str) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    if which == "fit":
        return [r for r in rows if not r.get("holdout_family")]
    if which == "holdout":
        return [r for r in rows if r.get("holdout_family")]
    return rows


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("csv", type=Path, help="canonical or merged CSV")
    parser.add_argument("reports", type=Path, help="directory of <pdb_id>.json")
    parser.add_argument("--rows", choices=["fit", "holdout", "all"], default="all")
    parser.add_argument(
        "--ph", type=float, action="append", help="pH (repeatable; default 7)"
    )
    parser.add_argument("--group", default="cluster_id", help="bootstrap unit")
    parser.add_argument("--reps", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument(
        "--verbose", action="store_true", help="list every changed site"
    )
    args = parser.parse_args(argv)

    calls, missing = site_calls(select_rows(args.csv, args.rows), args.reports)
    groups = len({c["row"][args.group] for c in calls})
    print(f"{len(calls)} sites, {groups} groups ({args.group}), {len(missing)} missing")
    for ph in args.ph or [7.0]:
        print(f"\npH {ph:g}")
        print(
            f"{'kind':5s} {'n':>4s} {'default':>8s} {'--ph':>5s} {'right':>6s} {'wrong':>6s} {'net':>4s}  95% CI"
        )
        kinds = sorted({c["kind"] for c in calls}) + ["all"]
        for kind in kinds:
            subset = [c for c in calls if kind in ("all", c["kind"])]
            t = tally(subset, ph)
            ci = bootstrap_net(subset, ph, args.group, args.reps, args.seed)
            ci_text = "n/a" if ci is None else f"[{ci[0]:+.0f}, {ci[1]:+.0f}]"
            print(
                f"{kind:5s} {t['n']:4d} {t['default']:8d} {t['ph']:5d}"
                f" {t['right']:6d} {t['wrong']:6d} {t['right'] - t['wrong']:+4d}  {ci_text}"
            )
        if args.verbose:
            for c in calls:
                applied = applied_state(c["kind"], c["predicted"], ph)
                if applied != DEFAULT_PROTONATED[c["kind"]]:
                    ok = applied == (c["experimental"] >= ph)
                    print(
                        f"  {'right' if ok else 'WRONG'} {c['row']['entry_id']}"
                        f" {c['row']['pdb_id']} {c['kind']}{c['row']['auth_seq_id']}"
                        f" exp={c['experimental']:g} pred={c['predicted']:.2f}"
                    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
