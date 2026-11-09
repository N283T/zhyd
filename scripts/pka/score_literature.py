"""Score zhyd pKa reports against the literature evaluation set.

The set (`eval/literature-2026-11.csv`) mixes point values with upper bounds
(`value_qualifier=upper_bound`, e.g. "pKa < 3.5"), so the usual MAE harness
does not apply. For each report directory this prints the point MAE, how many
upper bounds hold, and how often the predicted protonation state at `--ph`
matches the experimental one; that last number is what placement uses.

    python3 scripts/pka/score_literature.py scripts/pka/eval/literature-2026-11.csv \
        parity=reports_parity/ zhyd=reports_zhyd/ --exclude-condition 0mM_NaCl

Each directory holds `<pdb_id>.json` reports from `zhyd pka`.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from pkadr import match_site


def experimental_state(row: dict[str, str], ph: float) -> bool | None:
    """Protonated at `ph` per experiment; None when the bound leaves it open."""
    value = float(row["experimental_pka"])
    if row["value_qualifier"] == "upper_bound":
        return False if value <= ph else None
    return value > ph


def score(rows: list[dict[str, str]], reports: Path, ph: float) -> dict:
    cache: dict[str, dict] = {}
    errors: list[float] = []
    bounds_met = bounds = 0
    agree = states = 0
    missing: list[str] = []
    low: list[str] = []
    flips: list[str] = []
    for row in rows:
        pdb = row["pdb_id"]
        if pdb not in cache:
            path = reports / f"{pdb}.json"
            cache[pdb] = json.loads(path.read_text()) if path.exists() else {}
        site = match_site(cache[pdb], row) if cache[pdb] else None
        if site is None:
            missing.append(row["entry_id"])
            continue
        predicted = float(site["pka"]["predicted"])
        value = float(row["experimental_pka"])
        if site["confidence"] == "low":
            low.append(row["entry_id"])
        if row["value_qualifier"] == "upper_bound":
            bounds += 1
            bounds_met += predicted < value
        else:
            errors.append(abs(predicted - value))
        expected = experimental_state(row, ph)
        if expected is not None:
            states += 1
            if (predicted > ph) == expected:
                agree += 1
            else:
                flips.append(
                    f"{row['entry_id']} exp={value} pred={predicted:.2f} confidence={site['confidence']}"
                )
    return {
        "points": len(errors),
        "mae": sum(errors) / len(errors) if errors else float("nan"),
        "max_error": max(errors, default=float("nan")),
        "bounds_met": bounds_met,
        "bounds": bounds,
        "state_agree": agree,
        "states": states,
        "low_confidence": low,
        "state_mismatches": flips,
        "missing": missing,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("csv", type=Path)
    parser.add_argument(
        "reports", nargs="+", help="LABEL=DIR with <pdb_id>.json reports"
    )
    parser.add_argument("--ph", type=float, default=7.0)
    parser.add_argument(
        "--exclude-condition",
        action="append",
        default=[],
        help="drop rows of this condition (e.g. 0mM_NaCl to keep one L01 condition)",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="list state mismatches and low-confidence sites",
    )
    args = parser.parse_args(argv)

    with args.csv.open(newline="") as handle:
        rows = [
            r
            for r in csv.DictReader(handle)
            if r["condition"] not in args.exclude_condition
        ]
    print(
        f"{'set':12s} {'MAE':>5s} {'max':>5s} {'n':>3s}  bounds  pH{args.ph:g} state  low"
    )
    for spec in args.reports:
        label, sep, directory = spec.partition("=")
        if not sep:
            parser.error(f"expected LABEL=DIR, got {spec!r}")
        s = score(rows, Path(directory), args.ph)
        print(
            f"{label:12s} {s['mae']:5.2f} {s['max_error']:5.2f} {s['points']:3d}"
            f"  {s['bounds_met']:2d}/{s['bounds']:<2d}   {s['state_agree']:2d}/{s['states']:<2d}"
            f"      {len(s['low_confidence'])}"
        )
        if s["missing"]:
            print(f"  missing: {', '.join(s['missing'])}", file=sys.stderr)
        if args.verbose:
            for line in s["state_mismatches"]:
                print(f"  state mismatch: {line}")
            for entry in s["low_confidence"]:
                print(f"  low confidence: {entry}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
