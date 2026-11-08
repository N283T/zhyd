#!/usr/bin/env python3
"""Compare two zhyd pKa reports site by site (oracle vs candidate engine).

Used for the parity gate of the empirical engine
(docs/specs/zhyd-pka-empirical-engine.md): run `zhyd pka` once with the
oracle engine and once with the candidate on the same input, then compare.
Exit status is 1 when a tolerance is exceeded or the site sets differ.
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

TERMS = ("model", "burial", "hydrogen_bond", "electrostatic")
SiteKey = tuple[int, str, int, str, str]


def site_key(model: int, identity: dict) -> SiteKey:
    return (
        model,
        identity["chain_auth"],
        int(identity["auth_seq_id"]),
        identity.get("ins_code") or "",
        identity["site"],
    )


def load_sites(path: Path) -> dict[SiteKey, dict]:
    report = json.loads(path.read_text())
    sites: dict[SiteKey, dict] = {}
    for model in report["models"]:
        for site in model["sites"]:
            sites[site_key(int(model["model"]), site["id"])] = site["pka"]
    return sites


def compare(
    oracle: dict[SiteKey, dict], candidate: dict[SiteKey, dict]
) -> tuple[list[dict], list[SiteKey], list[SiteKey]]:
    rows = []
    for key in sorted(oracle.keys() & candidate.keys()):
        ref, got = oracle[key], candidate[key]
        row = {
            "model": key[0],
            "chain": key[1],
            "auth_seq_id": key[2],
            "ins_code": key[3],
            "site": key[4],
            "oracle": ref["predicted"],
            "candidate": got["predicted"],
            "delta": got["predicted"] - ref["predicted"],
        }
        for term in TERMS:
            row[f"delta_{term}"] = got["terms"][term] - ref["terms"][term]
        rows.append(row)
    only_oracle = sorted(oracle.keys() - candidate.keys())
    only_candidate = sorted(candidate.keys() - oracle.keys())
    return rows, only_oracle, only_candidate


def summarize(
    rows: list[dict], tolerance: float, max_tolerance: float, fraction: float
) -> tuple[bool, str]:
    if not rows:
        return False, "no common sites"
    deltas = [abs(r["delta"]) for r in rows]
    within = sum(d <= tolerance for d in deltas) / len(deltas)
    worst = max(deltas)
    ok = within >= fraction and worst <= max_tolerance
    text = (
        f"sites={len(rows)} within_{tolerance}={within:.4f} "
        f"max_abs_delta={worst:.4f} "
        + " ".join(
            f"max_abs_delta_{t}={max(abs(r[f'delta_{t}']) for r in rows):.4f}"
            for t in TERMS
        )
    )
    return ok, text


def main() -> int:
    parser = argparse.ArgumentParser(description=(__doc__ or "").splitlines()[0])
    parser.add_argument(
        "--oracle", required=True, type=Path, help="report from the reference engine"
    )
    parser.add_argument(
        "--candidate",
        required=True,
        type=Path,
        help="report from the engine under test",
    )
    parser.add_argument("-o", "--output", type=Path, help="per-site CSV")
    parser.add_argument("--tolerance", type=float, default=0.02)
    parser.add_argument("--max-tolerance", type=float, default=0.1)
    parser.add_argument(
        "--fraction",
        type=float,
        default=0.99,
        help="required share of sites within --tolerance",
    )
    args = parser.parse_args()

    rows, only_oracle, only_candidate = compare(
        load_sites(args.oracle), load_sites(args.candidate)
    )
    if args.output:
        with args.output.open("w", newline="") as handle:
            fields = [
                "model",
                "chain",
                "auth_seq_id",
                "ins_code",
                "site",
                "oracle",
                "candidate",
                "delta",
            ]
            fields += [f"delta_{t}" for t in TERMS]
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    ok, text = summarize(rows, args.tolerance, args.max_tolerance, args.fraction)
    print(text)
    for label, keys in (
        ("only_oracle", only_oracle),
        ("only_candidate", only_candidate),
    ):
        for key in keys:
            print(f"{label}: {key}")
    return 0 if ok and not only_oracle and not only_candidate else 1


if __name__ == "__main__":
    sys.exit(main())
