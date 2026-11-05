#!/usr/bin/env python3
"""Convert a PROPKA summary and zreduce pKa JSON to residue-level comparison CSV."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

SUMMARY = re.compile(
    r"^\s*(ASP|GLU|HIS|LYS|CYS|TYR|N\+|C-)\s+(-?\d+)\s+(\S+)\s+(-?\d+(?:\.\d+)?)(?:\s+.*)?$"
)


def parse_propka(path: Path) -> dict[tuple[str, int, str], float]:
    values: dict[tuple[str, int, str], float] = {}
    in_summary = False
    for line in path.read_text(errors="replace").splitlines():
        if "SUMMARY OF THIS PREDICTION" in line:
            in_summary = True
            continue
        if not in_summary:
            continue
        match = SUMMARY.match(line)
        if match:
            comp, seq, chain, value = match.groups()
            values[(comp, int(seq), chain)] = float(value)
    if not values:
        raise ValueError(f"no PROPKA summary rows found in {path}")
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--zhyd", required=True, type=Path, help="zreduce pka JSON report"
    )
    parser.add_argument("--propka", required=True, type=Path, help="PROPKA .pka output")
    parser.add_argument("-o", "--output", type=Path)
    args = parser.parse_args()
    report = json.loads(args.zhyd.read_text())
    oracle = parse_propka(args.propka)
    rows = []
    site_map = {"N_TERM": "N+", "C_TERM": "C-"}
    for model in report["models"]:
        for site in model["sites"]:
            identity = site["id"]
            kind = site_map.get(identity["site"], identity["site"])
            key = (kind, int(identity["auth_seq_id"]), identity["chain_auth"])
            if key not in oracle:
                continue
            predicted = float(site["pka"]["predicted"])
            rows.append(
                {
                    "model": model["model"],
                    "comp_id": kind,
                    "auth_seq_id": key[1],
                    "chain": key[2],
                    "zhyd_pka": predicted,
                    "propka_pka": oracle[key],
                    "delta": predicted - oracle[key],
                }
            )
    handle = (
        args.output.open("w", newline="") if args.output else __import__("sys").stdout
    )
    try:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "model",
                "comp_id",
                "auth_seq_id",
                "chain",
                "zhyd_pka",
                "propka_pka",
                "delta",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)
    finally:
        if args.output:
            handle.close()


if __name__ == "__main__":
    main()
