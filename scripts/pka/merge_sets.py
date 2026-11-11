#!/usr/bin/env python3
"""Merge canonical pKa CSVs into one fitting set with shared clusters.

Rows keep their entry ids (prefixed by set name) and experimental values.
Clusters are re-derived so that no protein family is split across folds:
two rows share a cluster when their original `cluster_id` is equal or when
their PDB entries share an RCSB 30 % sequence cluster. Holdout families
listed in `--keep-holdout` (default: snase) stay stress sets; every other
row becomes fitting data (`holdout_family` cleared).
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from external_set import load_clusters
from pkadr import CANONICAL_FIELDS


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, item: str) -> str:
        self.parent.setdefault(item, item)
        while self.parent[item] != item:
            self.parent[item] = self.parent[self.parent[item]]
            item = self.parent[item]
        return item

    def union(self, a: str, b: str) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.parent[max(ra, rb)] = min(ra, rb)


def merge(
    named_sets: list[tuple[str, list[dict[str, str]]]],
    clusters: dict[str, set[int]],
    keep_holdout: set[str],
) -> tuple[list[dict[str, str]], dict]:
    uf = UnionFind()
    rows = []
    for name, set_rows in named_sets:
        for row in set_rows:
            row = {f: row.get(f, "") for f in CANONICAL_FIELDS}
            row["entry_id"] = f"{name}:{row['entry_id']}"
            original = f"{name}:{row['cluster_id']}"
            uf.union(original, f"entry:{row['pdb_id'].upper()}")
            for index in clusters.get(row["pdb_id"].upper(), ()):
                uf.union(original, f"rcsb30:{index}")
            row["cluster_id"] = original
            if row["holdout_family"] not in keep_holdout:
                row["holdout_family"] = ""
            rows.append(row)
    for row in rows:
        row["cluster_id"] = uf.find(row["cluster_id"])
    summary = {
        "rows": len(rows),
        "fitting_rows": sum(not r["holdout_family"] for r in rows),
        "clusters": len({r["cluster_id"] for r in rows if not r["holdout_family"]}),
        "per_set": {name: len(s) for name, s in named_sets},
    }
    return rows, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="Merge canonical pKa CSVs.")
    parser.add_argument(
        "sets", nargs="+", help="NAME=PATH of canonical CSVs (pkadr.py layout)"
    )
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--keep-holdout", default="snase")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()

    named = []
    for spec in args.sets:
        name, path = spec.split("=", 1)
        with open(path, newline="") as f:
            named.append((name, list(csv.DictReader(f))))
    keep = {h for h in args.keep_holdout.split(",") if h}
    rows, summary = merge(named, load_clusters(args.clusters), keep)
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CANONICAL_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
