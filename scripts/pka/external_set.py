#!/usr/bin/env python3
"""Build an external pKa test set that shares no protein family with PKAD-R.

Inputs are downloaded separately and never committed:

- PKAD-3_V2.csv and time_stamped/blind_test.csv from the KaML-ESM repository
  (https://github.com/JanaShenLab/KaML-ESM, CC BY-NC 4.0);
- the PKAD-R JSON file (https://compbio.clemson.edu/pkad-r/);
- RCSB sequence clusters at 30 % identity
  (https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt).

A row is kept only if its PDB entry is not a PKAD-R structure (main or
alternative) and none of its polymer entities shares a 30 % cluster with any
entity of any PKAD-R structure. Entries missing from the cluster file are
dropped. Only wild-type rows of the six side-chain sites are kept; repeated
measurements of one residue are averaged. The output uses the canonical CSV
layout of `pkadr.py` with `holdout_family=external`, so `evaluate.py`
reports it as `stress:external`, never mixed with train/test.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

from pkadr import CANONICAL_FIELDS

ONE_LETTER = {"D": "ASP", "E": "GLU", "H": "HIS", "K": "LYS", "C": "CYS", "Y": "TYR"}
SITES = set(ONE_LETTER.values())
# PKAD-3 notes that keep the deposited sequence (cysteine redox states).
WILD_TYPE_NOTES = {"", "/", "reduced form", "oxidized form"}


def load_clusters(path: Path) -> dict[str, set[int]]:
    """PDB entry id -> set of 30 % cluster indices of its entities."""
    entry_clusters: dict[str, set[int]] = defaultdict(set)
    with path.open() as f:
        for index, line in enumerate(f):
            for entity in line.split():
                entry_clusters[entity.rsplit("_", 1)[0].upper()].add(index)
    return entry_clusters


def pkadr_entries(records: list[dict]) -> set[str]:
    entries = set()
    for record in records:
        entries.add(str(record["PDB"]).strip().upper())
        for alt in (
            str(record.get("Alternative_PDBs") or "").replace(";", ",").split(",")
        ):
            alt = alt.strip().upper()
            if len(alt) == 4:
                entries.add(alt)
    return entries


def candidate_rows(pkad3: list[dict], blind: list[dict]) -> list[dict]:
    """Normalize both sources to (pdb, chain, number, site, value, source)."""
    rows = []
    for r in pkad3:
        site = ONE_LETTER.get(r["Res_Name"].strip())
        if site is None or r["note"].strip() not in WILD_TYPE_NOTES:
            continue
        rows.append(
            {
                "pdb_id": r["PDB_ID"].strip().upper(),
                "chain": r["Chain"].strip(),
                "auth_seq_id": r["resid"].strip(),
                "site": site,
                "value": float(r["Expt_pKa"]),
                "source": "PKAD-3_V2",
            }
        )
    for r in blind:
        site = r["Res_Name"].strip().upper()
        if site not in SITES:
            continue
        rows.append(
            {
                "pdb_id": r["PDB_ID"].strip().upper(),
                "chain": r["Chain"].strip(),
                "auth_seq_id": r["Res_ID"].strip(),
                "site": site,
                "value": float(r["Expt_pKa"]),
                "source": "time_stamped",
            }
        )
    return [r for r in rows if r["auth_seq_id"].lstrip("-").isdigit()]


def build(
    pkad3: list[dict],
    blind: list[dict],
    pkadr: list[dict],
    clusters: dict[str, set[int]],
) -> tuple[list[dict[str, str]], dict]:
    excluded_entries = pkadr_entries(pkadr)
    excluded_clusters = set().union(*(clusters.get(e, set()) for e in excluded_entries))
    summary = {
        "candidates": 0,
        "pkadr_structure": 0,
        "homolog_of_pkadr": 0,
        "not_in_cluster_file": 0,
    }
    grouped: dict[tuple[str, str, str, str], list[dict]] = defaultdict(list)
    for row in candidate_rows(pkad3, blind):
        summary["candidates"] += 1
        entry = row["pdb_id"]
        if entry in excluded_entries:
            summary["pkadr_structure"] += 1
        elif entry not in clusters:
            summary["not_in_cluster_file"] += 1
        elif clusters[entry] & excluded_clusters:
            summary["homolog_of_pkadr"] += 1
        else:
            key = (entry, row["chain"], row["auth_seq_id"], row["site"])
            grouped[key].append(row)
    out = []
    for index, (key, rows) in enumerate(sorted(grouped.items()), start=1):
        entry, chain, number, site = key
        values = [r["value"] for r in rows]
        out.append(
            {
                "entry_id": f"ext{index}",
                "pdb_id": entry,
                "chain": chain,
                "auth_seq_id": number,
                "site": site,
                "experimental_pka": f"{sum(values) / len(values):.3f}",
                "classification": "Main",
                "mutant": "0",
                "cluster_id": f"rcsb30:{min(clusters[entry])}",
                "holdout_family": "external",
                "protein_name": "+".join(sorted({r["source"] for r in rows})),
            }
        )
    summary["sites"] = len(out)
    summary["entries"] = len({r["pdb_id"] for r in out})
    summary["clusters"] = len({r["cluster_id"] for r in out})
    return out, summary


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="") as f:
        return list(csv.DictReader(f))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build an external pKa test set disjoint from PKAD-R."
    )
    parser.add_argument("--pkad3", type=Path, required=True)
    parser.add_argument("--blind", type=Path)
    parser.add_argument("--pkadr", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args()

    pkadr_raw = args.pkadr.read_bytes()
    rows, summary = build(
        _read_csv(args.pkad3),
        _read_csv(args.blind) if args.blind else [],
        json.loads(pkadr_raw),
        load_clusters(args.clusters),
    )
    for name, path in (
        ("pkad3", args.pkad3),
        ("blind", args.blind),
        ("clusters", args.clusters),
    ):
        if path is not None:
            summary[f"{name}_sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    summary["pkadr_sha256"] = hashlib.sha256(pkadr_raw).hexdigest()
    with args.output.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=CANONICAL_FIELDS)
        writer.writeheader()
        writer.writerows(rows)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
