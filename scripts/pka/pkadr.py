#!/usr/bin/env python3
"""Build a reproducible PKAD-R evaluation table for `zreduce pka`.

Subcommands (run in order):

    normalize  PKAD-R JSON -> canonical experimental CSV with cluster ids
    fetch      download the referenced PDB entries (mmCIF) from RCSB
    predict    run `zreduce pka` on each entry and join predictions

The PKAD-R file itself is not redistributed. Download it from the official
site (https://compbio.clemson.edu/pkad-r/) and cite
doi:10.1142/S2737416525500164. Only the Python standard library is used.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import subprocess
import sys
import urllib.request
from pathlib import Path

SITE_NAMES = {
    "ASP": "ASP",
    "GLU": "GLU",
    "HIS": "HIS",
    "LYS": "LYS",
    "CYS": "CYS",
    "TYR": "TYR",
    "N-term": "N_TERM",
    "C-term": "C_TERM",
}

# Proteins kept out of train/test and reported only as a stress test.
HOLDOUT_FAMILIES = {"SNase": "snase"}

CANONICAL_FIELDS = [
    "entry_id",
    "pdb_id",
    "chain",
    "auth_seq_id",
    "site",
    "experimental_pka",
    "classification",
    "mutant",
    "cluster_id",
    "holdout_family",
    "protein_name",
]

PREDICTION_FIELDS = CANONICAL_FIELDS + [
    "comp_id",
    "predicted_pka",
    "model_pka",
    "burial",
    "hydrogen_bond",
    "electrostatic",
    "coupling",
    "burial_fraction",
    "confidence",
]

NUMERIC = re.compile(r"^-?\d+(\.\d+)?$")


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


def _pdb_ids(field: str) -> list[str]:
    """Extract 4-character PDB ids from 'XXXX.A, YYYY.B' style fields."""
    return [token[:4].upper() for token in re.findall(r"[0-9][A-Za-z0-9]{3}", field)]


def normalize(records: list[dict[str, str]]) -> tuple[list[dict[str, str]], dict]:
    """Return canonical rows plus a summary of what was excluded and why."""
    clusters = UnionFind()
    for rec in records:
        pdb = rec["PDB"].strip().upper()
        if len(pdb) != 4:
            continue
        clusters.find(pdb)
        for other in _pdb_ids(rec.get("Sequence_Identity_>_30%", "")):
            clusters.union(pdb, other)
        for other in _pdb_ids(rec.get("Alternative_PDBs", "")):
            clusters.union(pdb, other)
    # Entries naming the same protein belong together even without a link.
    by_name: dict[str, str] = {}
    for rec in records:
        pdb = rec["PDB"].strip().upper()
        name = rec["Protein_Name"].strip().lower()
        if len(pdb) != 4 or not name:
            continue
        if name in by_name:
            clusters.union(pdb, by_name[name])
        else:
            by_name[name] = pdb

    rows: list[dict[str, str]] = []
    excluded: dict[str, int] = {}

    def skip(reason: str) -> None:
        excluded[reason] = excluded.get(reason, 0) + 1

    for rec in records:
        pdb = rec["PDB"].strip().upper()
        value = rec["Expt_pKa"].strip()
        site = SITE_NAMES.get(rec["ResName"].strip())
        warning = rec.get("Warning", "").strip()
        if len(pdb) != 4:
            skip("no_pdb")
        elif site is None:
            skip("unsupported_residue")
        elif not NUMERIC.match(value):
            skip("censored_or_nonnumeric_pka")
        elif warning in ("pKa: range or ~", "ResID NOT exist"):
            skip(f"warning:{warning}")
        elif not rec["ResID_in_PDB"].strip().lstrip("-").isdigit():
            skip("no_pdb_residue_id")
        else:
            protein = rec["Protein_Name"].strip()
            rows.append(
                {
                    "entry_id": rec["Index"],
                    "pdb_id": pdb,
                    "chain": rec["Chain"].strip(),
                    "auth_seq_id": rec["ResID_in_PDB"].strip(),
                    "site": site,
                    "experimental_pka": value,
                    "classification": rec["pKa_Classification"].strip(),
                    "mutant": "1" if rec["Mut_Pos"].strip() else "0",
                    "cluster_id": f"pkadr:{clusters.find(pdb)}",
                    "holdout_family": HOLDOUT_FAMILIES.get(protein, ""),
                    "protein_name": protein,
                }
            )
    summary = {
        "records": len(records),
        "kept": len(rows),
        "excluded": dict(sorted(excluded.items())),
        "clusters": len({row["cluster_id"] for row in rows}),
        "pdb_entries": len({row["pdb_id"] for row in rows}),
    }
    return rows, summary


def fetch(pdb_ids: list[str], out_dir: Path) -> dict[str, str]:
    out_dir.mkdir(parents=True, exist_ok=True)
    digests: dict[str, str] = {}
    for pdb in sorted(set(pdb_ids)):
        path = out_dir / f"{pdb}.cif"
        if not path.exists():
            url = f"https://files.rcsb.org/download/{pdb}.cif.gz"
            try:
                with urllib.request.urlopen(url, timeout=60) as response:
                    path.write_bytes(gzip.decompress(response.read()))
            except Exception as exc:  # network errors are reported, not fatal
                print(f"warning: {pdb}: {exc}", file=sys.stderr)
                continue
        digests[pdb] = hashlib.sha256(path.read_bytes()).hexdigest()
    return digests


def _site_index(report: dict) -> dict[tuple[str, int, str], dict]:
    index: dict[tuple[str, int, str], dict] = {}
    for site in report["models"][0]["sites"]:
        ident = site["id"]
        key = (ident["chain_auth"], int(ident["auth_seq_id"]), ident["site"])
        index.setdefault(key, site)
    return index


def pka_command(
    cif: Path,
    out: Path,
    zreduce: Path,
    parameters: Path | None = None,
    engine: str | None = None,
    extra_args: list[str] | None = None,
) -> list[str]:
    """The `zreduce pka` command that writes the model-1 report of `cif`."""
    command = [str(zreduce), "pka", str(cif), "--model", "1", "-o", str(out)]
    if parameters is not None:
        command += ["--parameters", str(parameters)]
    if engine is not None:
        command += ["--engine", engine]
    return command + (extra_args or [])


def match_site(report: dict, row: dict[str, str]) -> dict | None:
    """The report site of a canonical row (chain-less rows: first match)."""
    index = _site_index(report)
    seq = int(row["auth_seq_id"])
    if row["chain"]:
        return index.get((row["chain"], seq, row["site"]))
    candidates = [v for k, v in index.items() if k[1:] == (seq, row["site"])]
    return candidates[0] if candidates else None


def joined_row(row: dict[str, str], site: dict) -> dict[str, str]:
    terms = site["pka"]["terms"]
    return {
        **row,
        "comp_id": row["site"],
        "predicted_pka": str(site["pka"]["predicted"]),
        "model_pka": str(terms["model"]),
        "burial": str(terms["burial"]),
        "hydrogen_bond": str(terms["hydrogen_bond"]),
        "electrostatic": str(terms["electrostatic"]),
        "coupling": str(terms["coupling"]),
        "burial_fraction": str(site["pka"]["burial_fraction"]),
        "confidence": site["confidence"],
    }


def predict(
    rows: list[dict[str, str]],
    structures: Path,
    zreduce: Path,
    out_dir: Path,
    parameters: Path | None = None,
    engine: str | None = None,
    extra_args: list[str] | None = None,
) -> tuple[list[dict[str, str]], dict[str, int]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    reports: dict[str, dict | None] = {}
    joined: list[dict[str, str]] = []
    status = {"matched": 0, "missing_structure": 0, "unmatched_site": 0}
    for row in rows:
        pdb = row["pdb_id"]
        if pdb not in reports:
            cif = structures / f"{pdb}.cif"
            out = out_dir / f"{pdb}.json"
            if not cif.exists():
                reports[pdb] = None
            else:
                command = pka_command(cif, out, zreduce, parameters, engine, extra_args)
                subprocess.run(command, check=True, stderr=subprocess.DEVNULL)
                reports[pdb] = json.loads(out.read_text())
        report = reports[pdb]
        if report is None:
            status["missing_structure"] += 1
            continue
        site = match_site(report, row)
        if site is None:
            status["unmatched_site"] += 1
            continue
        status["matched"] += 1
        joined.append(joined_row(row, site))
    return joined, status


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _write_csv(path: Path, fields: list[str], rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(__doc__ or "").strip().splitlines()[0]
    )
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("normalize")
    p.add_argument("pkadr_json", type=Path)
    p.add_argument("-o", "--output", type=Path, required=True)

    p = sub.add_parser("fetch")
    p.add_argument("canonical_csv", type=Path)
    p.add_argument("--structures", type=Path, required=True)

    p = sub.add_parser("predict")
    p.add_argument("canonical_csv", type=Path)
    p.add_argument("--structures", type=Path, required=True)
    p.add_argument("--zreduce", type=Path, default=Path("zig-out/bin/zreduce"))
    p.add_argument("--reports", type=Path, required=True)
    p.add_argument("--parameters", type=Path, help="zreduce pka parameter JSON")
    p.add_argument("--engine", help="zreduce pka engine (builtin, empirical, propka)")
    p.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        metavar="ARG",
        help="extra argument appended to the zreduce pka command (repeatable)",
    )
    p.add_argument("-o", "--output", type=Path, required=True)

    args = parser.parse_args()
    if args.command == "normalize":
        raw = args.pkadr_json.read_bytes()
        rows, summary = normalize(json.loads(raw))
        summary["source_sha256"] = hashlib.sha256(raw).hexdigest()
        _write_csv(args.output, CANONICAL_FIELDS, rows)
        print(json.dumps(summary, indent=2))
    elif args.command == "fetch":
        ids = [row["pdb_id"] for row in _read_csv(args.canonical_csv)]
        digests = fetch(ids, args.structures)
        manifest = args.structures / "manifest.json"
        manifest.write_text(json.dumps(digests, indent=2, sort_keys=True) + "\n")
        print(f"{len(digests)} structures; manifest: {manifest}")
    elif args.command == "predict":
        joined, status = predict(
            _read_csv(args.canonical_csv),
            args.structures,
            args.zreduce,
            args.reports,
            args.parameters,
            args.engine,
            args.extra_arg,
        )
        _write_csv(args.output, PREDICTION_FIELDS, joined)
        print(json.dumps(status, indent=2))


if __name__ == "__main__":
    main()
