#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "gemmi",
#     "numpy",
# ]
# ///
"""Compare hydrogen placement using nearest-neighbor coordinate matching.

Unlike name-based matching, this pairs each H in file A to its nearest H
in the same residue of file B. This avoids CCD naming convention mismatches
(e.g. zhyd H67 vs ChimeraX H3A in chlorophyll).

For standard residues with consistent naming, results should be identical
to name-based matching. For CCD-derived residues, this gives a more
accurate picture of actual placement quality.

Usage:
    ./scripts/chimerax_compare_h.py zhyd.cif chimerax.cif
    ./scripts/chimerax_compare_h.py zhyd.cif chimerax.cif --verbose
    ./scripts/chimerax_compare_h.py zhyd.cif chimerax.cif --summary-only
"""

import argparse
import sys
from dataclasses import dataclass, field
from math import sqrt

import gemmi
import numpy as np

SKIP_COMP_IDS = {"HOH", "UNK", "UNL"}


@dataclass
class AtomPos:
    name: str
    x: float
    y: float
    z: float
    altloc: str = ""

    def dist(self, other: "AtomPos") -> float:
        return sqrt(
            (self.x - other.x) ** 2 + (self.y - other.y) ** 2 + (self.z - other.z) ** 2
        )


def extract_h(path: str) -> tuple[dict[tuple, list[AtomPos]], dict[str, int]]:
    """Return ({(auth_chain, auth_seq, comp_id): [AtomPos, ...]}, {skipped_comp: count})"""
    doc = gemmi.cif.read(path)
    block = doc[0]
    table = block.find(
        "_atom_site.",
        [
            "type_symbol",
            "label_atom_id",
            "label_comp_id",
            "auth_asym_id",
            "auth_seq_id",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "label_alt_id",
            "pdbx_PDB_ins_code",
        ],
    )
    result: dict[tuple, list[AtomPos]] = {}
    skipped: dict[str, int] = {}
    for row in table:
        if row[0].strip() != "H":
            continue
        name = row[1].strip()
        comp = row[2].strip()
        chain = row[3].strip()
        seq = row[4].strip()
        ins = row[9].strip() if row[9] not in (".", "?") else ""
        alt = row[8].strip() if row[8] not in (".", "?") else ""
        try:
            x, y, z = float(row[5]), float(row[6]), float(row[7])
        except ValueError:
            continue
        if comp in SKIP_COMP_IDS:
            skipped[comp] = skipped.get(comp, 0) + 1
            continue
        key = (chain, seq, comp, ins)
        result.setdefault(key, []).append(AtomPos(name=name, x=x, y=y, z=z, altloc=alt))
    return result, skipped


def match_nearest(
    atoms_a: list[AtomPos], atoms_b: list[AtomPos]
) -> tuple[list[tuple[AtomPos, AtomPos, float]], list[AtomPos], list[AtomPos]]:
    """Match H atoms by nearest-neighbor within each residue.

    Greedy matching: sort all (A,B) pairs by distance, assign closest first.
    No name constraint — purely coordinate-based.
    """
    used_a: set[int] = set()
    used_b: set[int] = set()
    matched: list[tuple[AtomPos, AtomPos, float]] = []

    candidates: list[tuple[float, int, int]] = []
    for ia, aa in enumerate(atoms_a):
        for ib, bb in enumerate(atoms_b):
            if aa.altloc and bb.altloc and aa.altloc != bb.altloc:
                continue
            candidates.append((aa.dist(bb), ia, ib))
    candidates.sort()

    for d, ia, ib in candidates:
        if ia in used_a or ib in used_b:
            continue
        matched.append((atoms_a[ia], atoms_b[ib], d))
        used_a.add(ia)
        used_b.add(ib)

    unmatched_a = [a for i, a in enumerate(atoms_a) if i not in used_a]
    unmatched_b = [b for i, b in enumerate(atoms_b) if i not in used_b]
    return matched, unmatched_a, unmatched_b


@dataclass
class Stats:
    total_a: int = 0
    total_b: int = 0
    n_matched: int = 0
    n_only_a: int = 0
    n_only_b: int = 0
    dists: list[float] = field(default_factory=list)
    # Per-residue type stats
    type_dists: dict[str, list[float]] = field(default_factory=dict)
    res_diffs: list[tuple[str, int, int, list[str], list[str]]] = field(
        default_factory=list
    )
    skipped: dict[str, int] = field(default_factory=dict)


# Standard residue types for grouping
STANDARD_AA = {
    "ALA",
    "ARG",
    "ASN",
    "ASP",
    "CYS",
    "GLN",
    "GLU",
    "GLY",
    "HIS",
    "ILE",
    "LEU",
    "LYS",
    "MET",
    "PHE",
    "PRO",
    "SER",
    "THR",
    "TRP",
    "TYR",
    "VAL",
}
STANDARD_NUC = {"DA", "DC", "DG", "DT", "A", "C", "G", "U"}


def classify_residue(comp_id: str) -> str:
    if comp_id in STANDARD_AA:
        return "standard_aa"
    if comp_id in STANDARD_NUC:
        return "standard_nuc"
    return "ccd_other"


def compare_all(h_a: dict, h_b: dict) -> Stats:
    stats = Stats()
    stats.total_a = sum(len(v) for v in h_a.values())
    stats.total_b = sum(len(v) for v in h_b.values())

    all_keys = sorted(set(h_a.keys()) | set(h_b.keys()))

    for key in all_keys:
        aa = h_a.get(key, [])
        bb = h_b.get(key, [])
        comp_id = key[2]
        rtype = classify_residue(comp_id)

        if not aa and bb:
            stats.n_only_b += len(bb)
            label = f"{key[0]}/{key[2]} {key[1]}"
            stats.res_diffs.append((label, 0, len(bb), [], [b.name for b in bb]))
            continue
        if aa and not bb:
            stats.n_only_a += len(aa)
            label = f"{key[0]}/{key[2]} {key[1]}"
            stats.res_diffs.append((label, len(aa), 0, [a.name for a in aa], []))
            continue

        matched, only_a, only_b = match_nearest(aa, bb)
        stats.n_matched += len(matched)
        stats.n_only_a += len(only_a)
        stats.n_only_b += len(only_b)
        for _, _, d in matched:
            stats.dists.append(d)
            stats.type_dists.setdefault(rtype, []).append(d)

        if only_a or only_b or len(aa) != len(bb):
            label = f"{key[0]}/{key[2]} {key[1]}"
            stats.res_diffs.append(
                (
                    label,
                    len(aa),
                    len(bb),
                    [a.name for a in only_a],
                    [b.name for b in only_b],
                )
            )

    return stats


def print_report(
    stats: Stats, label_a: str, label_b: str, verbose: bool, summary_only: bool
):
    print(f"\n{'=' * 60}")
    skip_label = ", ".join(sorted(SKIP_COMP_IDS))
    print(f"Hydrogen Placement Comparison (coordinate-based, skipped: {skip_label})")
    print(f"  A: {label_a}")
    print(f"  B: {label_b}")
    print(f"{'=' * 60}")

    if stats.skipped:
        parts = [f"{comp}={n}" for comp, n in sorted(stats.skipped.items())]
        print(f"  Skipped H atoms (both files): {', '.join(parts)}")

    print(
        f"\n  Total H:  A={stats.total_a}  B={stats.total_b}  diff={stats.total_a - stats.total_b:+d}"
    )
    print(
        f"  Matched:  {stats.n_matched}  only-A: {stats.n_only_a}  only-B: {stats.n_only_b}"
    )

    if stats.dists:
        _print_dist_stats("All residues", stats.dists)

        # Per-type breakdown
        for rtype in ["standard_aa", "standard_nuc", "ccd_other"]:
            if rtype in stats.type_dists:
                labels = {
                    "standard_aa": "Standard AA",
                    "standard_nuc": "Standard nucleotides",
                    "ccd_other": "CCD-derived (other)",
                }
                _print_dist_stats(labels[rtype], stats.type_dists[rtype])

    if summary_only:
        if stats.res_diffs:
            print(f"\n  Residues with count differences: {len(stats.res_diffs)}")
        return

    only_a_res = [
        (r, ca, cb, oa, ob) for r, ca, cb, oa, ob in stats.res_diffs if cb == 0
    ]
    only_b_res = [
        (r, ca, cb, oa, ob) for r, ca, cb, oa, ob in stats.res_diffs if ca == 0
    ]
    count_diff = [
        (r, ca, cb, oa, ob)
        for r, ca, cb, oa, ob in stats.res_diffs
        if ca > 0 and cb > 0
    ]

    if only_a_res:
        print(f"\n  Residues with H only in A ({len(only_a_res)}):")
        for r, ca, _, oa, _ in only_a_res[:20]:
            print(f"    {r:25s}  A={ca}  [{','.join(oa[:5])}]")

    if only_b_res:
        print(f"\n  Residues with H only in B ({len(only_b_res)}):")
        for r, _, cb, _, ob in only_b_res[:20]:
            print(f"    {r:25s}  B={cb}  [{','.join(ob[:5])}]")

    if count_diff:
        print(f"\n  Residues with H count differences ({len(count_diff)}):")
        for r, ca, cb, oa, ob in count_diff[:50]:
            extra = ""
            if oa:
                extra += f"  only-A: {','.join(oa)}"
            if ob:
                extra += f"  only-B: {','.join(ob)}"
            print(f"    {r:25s}  A={ca:3d}  B={cb:3d}  diff={ca - cb:+d}{extra}")

    if verbose and stats.dists:
        bad = sum(1 for d in stats.dists if d >= 1.0)
        if bad:
            print(f"\n  Pairs with dist >= 1.0A: {bad}")


def _print_dist_stats(label: str, dists: list[float]):
    dists_sorted = sorted(dists)
    n = len(dists_sorted)
    if n == 0:
        return
    mean_d = sum(dists_sorted) / n
    median_d = dists_sorted[n // 2]
    p95 = dists_sorted[int(n * 0.95)]
    max_d = dists_sorted[-1]
    rmsd = sqrt(sum(d * d for d in dists_sorted) / n)
    n_close = sum(1 for d in dists_sorted if d < 0.1)
    n_mid = sum(1 for d in dists_sorted if 0.1 <= d < 0.5)
    n_far = sum(1 for d in dists_sorted if 0.5 <= d < 1.0)
    n_bad = sum(1 for d in dists_sorted if d >= 1.0)

    print(f"\n  {label} ({n} pairs):")
    print(
        f"    RMSD={rmsd:.3f}A  mean={mean_d:.3f}A  median={median_d:.3f}A  P95={p95:.3f}A  max={max_d:.3f}A"
    )
    print(f"    <0.1A: {n_close:5d} ({100 * n_close / n:5.1f}%)")
    print(f"    0.1-0.5A: {n_mid:3d} ({100 * n_mid / n:5.1f}%)")
    print(f"    0.5-1.0A: {n_far:3d} ({100 * n_far / n:5.1f}%)")
    print(f"    >=1.0A: {n_bad:5d} ({100 * n_bad / n:5.1f}%)")


def main():
    parser = argparse.ArgumentParser(
        description="Compare H placement (coordinate-based)"
    )
    parser.add_argument("file_a", help="First mmCIF (e.g. zhyd)")
    parser.add_argument("file_b", help="Second mmCIF (e.g. ChimeraX)")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--summary-only", "-s", action="store_true")
    args = parser.parse_args()

    h_a, skipped_a = extract_h(args.file_a)
    h_b, skipped_b = extract_h(args.file_b)
    stats = compare_all(h_a, h_b)
    all_skipped: dict[str, int] = {}
    for comp in {*skipped_a, *skipped_b}:
        all_skipped[comp] = skipped_a.get(comp, 0) + skipped_b.get(comp, 0)
    stats.skipped = all_skipped
    print_report(stats, args.file_a, args.file_b, args.verbose, args.summary_only)


if __name__ == "__main__":
    main()
