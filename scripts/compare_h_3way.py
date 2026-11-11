#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "gemmi",
# ]
# ///
"""Three-way comparison: zhyd vs reduce2 vs ChimeraX.

Compares hydrogen placement across all 11 test structures using
distance-based matching with equivalent-H awareness.

Note: reduce2 uses different default bond lengths than ChimeraX/zhyd,
causing a systematic ~0.12A offset. This affects <0.1A% but not >=1.0A%
(which reflects orientation/flip differences).

Usage:
    ./scripts/compare_h_3way.py
    ./scripts/compare_h_3way.py --structures AF-P0A9J6
    ./scripts/compare_h_3way.py --verbose
    ./scripts/compare_h_3way.py --detail AF-P0A9J6  # per-residue detail
"""

import argparse
import os
import sys
from dataclasses import dataclass, field
from math import sqrt

import gemmi

SKIP_COMP_IDS = {"HOH", "UNK", "UNL"}

PROJECT_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_DIR = os.path.join(PROJECT_ROOT, "examples", "result")
CHIMERAX_DIR = os.path.join(RESULT_DIR, "chimerax_addh")
REDUCE2_DIR = os.path.join(RESULT_DIR, "reduce2")

STRUCTURES = [
    ("AF-C1P619-F1-model_v6", "AF model (16 res)"),
    ("AF-P0A9J6-F1-model_v6", "AF model (309 res)"),
    ("AF-P0DSH8-F1-model_v6", "AF model (16 res)"),
    ("AF-P22523-F1-model_v6", "AF model (1486 res)"),
    ("AF-P76347-F1-model_v6", "AF model (2339 res)"),
    ("1rqf", "PDB (4-chain, gaps)"),
    ("2cf8", "PDB (thrombin)"),
    ("2hnt", "PDB (gaps+inscodes)"),
    ("3rk2", "PDB (SNARE)"),
    ("6fys", "PDB (inscodes+altlocs)"),
    ("fold_test2_model_0", "PDB (DNA/RNA/glycan)"),
]

# Equivalent H pairs (same as compare_h.py)
EQUIV_PAIRS: list[tuple[str, str]] = [
    ("HB2", "HB3"),
    ("HG2", "HG3"),
    ("HD2", "HD3"),
    ("HE2", "HE3"),
    ("HZ2", "HZ3"),
    ("HA2", "HA3"),
    ("HG12", "HG13"),
    ("HD12", "HD13"),
    ("HH11", "HH12"),
    ("HH21", "HH22"),
    ("HD21", "HD22"),
    ("HE21", "HE22"),
    ("H21", "H22"),
    ("H41", "H42"),
    ("H61", "H62"),
    ("H71", "H72"),
    ("H72", "H73"),
    ("H71", "H73"),
    ("H2'", "H2''"),
    ("H5'", "H5''"),
    ("HG21", "HG22"),
    ("HG22", "HG23"),
    ("HG21", "HG23"),
    ("HB1", "HB2"),
    ("HB2", "HB3"),
    ("HB1", "HB3"),
    ("H1", "H2"),
    ("H2", "H3"),
    ("H1", "H3"),
    ("HE1", "HE2"),
    ("HE2", "HE3"),
    ("HE1", "HE3"),
]


def _build_equiv_map() -> dict[str, set[str]]:
    m: dict[str, set[str]] = {}
    for a, b in EQUIV_PAIRS:
        m.setdefault(a, {a}).add(b)
        m.setdefault(b, {b}).add(a)
    changed = True
    while changed:
        changed = False
        for v in m.values():
            for name in list(v):
                if name in m:
                    before = len(v)
                    v.update(m[name])
                    if len(v) > before:
                        changed = True
    return m


EQUIV_MAP = _build_equiv_map()


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


@dataclass
class CompareResult:
    n_matched: int = 0
    n_only_a: int = 0
    n_only_b: int = 0
    total_a: int = 0
    total_b: int = 0
    dists: list[float] = field(default_factory=list)
    # Per-residue details: (reskey, atom_a_name, atom_b_name, dist)
    details: list[tuple[str, str, str, float]] = field(default_factory=list)


def extract_h(path: str) -> tuple[dict[tuple, list[AtomPos]], dict[str, int]]:
    """Return ({residue_key: [AtomPos, ...]}, {skipped_comp: count})"""
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
        if comp in SKIP_COMP_IDS:
            skipped[comp] = skipped.get(comp, 0) + 1
            continue
        chain = row[3].strip()
        seq = row[4].strip()
        ins = row[9].strip() if row[9] not in (".", "?") else ""
        alt = row[8].strip() if row[8] not in (".", "?") else ""
        try:
            x, y, z = float(row[5]), float(row[6]), float(row[7])
        except ValueError:
            continue
        key = (chain, seq, comp, ins)
        result.setdefault(key, []).append(AtomPos(name=name, x=x, y=y, z=z, altloc=alt))
    return result, skipped


def match_atoms(
    atoms_a: list[AtomPos],
    atoms_b: list[AtomPos],
) -> tuple[list[tuple[AtomPos, AtomPos, float]], list[AtomPos], list[AtomPos]]:
    candidates: list[tuple[float, int, int]] = []
    for ia, aa in enumerate(atoms_a):
        for ib, bb in enumerate(atoms_b):
            if aa.altloc and bb.altloc and aa.altloc != bb.altloc:
                continue
            candidates.append((aa.dist(bb), ia, ib))
    candidates.sort()

    used_a: set[int] = set()
    used_b: set[int] = set()
    matched: list[tuple[AtomPos, AtomPos, float]] = []
    for d, ia, ib in candidates:
        if ia in used_a or ib in used_b:
            continue
        matched.append((atoms_a[ia], atoms_b[ib], d))
        used_a.add(ia)
        used_b.add(ib)

    unmatched_a = [a for i, a in enumerate(atoms_a) if i not in used_a]
    unmatched_b = [b for i, b in enumerate(atoms_b) if i not in used_b]
    return matched, unmatched_a, unmatched_b


def compare(h_a: dict, h_b: dict) -> CompareResult:
    result = CompareResult()
    result.total_a = sum(len(v) for v in h_a.values())
    result.total_b = sum(len(v) for v in h_b.values())

    all_keys = sorted(set(h_a.keys()) | set(h_b.keys()))
    for key in all_keys:
        aa = h_a.get(key, [])
        bb = h_b.get(key, [])

        if not aa and bb:
            result.n_only_b += len(bb)
            continue
        if aa and not bb:
            result.n_only_a += len(aa)
            continue

        matched, only_a, only_b = match_atoms(aa, bb)
        result.n_matched += len(matched)
        result.n_only_a += len(only_a)
        result.n_only_b += len(only_b)
        result.dists.extend(d for _, _, d in matched)

        reskey = f"{key[0]}/{key[2]} {key[1]}"
        for atom_a, atom_b, d in matched:
            result.details.append((reskey, atom_a.name, atom_b.name, d))

    return result


def dist_buckets(dists: list[float]) -> tuple[float, float, float, float]:
    """Return (<0.1A%, 0.1-0.5A%, 0.5-1.0A%, >=1.0A%) for matched pairs."""
    n = len(dists)
    if n == 0:
        return (0, 0, 0, 0)
    b1 = sum(1 for d in dists if d < 0.1) / n * 100
    b2 = sum(1 for d in dists if 0.1 <= d < 0.5) / n * 100
    b3 = sum(1 for d in dists if 0.5 <= d < 1.0) / n * 100
    b4 = sum(1 for d in dists if d >= 1.0) / n * 100
    return (b1, b2, b3, b4)


def find_file(stem: str, tool: str) -> str | None:
    """Find output file for a structure and tool."""
    if tool == "zhyd":
        p = os.path.join(RESULT_DIR, f"{stem}-zhyd.cif")
        return p if os.path.exists(p) else None
    elif tool == "reduce2":
        p = os.path.join(REDUCE2_DIR, f"{stem}_reduce2.cif")
        return p if os.path.exists(p) else None
    elif tool == "chimerax":
        for ext in ("_addh.cif.gz", "_addh.cif"):
            p = os.path.join(CHIMERAX_DIR, f"{stem}{ext}")
            if os.path.exists(p):
                return p
        return None
    return None


def print_summary_table(rows: list[dict]):
    """Print the main comparison table."""
    print(f"\n{'=' * 120}")
    skip_label = ", ".join(sorted(SKIP_COMP_IDS))
    print(f"Three-way H placement comparison: Z vs R vs C (skipped: {skip_label})")
    print(f"{'=' * 120}")

    # Header
    print(f"\n{'Structure':30s} │ {'Z vs C':^25s} │ {'R vs C':^25s} │ {'Z vs R':^25s}")
    print(
        f"{'':30s} │ {'<0.1A  >=1.0A  H-diff':^25s} │ {'<0.1A  >=1.0A  H-diff':^25s} │ {'<0.1A  >=1.0A  H-diff':^25s}"
    )
    print(f"{'─' * 30}─┼{'─' * 27}┼{'─' * 27}┼{'─' * 27}")

    for row in rows:
        stem = row["stem"]
        name = stem[:30]
        parts = []
        for pair in ["zc", "rc", "zr"]:
            r = row.get(pair)
            if r is None:
                parts.append(f"{'N/A':^25s}")
            else:
                b1, _, _, b4 = dist_buckets(r.dists)
                hdiff = r.total_a - r.total_b
                parts.append(f"{b1:5.1f}% {b4:5.1f}% {hdiff:+5d}")
        print(f"{name:30s} │ {parts[0]:^25s} │ {parts[1]:^25s} │ {parts[2]:^25s}")

    # Totals
    print(f"{'─' * 30}─┼{'─' * 27}┼{'─' * 27}┼{'─' * 27}")
    for pair, label in [
        ("zc", "Z vs C total"),
        ("rc", "R vs C total"),
        ("zr", "Z vs R total"),
    ]:
        all_dists = []
        total_a = total_b = 0
        for row in rows:
            r = row.get(pair)
            if r:
                all_dists.extend(r.dists)
                total_a += r.total_a
                total_b += r.total_b
        if all_dists:
            b1, b2, b3, b4 = dist_buckets(all_dists)
            n = len(all_dists)
            median = sorted(all_dists)[n // 2]
            print(
                f"  {label}: {n} matched, median={median:.3f}A, "
                f"<0.1A={b1:.1f}%, 0.1-0.5A={b2:.1f}%, 0.5-1.0A={b3:.1f}%, >=1.0A={b4:.1f}%, "
                f"H-diff={total_a - total_b:+d}"
            )


def print_detail(stem: str, rows: list[dict]):
    """Print per-residue detail for >=1.0A pairs."""
    for row in rows:
        if row["stem"] != stem:
            continue

        for pair, label in [
            ("zc", "zhyd vs ChimeraX"),
            ("rc", "reduce2 vs ChimeraX"),
            ("zr", "zhyd vs reduce2"),
        ]:
            r = row.get(pair)
            if r is None:
                continue
            bad = [(rk, na, nb, d) for rk, na, nb, d in r.details if d >= 1.0]
            if not bad:
                continue
            print(f"\n  {label}: {len(bad)} pairs >= 1.0A")
            bad.sort(key=lambda x: -x[3])
            for rk, na, nb, d in bad[:50]:
                print(f"    {rk:25s}  {na:5s} vs {nb:5s}  {d:.3f}A")
            if len(bad) > 50:
                print(f"    ... and {len(bad) - 50} more")
        break


def main():
    parser = argparse.ArgumentParser(description="Three-way H comparison")
    parser.add_argument(
        "--structures", "-s", help="Comma-separated structure stems to include"
    )
    parser.add_argument(
        "--verbose", "-v", action="store_true", help="Show per-structure detail"
    )
    parser.add_argument(
        "--detail", "-d", help="Show >=1.0A detail for a specific structure"
    )
    args = parser.parse_args()

    structures = STRUCTURES
    if args.structures:
        filters = [s.strip() for s in args.structures.split(",")]
        structures = [(s, d) for s, d in STRUCTURES if any(f in s for f in filters)]

    rows = []
    for stem, desc in structures:
        row = {"stem": stem, "desc": desc}
        files = {
            "z": find_file(stem, "zhyd"),
            "r": find_file(stem, "reduce2"),
            "c": find_file(stem, "chimerax"),
        }

        h_cache = {}
        for key, path in files.items():
            if path:
                h_cache[key], _ = extract_h(path)

        for pair, a_key, b_key in [
            ("zc", "z", "c"),
            ("rc", "r", "c"),
            ("zr", "z", "r"),
        ]:
            if a_key in h_cache and b_key in h_cache:
                row[pair] = compare(h_cache[a_key], h_cache[b_key])

        rows.append(row)

    print_summary_table(rows)

    if args.verbose:
        for row in rows:
            print_detail(row["stem"], rows)

    if args.detail:
        print_detail(args.detail, rows)


if __name__ == "__main__":
    main()
