#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "gemmi",
# ]
# ///
"""Compare hydrogen placement between two mmCIF files with equivalent-H awareness.

Methylene H pairs (HB2/HB3, HG2/HG3, etc.) and symmetric groups (HH11/HH12,
HH21/HH22, NH2 HD21/HD22, HE21/HE22) are matched by minimum distance rather
than by name, since naming conventions differ between tools.

Usage:
    ./scripts/compare_h_v2.py zhyd.cif chimerax.cif
    ./scripts/compare_h_v2.py zhyd.cif chimerax.cif --verbose
    ./scripts/compare_h_v2.py zhyd.cif chimerax.cif --summary-only
"""

import argparse
from dataclasses import dataclass, field
from math import sqrt

import gemmi

SKIP_COMP_IDS = {"HOH", "UNK", "UNL"}

# Pairs of equivalent hydrogen names that may be swapped between tools.
# These are methylene pairs, symmetric amine/guanidinium pairs, etc.
EQUIV_PAIRS: list[tuple[str, str]] = [
    # Methylene CH2
    ("HB2", "HB3"),
    ("HG2", "HG3"),
    ("HD2", "HD3"),
    ("HE2", "HE3"),
    ("HZ2", "HZ3"),
    ("HA2", "HA3"),
    ("HG12", "HG13"),
    ("HD12", "HD13"),
    # Symmetric groups (guanidinium, amide, etc.)
    ("HH11", "HH12"),
    ("HH21", "HH22"),
    ("HD21", "HD22"),
    ("HE21", "HE22"),
    # Nucleotide base NH2 (symmetric amine)
    ("H21", "H22"),  # guanine/DG NH2
    ("H41", "H42"),  # cytosine/DC NH2
    ("H61", "H62"),  # adenine/DA NH2
    # Nucleotide methyl (thymine/DT)
    ("H71", "H72"),
    ("H72", "H73"),
    ("H71", "H73"),
    # Nucleotide sugar methylene (deoxyribose H2'/H2'', ribose H5'/H5'')
    ("H2'", "H2''"),
    ("H5'", "H5''"),
    # Methyl (any permutation)
    ("HG21", "HG22"),
    ("HG22", "HG23"),
    ("HG21", "HG23"),
    ("HB1", "HB2"),
    ("HB2", "HB3"),
    ("HB1", "HB3"),
    # Terminal NH3+
    ("H1", "H2"),
    ("H2", "H3"),
    ("H1", "H3"),
    # Methyl on MET etc
    ("HE1", "HE2"),
    ("HE2", "HE3"),
    ("HE1", "HE3"),
]


def _normalize_h_name(name: str) -> str:
    """Normalize H atom name to mmCIF convention.

    PDB legacy uses * for prime (e.g. H1*, H2*1), mmCIF uses ' (e.g. H1', H2'1).
    reduce2 (cctbx) outputs * style; zhyd/ChimeraX output ' style.
    Also handles the digit-after-star pattern: H2*1 -> H2', H2*2 -> H2''.
    """
    if "*" not in name:
        return name
    # H5*1 -> H5', H5*2 -> H5''
    # H1* -> H1'
    if name.endswith("*"):
        return name[:-1] + "'"
    # Pattern like H2*1, H2*2: * followed by digit
    idx = name.index("*")
    suffix_digit = name[idx + 1 :]
    if suffix_digit == "1":
        return name[:idx] + "'"
    elif suffix_digit == "2":
        return name[:idx] + "''"
    else:
        # Fallback: just replace * with '
        return name.replace("*", "'")


def _build_equiv_map() -> dict[str, set[str]]:
    """Build name -> set of equivalent names."""
    m: dict[str, set[str]] = {}
    for a, b in EQUIV_PAIRS:
        m.setdefault(a, {a}).add(b)
        m.setdefault(b, {b}).add(a)
    # Transitive closure
    changed = True
    while changed:
        changed = False
        for k, v in m.items():
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


def extract_h(path: str) -> tuple[dict[tuple, list[AtomPos]], dict[str, int]]:
    """Return ({(auth_chain, auth_seq, comp_id, ins): [AtomPos, ...]}, {skipped_comp: count})

    Uses auth_asym_id + auth_seq_id as the residue key for cross-tool
    compatibility: label_seq_id for non-polymer residues varies between
    tools (e.g. zhyd="1" vs ChimeraX="."), and label_asym_id may
    differ for glycan chains.
    """
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
        name = _normalize_h_name(row[1].strip())
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


def are_equivalent(name_a: str, name_b: str) -> bool:
    """Check if two H atom names are equivalent (swappable)."""
    if name_a == name_b:
        return True
    equivs = EQUIV_MAP.get(name_a)
    return equivs is not None and name_b in equivs


def match_atoms(
    atoms_a: list[AtomPos], atoms_b: list[AtomPos]
) -> tuple[list[tuple[AtomPos, AtomPos, float]], list[AtomPos], list[AtomPos]]:
    """Match H atoms between A and B by minimum distance within each residue.

    Strategy: greedily match each A atom to the nearest unmatched B atom
    (respecting altloc constraints). This handles naming convention differences
    (e.g. HB2/HB3 swapped) and methyl permutations automatically.

    Returns: (matched_pairs, only_a, only_b)
    """
    used_b: set[int] = set()
    matched: list[tuple[AtomPos, AtomPos, float]] = []
    unmatched_a: list[AtomPos] = []

    # Sort candidates by distance to get greedy-optimal matches
    candidates: list[tuple[float, int, int]] = []  # (dist, idx_a, idx_b)
    for ia, aa in enumerate(atoms_a):
        for ib, bb in enumerate(atoms_b):
            if aa.altloc and bb.altloc and aa.altloc != bb.altloc:
                continue
            candidates.append((aa.dist(bb), ia, ib))
    candidates.sort()

    used_a: set[int] = set()
    used_b_set: set[int] = set()
    for d, ia, ib in candidates:
        if ia in used_a or ib in used_b_set:
            continue
        matched.append((atoms_a[ia], atoms_b[ib], d))
        used_a.add(ia)
        used_b_set.add(ib)

    unmatched_a = [a for i, a in enumerate(atoms_a) if i not in used_a]
    unmatched_b = [b for i, b in enumerate(atoms_b) if i not in used_b_set]
    return matched, unmatched_a, unmatched_b


@dataclass
class Stats:
    total_a: int = 0
    total_b: int = 0
    n_matched: int = 0
    n_only_a: int = 0
    n_only_b: int = 0
    dists: list[float] = field(default_factory=list)
    res_diffs: list[tuple[str, int, int, list[str], list[str]]] = field(
        default_factory=list
    )
    skipped: dict[str, int] = field(default_factory=dict)


def compare_all(h_a: dict, h_b: dict) -> Stats:
    stats = Stats()
    stats.total_a = sum(len(v) for v in h_a.values())
    stats.total_b = sum(len(v) for v in h_b.values())

    all_keys = sorted(set(h_a.keys()) | set(h_b.keys()))

    for key in all_keys:
        aa = h_a.get(key, [])
        bb = h_b.get(key, [])

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

        matched, only_a, only_b = match_atoms(aa, bb)
        stats.n_matched += len(matched)
        stats.n_only_a += len(only_a)
        stats.n_only_b += len(only_b)
        stats.dists.extend(d for _, _, d in matched)

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
    print(f"Hydrogen Placement Comparison (equiv-aware, skipped: {skip_label})")
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
        dists = sorted(stats.dists)
        n = len(dists)
        mean_d = sum(dists) / n
        median_d = dists[n // 2]
        p95 = dists[int(n * 0.95)]
        max_d = dists[-1]
        n_close = sum(1 for d in dists if d < 0.1)
        n_mid = sum(1 for d in dists if 0.1 <= d < 0.5)
        n_far = sum(1 for d in dists if 0.5 <= d < 1.0)
        n_bad = sum(1 for d in dists if d >= 1.0)

        print(f"\n  Position differences ({n} matched pairs):")
        print(
            f"    mean={mean_d:.3f}A  median={median_d:.3f}A  P95={p95:.3f}A  max={max_d:.3f}A"
        )
        print(f"    <0.1A: {n_close:5d} ({100 * n_close / n:5.1f}%)")
        print(f"    0.1-0.5A: {n_mid:3d} ({100 * n_mid / n:5.1f}%)")
        print(f"    0.5-1.0A: {n_far:3d} ({100 * n_far / n:5.1f}%)")
        print(f"    >=1.0A: {n_bad:5d} ({100 * n_bad / n:5.1f}%)")

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
        if len(only_a_res) > 20:
            print(f"    ... and {len(only_a_res) - 20} more")

    if only_b_res:
        print(f"\n  Residues with H only in B ({len(only_b_res)}):")
        for r, _, cb, _, ob in only_b_res[:20]:
            print(f"    {r:25s}  B={cb}  [{','.join(ob[:5])}]")
        if len(only_b_res) > 20:
            print(f"    ... and {len(only_b_res) - 20} more")

    if count_diff:
        print(f"\n  Residues with H count differences ({len(count_diff)}):")
        for r, ca, cb, oa, ob in count_diff[:50]:
            extra = ""
            if oa:
                extra += f"  only-A: {','.join(oa)}"
            if ob:
                extra += f"  only-B: {','.join(ob)}"
            print(f"    {r:25s}  A={ca:3d}  B={cb:3d}  diff={ca - cb:+d}{extra}")
        if len(count_diff) > 50:
            print(f"    ... and {len(count_diff) - 50} more")

    if verbose and stats.dists:
        bad = sum(1 for d in stats.dists if d >= 1.0)
        if bad:
            print(f"\n  Pairs with dist >= 1.0A: {bad}")


def main():
    parser = argparse.ArgumentParser(description="Compare H placement (equiv-aware)")
    parser.add_argument("file_a", help="First mmCIF (e.g. zhyd)")
    parser.add_argument("file_b", help="Second mmCIF (e.g. ChimeraX)")
    parser.add_argument("--verbose", "-v", action="store_true")
    parser.add_argument("--summary-only", "-s", action="store_true")
    args = parser.parse_args()

    h_a, skipped_a = extract_h(args.file_a)
    h_b, skipped_b = extract_h(args.file_b)
    stats = compare_all(h_a, h_b)
    # Merge skip counts from both files for reporting
    all_skipped: dict[str, int] = {}
    for comp, n in {**skipped_a, **skipped_b}.items():
        all_skipped[comp] = skipped_a.get(comp, 0) + skipped_b.get(comp, 0)
    stats.skipped = all_skipped
    print_report(stats, args.file_a, args.file_b, args.verbose, args.summary_only)


if __name__ == "__main__":
    main()
