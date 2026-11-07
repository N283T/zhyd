# zhyd empirical pKa engine — specification

Status: draft r2 (2026-11-07); parity mode implemented and accepted with
documented deviations (see [Parity status](#parity-status-r2)); zhyd mode in
progress. Target: an MIT-licensed protein pKa engine in
zreduce that first reproduces the PROPKA 3 empirical model ("parity mode") and
then improves on it ("zhyd mode").

This document describes *what* the model computes: equations, conditions and
parameter values. It deliberately contains no code from PROPKA or zpropka.
Implementers work from this text, the papers below, and black-box oracle
output only (see [Clean-room protocol](#clean-room-protocol)).

## References

- [P3.0] Olsson, Søndergaard, Rostkowski, Jensen. PROPKA3: consistent
  treatment of internal and surface residues in empirical pKa predictions.
  *JCTC* 7, 525–537 (2011). doi:10.1021/ct100578z
- [P3.1] Søndergaard, Olsson, Rostkowski, Jensen. Improved treatment of
  ligands and coupling effects in empirical calculation and rationalization
  of pKa values. *JCTC* 7, 2284–2295 (2011). doi:10.1021/ct200133y
- [P1] Li, Robertson, Jensen. Very fast empirical prediction and
  rationalization of protein pKa values. *Proteins* 61, 704–721 (2005).
  doi:10.1002/prot.20660
- Parameter values: the default parameter set distributed with PROPKA 3.5.1.
  Numerical parameters are treated as facts about the published model;
  they are restated here in tables of our own layout.

## Clean-room protocol

1. The spec author has read PROPKA 3.5.1 and written zpropka (LGPL-2.1).
   Implementers of this engine **must not** open the PROPKA or zpropka
   sources. They get this document, the papers and the oracle's output.
2. The oracle is used only as a black box: `zreduce pka --engine propka`
   (zreduce built with `-Dpropka=true`) and its JSON report.
3. When the engine disagrees with the oracle, the implementer reports the
   site, the terms that differ and the input. The spec author amends this
   document *in words* (new rule, corrected constant, clarified condition),
   adds a row to [Revision log](#revision-log), and the implementer changes
   the code from the amended text.
4. No code, identifiers, comments or file layout are carried over from
   PROPKA/zpropka. Unit-test expectations are derived from the formulas in
   this document (hand calculation), not copied from oracle output; oracle
   comparisons are a separate, data-driven test layer.

This is a process separation to keep the MIT engine an independent
implementation of a published model. It is not legal advice.

### Other sources

Papers and implementations beyond PROPKA may be used freely for ideas and
for zhyd-mode improvements, according to their licenses:

| Source | License | Use |
|---|---|---|
| Papers (any method) | — | ideas, equations, parameters with citation |
| pKAI (Reis et al. 2022), [repo](https://github.com/bayer-science-for-a-better-life/pKAI) | MIT | code may be adapted with attribution |
| PDB2PQR / APBS | BSD-3-Clause (recent releases; check the LICENSE of the version read) | code may be adapted with attribution |
| PyMCCE | MIT | code may be adapted with attribution |
| PypKa | LGPL-3.0 (DelPhi backend proprietary) | papers / black-box oracle only, like PROPKA |
| MCCE4 | not confirmed | read its LICENSE before looking at code |

Copyleft or unlicensed code is treated like PROPKA: papers and black-box
output only. Adapted permissive code keeps its copyright notice in the
source file and is listed in the revision log.

## Scope

**In scope (r1):** proteins, one model at a time.

| Site | Atom(s) defining the site | Titratable |
|---|---|---|
| ASP, GLU carboxyl | CG / CD and its two oxygens | yes |
| C-terminal carboxyl | OXT and the other oxygen on C | yes |
| HIS imidazole | ring CG ND1 CD2 CE1 NE2 | yes |
| CYS thiol | SG (not in a disulfide) | yes |
| TYR phenol | OH | yes |
| LYS amine | NZ | yes |
| ARG guanidinium | CZ and NE, NH1, NH2 | yes |
| N-terminal amine | N of the first residue | yes |
| SER/THR hydroxyl (`ROH`) | OG / OG1 | no (H-bond partner) |
| ASN/GLN amide (`AMD`) | CG / CD, its O and N | no (H-bond partner) |
| TRP indole NH | NE1 | no (H-bond partner) |
| Backbone NH (`BBN`) | N of every non-PRO residue that is not the N-terminus | no |
| Backbone C=O (`BBC`) | C bonded to exactly one O | no |
| Ion | HETATM residue in the [ion table](#ions) | no (fixed charge) |

**Out of scope (r1):** ligand groups (PROPKA's Sybyl-type ligand typing),
nucleic acids, isoelectric point and folding free energy profiles. Covalent
coupling between titratable groups does occur in proteins (r2) and is
specified in [Covalent coupling](#covalent-coupling).

## Pipeline

For each model:

1. Build conformations ([Conformations](#conformations)).
2. Per conformation: detect sites and partner groups; add internal hydrogens
   ([Internal hydrogens](#internal-hydrogens)).
3. Desolvation for every titratable group and every ion
   ([Desolvation](#desolvation)).
4. Backbone hydrogen bonds ([Backbone](#backbone-hydrogen-bonds)).
5. Ion Coulomb terms ([Ions](#ions)).
6. Local backbone reorganisation term for ASP/GLU
   ([Reorganisation](#local-backbone-reorganisation-term)).
7. Side-chain pair interactions, non-iterative pairs first, then the
   iterative pairs ([Pair interactions](#side-chain-pair-interactions),
   [Iterative scheme](#iterative-scheme)).
8. Total pKa ([Total](#total-pka)).
9. Average over conformations.

The reported terms map onto the zhyd schema as: `burial` = desolvation +
reorganisation; `hydrogen_bond` = side-chain + backbone H-bond determinants;
`electrostatic` = Coulomb determinants (pair and ion); `coupling` = 0.

## Input preparation

- Heavy atoms only from the input; existing hydrogens are discarded.
- Residues named HOH, H2O, SO4, PO4, PEG, EPE, TRS are removed before
  anything else (parity). zhyd mode additionally removes every residue zreduce
  classifies as water.
- **Disulfide:** two CYS SG atoms closer than 2.5 Å. Both CYS are not
  titratable and take no part in pair interactions (they still count as
  atoms for desolvation of others). Their reported value is "not titrating";
  zreduce lists them in `excluded_sites`.
- **N-terminus (parity):** the N atom of the first ATOM-record residue of each
  chain (after a chain start, model start or TER). Chain breaks inside a chain
  do not create an N-terminus. HETATM residues never form one.
- **C-terminus (parity):** only a residue that has an OXT atom. Without OXT
  there is no C-terminal site.
- Bonds (r2): two heavy atoms are bonded when d ≤ 2.0 Å (S–S: d < 2.5 Å),
  regardless of residue topology; distorted geometry therefore changes the
  bond graph. Bond graphs are used for (a) the 4-bond exclusion of side-chain
  H-bonds, (b) placing internal hydrogens, (c) carboxyl and backbone C=O
  atoms (see [Group geometry](#group-geometry)) and (d) covalent coupling.
- **Parity input (r2):** what zreduce's glue passes to the oracle: heavy
  atoms only, coordinates rounded to 3 decimals, only standard polymer
  residues as ATOM records, zreduce waters dropped, chains identified by the
  auth chain id.
- **N-terminus details (r2):** chain-start detection also restarts after a
  residue carrying OXT. Every residue sharing the residue number of the
  chain's N-terminal residue (insertion-code variants) is an N-terminus too.
  UNK residues receive termini and backbone groups like standard residues.

## Conformations

Parity mode follows alternate locations: each altloc letter defines one
conformation containing the blank-altloc atoms plus the atoms with that
letter (a blank-only structure has exactly one conformation). All steps run
independently per conformation. The reported pKa and every reported term of
a site is the arithmetic mean over conformations; sites are identified by
(chain, residue number, insertion code, site type).

Refinements (r2), each inferred from oracle reports:

- **Uneven alternates:** an atom with no alternate for letter L uses its
  first listed alternate in conformation L.
- **Not titrating:** a site that does not titrate in a conformation (for
  example through covalent coupling) enters the mean as pKa 99.99 with all
  terms 0. A site that does not titrate in the first conformation is reported
  as not titrating.
- **Insertion codes:** conformations are matched by (chain, residue number,
  site type), ignoring the insertion code. In every conformation after the
  first, only the last residue of each (chain, number) family is present. The
  mean runs over all conformations; a site absent from one counts as 0 there.
- **Shared identities:** when one conformation holds several groups with the
  same identity, the first in input order supplies the value, and every
  residue carrying that identity reports it.

## Internal hydrogens

The engine adds its own hydrogens for the energy terms only; they are never
written out and are independent of zreduce's placement and bond-length
policy. Hydrogens are needed on:

| Atom | Count | Geometry |
|---|---|---|
| backbone N (BBN) | 1 | bisector |
| HIS ND1 and NE2 | 1 each (charged form, both N protonated) | bisector |
| ARG NE | 1 | bisector |
| ARG NH1, NH2 | 2 each | planar pair |
| ASN ND2 / GLN NE2 | 2 | planar pair |
| TRP NE1 | 1 | bisector |

N–H length is 1.01 Å.

- **Bisector:** for N with heavy neighbours A and B, the H lies along
  −(â + b̂), where â and b̂ are unit vectors from N towards A and B. A and B
  are the first two bonded heavy neighbours of N in input order (r2; matters
  when mixed alternates give N a third bond). A backbone N gets an H only when
  it has exactly two bonded heavy neighbours.
- **Planar pair:** for N bonded to one trigonal atom X with other neighbours
  Y₁ and Y₂, let u = unit(X − N) and p = the unit component of (Y₁ − Y₂)
  perpendicular to u. The two H are at N + 1.01·(−cos 60°·u ± sin 60°·p),
  i.e. 120° from N→X (r2; this replaces the less precise "plane of N, X and
  X's neighbours" wording, which fails for non-coplanar input).

## Group geometry

Each group has a **centre** and two **interaction-atom lists**: one used when
the partner is an acid or neutral group ("acid list") and one used when the
partner is a base ("base list"). Bases are ARG, LYS, HIS and the N-terminus.

| Group (type) | Centre | Acid list | Base list |
|---|---|---|---|
| ASP/GLU (`COO`) | mean of the two O | 2 O | 2 O |
| C-terminus (`COO`) | mean of O, OXT | O, OXT | O, OXT |
| HIS (`HIS`) | mean of the 5 ring atoms | ND1, NE2, HD1, HE2 | ND1, NE2 |
| CYS (`CYS`) | SG | SG | SG |
| TYR (`TYR`) | OH | OH | OH |
| LYS (`LYS`) | NZ | NZ | NZ |
| ARG (`ARG`) | CZ | NE, NH1, NH2 + 5 H | NE, NH1, NH2 |
| N-terminus (`N+`) | N | N | N |
| SER/THR (`ROH`) | OG / OG1 | same atom | same atom |
| ASN/GLN (`AMD`) | mean of O and N | N + 2 H | O |
| TRP (`TRP`) | NE1 | NE1, HE1 | NE1 |
| Backbone NH (`BBN`) | N | N, H | N, H |
| Backbone CO (`BBC`) | C | O | O |

The *group atom* (used for the residue-exclusion and bond-distance rules) is
CG (ASP), CD (GLU), OXT (C-terminus), CG (HIS), SG, OH, NZ, CZ (ARG), N
(N-terminus), OG/OG1, CG/CD (amides), NE1, backbone N, backbone C.

Carboxyl and C=O atoms come from the bond graph (r2): an ASP/GLU carboxyl
uses only the oxygens bonded to its group carbon. If none is bonded, the
centre is the group atom and the group has no interaction atoms. A backbone
C=O exists only for a C with exactly one bonded O.

Distances "between groups" are between centres. "Closest pair" means the
minimum distance over the cross product of the two interaction-atom lists.

## Model values and charges

| Site | Model pKa | Charge q |
|---|---|---|
| C-terminus | 3.20 | −1 |
| ASP | 3.80 | −1 |
| GLU | 4.50 | −1 |
| HIS | 6.50 | +1 |
| CYS | 9.00 | −1 |
| TYR | 10.00 | −1 |
| LYS | 10.50 | +1 |
| ARG | 12.50 | +1 |
| N-terminus | 8.00 | +1 |
| ROH, AMD, TRP, BBN, BBC | — | 0 |

q is the charge of the *charged* form, used as the sign of every shift.

## Desolvation

For group G with centre c (titratable groups and ions):

- Atom set: all heavy atoms of the conformation (protein and HETATM) except
  atoms whose (residue number, chain) equals that of G's group atom.
  *Parity quirk:* the insertion code and residue name are not compared.
- Relative volume v(a): carbon named CA or C → 1.40; any other carbon → 2.64;
  N 1.06; O 1.00; S 1.66; F 0.90; Cl 1.53; P 1.66; any other element 1.00.
- With d = |a − c|:
  - V = Σ v(a) / max(2.75⁴, d⁴) over atoms with d < 20 Å.
  - N = number of atoms with d < 15 Å.
- Burial weight w = clamp((N − 280) / (560 − 280), 0, 1). This is the
  reported `burial_fraction`.
- Surface scaling s = 1 − (1 − 0.25)(1 − w).
- Desolvation shift ΔpKa_desolv = q · (−13.0) · max(0, V − 0.0) · s.

Pair weight for two groups with counts N₁, N₂:
w₁₂ = clamp((N₁ + N₂ − 560) / (1120 − 560), 0, 1).

## Distance and angle factors

H-bond profile for distance d with cutoffs (c₁, c₂) and maximum Δ:

- f(d) = 1 if d < c₁; 0 if d > c₂; otherwise 1 − (d − c₁)/(c₂ − c₁).
- E_hb = |Δ · f(d) · F_angle|.

Angle factor F_angle for donor heavy atom D, its hydrogen H and acceptor A:
the cosine between the D→H direction and the H→A direction
(1 for a linear D–H···A). *Parity quirk:* the absolute value above means a
negative cosine still yields a positive magnitude.

Coulomb profile for distance r and pair weight w:

- ε = 160 − (160 − 30) · w
- r′ = max(r, 4.0)
- t = clamp((r′ − 10.0) / (4.0 − 10.0), 0, 1)
- E_coul = |244.12 / (ε · r′) · t|

## Backbone hydrogen bonds

For every titratable group T (its **acid list** is always used here) and
every backbone group B:

1. Closest pair (b, t) between B's list and T's acid list, distance d.
2. Parameters by (B type, T type); no parameters → no interaction:

   | Backbone | Site type | Δ | c₁ | c₂ |
   |---|---|---|---|---|
   | BBN (N–H donor) | COO (ASP, GLU, C-terminus) | 0.85 | 2.00 | 3.00 |
   | BBN | CYS | 0.85 | 3.00 | 4.00 |
   | BBN | TYR | 0.85 | 2.20 | 3.20 |
   | BBC (C=O acceptor) | HIS | 0.85 | 2.00 | 3.00 |

3. Only if d < c₂:
   - BBN: if b is the backbone H, F_angle = angle factor (N, H, t);
     if b is the N, F_angle = 0.
   - BBC with T of type HIS: if t is a hydrogen, F_angle = angle factor
     (t's N, t, b's O); if t is a nitrogen, F_angle = 0.
   - Otherwise F_angle = 1.
4. If F_angle > 0.001, add backbone determinant q_T · E_hb(d, Δ, c₁, c₂,
   F_angle) to T. Note that d is the closest-pair distance, not
   recomputed from the hydrogen.

## Ions

| Residue | Charge | Residue | Charge |
|---|---|---|---|
| MG | +2 | CD | +2 |
| CA | +2 | FE | +3 |
| ZN | +2 | SR | +2 |
| NA | +1 | CU | +2 |
| CL | −1 | IOD | −1 |
| MN | +2 | HG | +2 |
| K | +1 | BR | −1 |
| CO | +2 | NI | +2 |
| FE2 | +2 | 1P, 2P, 1N, 2N | +1, +2, −1, −2 |

For every titratable group T and every ion I with centre distance r < 10 Å:
Coulomb determinant on T of −q_I · E_coul(r, w_TI), where w_TI is the pair
weight. No burial threshold applies to ions.

## Local backbone reorganisation term

For ASP and GLU (not C-terminus) only. For each backbone C=O group with O at
distance d from T's centre (d = |centre − O|):

- cosθ = cosine between the C→O direction and the O→centre direction.
- If d < 6.0 and cosθ > 0.001: add 0.8 · min(1, 1 − (d − 3.0)/(6.0 − 3.0)).

Sum × w_T (T's burial weight) is the reorganisation shift. It is reported
inside `burial`.

## Side-chain pair interactions

Candidate pairs: every unordered pair of side-chain groups (titratable sites,
ROH, AMD, TRP; not backbone groups, not disulfide CYS) with centre distance
r < 10 Å. Groups are ordered by their group atom's position in the input; in
a pair (G₁, G₂), G₁ is the *later* group. Non-titratable groups have no
desolvation step, so their atom count N is 0 in every pair-weight and
buried-pair formula. The interaction class comes from the matrix (symmetric;
I = iterative, N = non-iterative, – = none):

|     | CYS | N+ | HIS | LYS | AMD | COO | ARG | TRP | ROH | TYR |
|---|---|---|---|---|---|---|---|---|---|---|
| CYS | I | N | I | N | N | I | N | N | N | N |
| N+  | N | I | N | N | – | N | N | – | – | I |
| HIS | I | N | I | N | N | I | N | – | – | I |
| LYS | N | N | N | I | – | N | N | – | – | I |
| AMD | N | – | N | – | – | N | – | – | – | N |
| COO | I | N | I | N | N | I | N | N | N | N |
| ARG | N | N | N | N | – | N | I | – | – | N |
| TRP | N | – | – | – | – | N | – | – | – | N |
| ROH | N | – | – | – | – | N | – | – | – | N |
| TYR | N | I | I | I | N | N | N | N | N | I |

### Side-chain H-bond magnitude h(G₁, G₂)

1. Interaction atoms: each group uses its base list if the *other* group is a
   base, else its acid list. Closest pair (a₁, a₂), distance d.
2. Cutoffs (c₁, c₂) by the two group types (symmetric; default 3.0/4.0);
   Δ = 0.85 for every pair:

   | Pair | c₁ | c₂ | Pair | c₁ | c₂ |
   |---|---|---|---|---|---|
   | COO–COO | 2.50 | 3.50 | HIS–AMD | 2.00 | 3.00 |
   | COO–ARG | 1.85 | 2.85 | HIS–TYR | 2.00 | 3.00 |
   | COO–LYS | 2.85 | 3.85 | CYS–CYS | 3.00 | 5.00 |
   | COO–HIS | 2.00 | 3.00 | CYS–TRP | 2.50 | 3.50 |
   | COO–AMD | 2.00 | 3.00 | CYS–ROH | 3.50 | 4.50 |
   | COO–TRP | 2.00 | 3.00 | CYS–AMD | 2.50 | 3.50 |
   | COO–ROH | 2.65 | 3.65 | CYS–TYR | 3.50 | 4.50 |
   | COO–TYR | 2.65 | 3.65 | CYS–N+ | 3.00 | 4.50 |
   | COO–N+ | 2.85 | 3.85 | TYR–TYR | 3.50 | 4.50 |
   | ARG–CYS | 2.50 | 4.00 | TYR–N+ | 3.00 | 4.50 |
   | ARG–TYR | 2.50 | 4.00 | TYR–AMD | 2.50 | 3.50 |
   |  |  |  | TYR–TRP | 2.50 | 3.50 |
   |  |  |  | TYR–ROH | 3.50 | 4.50 |

3. No interaction if d ≥ c₂.
4. No interaction if the two group atoms are within 4 bonds of each other.
5. Angle: if G₂ is angular (HIS, ARG, AMD, TRP) and a₂ is a hydrogen,
   F_angle = angle factor (a₂'s N, a₂, a₁) and d = |a₁ − a₂|; if a₂ is a
   heavy atom, F_angle = 0. Otherwise the same test on G₁. Otherwise 1.
6. Special cases, replacing the general value:
   - **COO–ARG:** sum over the two closest *disjoint* atom pairs (take the
     closest pair, remove both atoms, take the next closest). Each pair uses
     the COO–ARG cutoffs and an angle factor computed from (ARG atom's first
     heavy neighbour → ARG atom, ARG atom → COO atom). For an ARG H the first
     heavy neighbour is its N. For an ARG N, the bonded heavy neighbour that
     gives the largest |cos| is used (r2: NE can use CD or CZ; neither fixed
     choice reproduces the oracle).
   - **COO–COO:** F_angle = 1, value × (1 + w₁₂).
   - **COO–HIS, CYS–HIS:** 1.60 if the pair is *buried*, else general value.
   - **CYS–CYS:** 3.60 if buried, else general value.
   - *Buried pair:* not (N₁ + N₂ ≤ 900 and (N₁ ≤ 400 or N₂ ≤ 400)).
7. Otherwise h = E_hb(d, 0.85, c₁, c₂, F_angle).

### Coulomb magnitude c(G₁, G₂)

Only when both groups are titratable, r ≤ 10 Å and N₁ + N₂ ≥ 280:
c = E_coul(r, w₁₂). Otherwise 0.

### Non-iterative pairs (N)

H-bond part (when h > 0):
- equal charges q₁ = q₂: the group with the lower model pKa gets −h, the
  other +h;
- otherwise each group i gets q_i · h (a neutral partner receives 0).

Coulomb part (when c > 0):
- both acids: the group with the higher model pKa gets +c;
- both bases: the group with the lower model pKa gets −c;
- acid–base: each group i gets q_i · c.

H-bond parts are side-chain determinants; Coulomb parts are Coulomb
determinants. Ties in model pKa: for the H-bond part G₁ gets +h and G₂ −h;
for the Coulomb part the shift goes to G₂.

### Iterative scheme

Applies to the I pairs with h > 0 or c > 0. Each group taking part starts
from its non-iterative value

p⁰ = model + desolvation + reorganisation + all determinants assigned so far
(backbone, ion, non-iterative pairs).

Each pair keeps a persistent correction (a₁, a₂), initially (0, 0). One
iteration, using p from the previous iteration (pᵒˡᵈ, initially p⁰):

- **Both acids:** with D = c + 2h, compare pᵒˡᵈ₁ + a₁ against pᵒˡᵈ₂ + a₂
  (strictly greater wins; ties go to group 2). Reset (a₁, a₂) = (0, 0). The
  winner w gets +h (side-chain) and +c (Coulomb), the loser −h (side-chain),
  and a_w = −D.
- **Both bases:** with D = −(c + 2h), compare pᵒˡᵈ₁ + a₁ + D against
  pᵒˡᵈ₂ + a₂ + D (strictly smaller wins; ties go to group 2). Reset (a₁, a₂).
  The winner gets −h and −c, the loser +h, and a_w = −D.
- **Otherwise** (acid–base or a neutral partner): kᵢ = pᵒˡᵈᵢ + aᵢ +
  qᵢ(c + h). Reset (a₁, a₂). Apply when (q₁ = −1 and k₁ < k₂) or
  (q₁ = +1 and k₁ > k₂). When applied: if c > 0.005, each group i gets
  qᵢ·c (Coulomb) and aᵢ −= qᵢ·c; if h > 0.005, each group i gets qᵢ·h
  (side-chain) and aᵢ −= qᵢ·h.

pⁿᵉʷ = p⁰ + the contributions of this iteration. Stop when pⁿᵉʷ = pᵒˡᵈ for
every group (exact equality), or after 10 iterations. The contributions of
the final iteration whose magnitude exceeds 0.005 become the groups'
determinants. Pair order: groups in input atom order; pair (i, j) with
j < i, iterating i ascending then j ascending.

## Total pKa

pKa = model + desolvation + reorganisation + Σ side-chain + Σ backbone +
Σ Coulomb determinants.

## Covalent coupling

(r2) Two titratable groups whose group atoms are within three bonds of each
other are covalently coupled. One of them stops titrating and drops out of all
pair interactions; it still counts for desolvation of others.

- Acid with base: the acid drops out.
- Two bases: the one with the lower pKa so far drops out.
- Two acids: the one with the higher pKa so far drops out.
- Two N-termini never couple.

In proteins this removes ASP/HIS side chains on N-terminal residues and
groups in clashing NMR models.

## Non-covalent coupling flag

PROPKA marks pairs of titratable sites whose protonation states are coupled
(strong mutual interaction, similar intrinsic pKa, small free-energy change
on swapping). The flag does not change the reported pKa. r1 reports
`coupled = false`; the detection rule will be specified in a later revision.

## Known parity quirks

These are reproduced in parity mode and corrected in zhyd mode:

| # | Parity behaviour | zhyd mode |
|---|---|---|
| Q1 | Same-residue exclusion ignores insertion code | compare full residue identity |
| Q2 | Negative angle cosine counts as a positive H-bond | F_angle = max(0, cos) |
| Q3 | No C-terminus without OXT | chain-last polymer residue is a C-terminus |
| Q4 | N-terminus only after TER / chain start; chain breaks ignored | follows zreduce chain-break detection (`--nterm` policy) |
| Q5 | Ions by residue name only | metals by element and CCD; coordination-aware exclusion (already in zreduce) |
| Q6 | Crystal waters count as burial atoms (if not removed) | waters removed |
| Q7 | Averaging 99.99 for conformations where a site does not titrate yields values like 52.96 | average over titrating conformations only; values outside the physical range flagged `pka_outside_physical_range` |

## Implementation layout

- `src/pka_empirical.zig` (public entry) with submodules under
  `src/pka_empirical/`: sites, hydrogens, desolvation, interactions,
  iterative, engine.
- Engine option `mode: parity | zhyd`; parameter tables as one versioned
  struct (`zhyd-empirical-v1-parity`, later `zhyd-empirical-v1`).
- CLI: `zreduce pka --engine empirical`. When it meets the acceptance gates it
  replaces the current `builtin` (v0.2) engine.

## Validation

Tooling: `scripts/pka/compare_engines.py` joins two zhyd pKa reports (oracle
and candidate) by site identity and reports per-site differences in the
predicted value and each term.

Gates, in order:

1. **Unit tests** per section, with expectations computed by hand from this
   document (single atoms, two-group fixtures, a three-group iterative case).
2. **Parity (G1):** on the 185-structure PKAD-R set and the zreduce test
   structures, parity mode vs `--engine propka` (waters dropped, which matches
   the parity water rule) on the same input file: |Δ predicted| ≤ 0.02 for ≥ 99 % of sites, ≤ 0.1 for all
   sites, and the same set of reported sites. Every exceedance is either
   fixed or recorded here as a deviation. Parity is a check that the
   implementation is faithful, not a measure of accuracy: PROPKA's own error
   against experiment is about 0.7 pK units. Residual deviations well below
   that, explained or not, do not block zhyd mode once they are recorded.
3. **Accuracy (G2):** zhyd mode on the cluster split (`evaluate.py`): test MAE
   no worse than PROPKA's and below the null model; SNase stress set
   reported separately. Each zhyd-mode correction (Q1–Q7, recalibration) is
   evaluated on its own.

### Parity status (r2)

Parity mode against the oracle (`compare_engines.py`, same input file):

| Set | Files | Sites | ≤ 0.02 | > 0.1 | Max | Same site set |
|---|---|---|---|---|---|---|
| PKAD-R, protein only (non-ion HETATM and nucleotides stripped) | 162 | 19,598 | 99.84 % | 19 | 0.825 | yes |
| same without 4HHB | 161 | 19,424 | 99.98 % | 2 | 0.177 | yes |
| PKAD-R as deposited (ligands present) | 162 | 19,502 | 95.01 % | 813 | 6.66 | yes |
| zreduce test/example structures, protein only | 8 | 1,188 | 98.74 % | 5 | 0.84 | yes |

The model term matches on every site. The "≤ 0.1 for all sites" bound fails
on 19 sites (17 in 4HHB); G1 is accepted with these recorded as deviations,
and further work is judged by G2. Unexplained deviations:

- **D1** 4HHB: badly distorted geometry (HIS ring bond 0.67 Å); no choice of
  ring atoms reproduces the oracle's HIS centre.
- **D2** 1LNI B14 GLU / B17 ASP (±0.18): the oracle's acid–acid winner flips
  only when two third-party charges are both present; the iterative scheme as
  written cannot produce this.
- **D3** 3EBX A6 HIS (0.08): an ASN H-bond is absent in the oracle's second
  conformation.
- **D4** Backbone NH after a chain break where N has one bonded atom: the
  oracle builds an H (109.5° from N→CA, orientation rule unknown); parity
  mode skips it.
- **D5** 1MEK models 8 and 29: the acid–acid coupling choice is opposite to
  [Covalent coupling](#covalent-coupling).
- **D6** Insertion-code families (2HNT, 2CF8, 6FYS): backbone H-bond
  residuals up to 0.84.
- **D7** Ligand and nucleic-acid groups: out of scope (943 of the 974
  deposited-set outliers above 0.02 disappear when they are stripped).

## Revision log

| Rev | Date | Change |
|---|---|---|
| r1 | 2026-11-06 | Initial protein-only specification |
| r1.1 | 2026-11-06 | Allow other papers and permissively licensed implementations |
| r2 | 2026-11-07 | Amendments from the first parity implementation (black-box oracle probes): distance-based bonds and carboxyl atoms, covalent coupling, conformation refinements (uneven alternates, not-titrating averaging, insertion codes, shared identities), planar-pair and bisector H rules, COO–ARG neighbour rule, N-terminus details; Q7 reworded; parity status and deviations D1–D7 |
