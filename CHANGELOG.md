# Changelog

All notable changes to zhyd will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/),
and this project adheres to [Semantic Versioning](https://semver.org/).

## [Unreleased]

### Changed
- **Atom-pair scores are reduce's dot score, integrated in closed form**
  (`scorer.pairDotScore`). The old pair score judged the whole pair by its
  center distance: an H-bond whose center overlap passed the 0.6 A allowance
  (H...O below about 1.85 A, the usual OH H-bond) became a full bump of -3.5.
  reduce scores dots on the H surface and only penalizes the overlap beyond
  the allowance. Scores are now in reduce's units (square Angstroms), which
  the flip penalties already assumed. On the 100-structure set H-bonds go
  from 143,750 to 152,174 (reduce -NUClear: 151,794); Ser OG H-bonds go from
  2,203 to 4,783 (4,828) and Thr OG1 from 3,021 to 5,515 (5,034)
- Atoms within three bonds of a mover atom are no longer scored against it,
  as in reduce. Bonds are inferred from covalent radii. The 1-4 terms changed
  with the torsion and pushed OH groups away from orientations eclipsing CA
- Two charged atoms only H-bond when their charges are opposite, as in reduce
- **Methyl groups are no longer optimized.** They stay staggered as placed.
  The score has no torsion term, so optimization turned 21 % of methyls
  eclipsed (24,425 of 116,224 on the 100-structure benchmark set), and their
  movers (about 65 % of all movers) joined most of a large protein into one
  clique that only got the greedy search. `methyl_rotator` and
  `aromatic_methyl` mover kinds are gone; fix-file entries for methyl groups
  now match no mover
- The interaction graph joins two movers only when some atom of one, in any
  orientation, is within scoring range (radius sum + 0.5 A) of some atom of
  the other, instead of within 6 A at the current positions. Components are
  therefore independent and are optimized exactly; the brute-force limit is
  raised from 100,000 to 1,000,000 states. On the 100-structure set, greedy
  cliques drop to 19 and 9qr3 goes from one 8,933-mover clique to 588
  brute-force cliques and 1,706 singletons
- NH3+ rotators use 0/30/60/90 degree orientations (fine search +-15 degrees
  in 5 degree steps). The old 0/+60/-60 set had two identical orientations
  under the group's 3-fold symmetry. NH3+ fix-file indices change meaning
- **Renamed from zreduce to zhyd** (repository, executable, Zig module and
  package). The scope now covers protonation states and pKa, not only
  reduce-style hydrogen placement. No `zreduce` alias is kept. The mmCIF
  output's log category is now `_zhyd_log` (was `_zreduce_log`) and its
  data block `data_ZHYD`; `scripts/pka` options `--zreduce` are now `--zhyd`

### Added
- `evaluate.py` reports error by experimental-shift bin, a fixed reference
  null model, protonation-state accuracy at pH 7 and how many state flips
  (sites whose state differs from the reference pKa) are predicted
- `CITATION.cff` and `THIRD_PARTY_NOTICES.md`
- `zhyd pka <dir>` batch mode: one `<stem>.json` per structure file,
  byte-identical to single-file runs, in parallel (`-j N`; PROPKA engine on one
  thread), with a per-file `--jsonl` log; failed files do not stop the batch
- `--timing` for `run`, `batch` and `pka`: per-stage timing breakdown (read,
  parse, chemistry, place, optimize, validate, pKa, write), wall time, and file
  count/throughput in batch modes, on stderr
- Generic `ThreadPool` (`src/thread_pool.zig`, adapted from zsasa); `batch`
  now runs its workers on it
- `zhyd-empirical-v1.2`: CYS, TYR and HIS parameters refitted on PKAD-R merged
  with the external set (`scripts/pka/merge_sets.py`,
  `calibrate_empirical.py --rows all --nested-folds`), kept by nested
  cluster CV; fit record `scripts/pka/params/zhyd-empirical-v1.2.json`
- `scripts/pka/external_set.py` builds an external pKa test set from PKAD-3_V2
  and the KaML time-stamped set with no PKAD-R structure or 30 % sequence
  homolog; results for all engines are in the empirical-engine spec
- `zhyd pka` without PROPKA now defaults to the empirical engine in zhyd
  mode (`zhyd-empirical-v1.2`); `--engine builtin` keeps the v0.2 model,
  which `--model-only` selects implicitly
- Experimental, independently implemented staged protein pKa predictor for
  ASP/GLU/HIS/LYS/CYS/TYR and peptide termini (`zhyd pka`), with explicit
  model, burial, hydrogen-bond, electrostatic, and mean-field coupling terms
- Versioned `zhyd.protonation.v1` JSON schema and `run --state-json` sidecar for
  applied residue states, formal charges, and provenance
  - Disulfide-bonded cysteines are excluded from pKa sites and reported as
    `disulfide_bonded` (source `covalent_bond`) in the state report
  - Disulfide and metal-coordinated residues are listed in `excluded_sites`
    with a reason; `run --state-json` warns on metal ligands
  - Crystal waters no longer count toward burial or hydrogen-bond terms
  - Parameter set `zhyd-pka-protein-v0.2`, fitted on the PKAD-R training split
    (`scripts/pka/params/`); `pka --parameters PATH` loads a custom set
  - PKAD-R normalization, structure fetch, prediction join, and calibration
    scripts (`scripts/pka/pkadr.py`, `scripts/pka/calibrate.py`)
- Optional PROPKA 3.5.1 engine for `zhyd pka` via the LGPL-2.1 zpropka
  package (`zig build -Dpropka=true`, `--engine propka|builtin`,
  `--propka-keep-waters`); output matches `propka3` on 185 PKAD-R structures
  (14,291 sites, max |Δ| 0.005 = PROPKA's print rounding). The report records
  the engine; ARG sites are reported by this engine
- CYS bonded to a non-cysteine heavy atom is excluded as `covalently_bonded`;
  metal ligands are classified before disulfides, so `_struct_conn` metalc
  flags no longer make a metal-bound CYS look like a cystine
- Values outside -5..20 are kept but flagged `pka_outside_physical_range`
- External PROPKA oracle and cluster-split PKAD-R evaluation harness under
  `scripts/pka/`
- Specification for an MIT empirical pKa engine that reproduces the PROPKA 3
  protein model from a written description (clean-room protocol, parity
  quirks, acceptance gates): `docs/specs/zhyd-pka-empirical-engine.md`;
  `scripts/pka/compare_engines.py` compares two zhyd pKa reports site by site
- Experimental MIT empirical pKa engine, parity mode (`zhyd pka --engine
  empirical`, `src/pka_empirical/`), implemented from that specification
  with PROPKA used only as a black-box oracle; 99.84 % of 19,598 protein
  sites on 185 PKAD-R structures within 0.02 of PROPKA (spec r2 records the
  inferred rules and remaining deviations)
- zhyd mode of the empirical engine (`--empirical-mode zhyd`): corrections
  Q1–Q8 of the parity quirks, each behind its own switch
  (`--empirical-fixes Q1,..|none` for ablation), all adopted (spec r4), and
  a recalibrated table `zhyd-empirical-v1.1` (CYS and TYR model pKa and
  desolvation scale, fitted on the PKAD-R train split with
  leave-one-cluster-out CV; test MAE 1.170 → 0.712). Parity remains the
  default and its output is unchanged
- `zhyd pka --engine empirical --empirical-mode zhyd --parameters PATH`:
  per-site model pKa and desolvation scale overrides (JSON with its own
  version; unknown keys rejected)
- `scripts/pka/calibrate_empirical.py`: grid fit with leave-one-cluster-out
  CV for one site kind of the empirical engine
- `pkadr.py predict --extra-arg ARG` (repeatable) passes extra options to
  `zhyd pka`
- SDF/MOL V2000 parser (`src/sdf.zig`) for providing ligand topology for compounds not in the CCD (#257)
  - `-s, --sdf PATH` flag on `run` and `batch` subcommands
  - Multiple molecules per file (separated by `$$$$`), molecule name used as comp_id
  - `M  CHG` charge parsing, V2000 atom/bond blocks
- Distance-based bond inference fallback (`src/place/distance_derive.zig`) for residues with no dictionary entry (#257)
  - Automatic last-resort when no CCD, SDF, or inline dict entry exists
  - Bond detection via covalent radii (1.3× tolerance, 2.5 Å cap)
  - Valence-rule bond order promotion (C=O > C=N > C=C priority)
  - Synthetic hydrogen generation for `ccd_derive.derivePlans()` consumption
  - Diagnostic output: `note: N residues used distance-based bond inference`
- Placement fallback chain: hardcoded → inline dict → SDF dict → CCD → distance-based inference (#257)
- `--nterm MODE` flag on `run` and `batch` to control N-terminal protonation (#251)
  - `auto` (default): NH3+/NH2+ on chain-first residues only, break-amide on gaps (matches ChimeraX addh, unchanged from prior behavior)
  - `aggressive`: NH3+/NH2+ on both chain-first and chain-break residues (matches reduce2 `first_in_chain`)
  - `neutral`: non-PRO real N-termini get a neutral NH2 (H2/H3 at gauche ±60°, no positive charge flag). PRO is an exception: it keeps the NH2+ geometry and positive charge flag because its secondary amine is protonated at physiological pH.
- `placeNtermNH2Neutral` in `src/place/terminal.zig` for neutral primary amine placement
- `NtermMode` enum re-exported from `zhyd.place` for library callers, with `treatsBreakAsNterm` and `placesPositiveCharge` predicate methods consolidating mode-dependent dispatch

### Fixed
- Fix-file entries that match no mover (a typo, a methyl group, a residue
  missing from the model) now log a warning instead of being ignored silently
- The run summary no longer calls H skipped for a missing reference atom a
  "potential plan bug"; truncated side chains are the usual cause, and reduce
  skips the same H
- Optimization was slow on structures with large cliques: brute force rescored
  every clique mover for every orientation combination, and every mover score
  scanned all movers in the model. Brute force now tabulates per-mover and
  per-pair scores once, and each mover keeps a precomputed list of movers it
  can interact with. Output is unchanged (byte-identical on the 100-structure
  benchmark set except one Lys NH3+ whose three H labels permute); 8edu goes
  from 54 s to 0.3 s and the whole set from 456 s to 13 s of CPU time
- Absent His H atoms parked at `ABSENT_H_POS` were scored against each other
  and counted in mover bounding radii and the interaction graph; they are now
  skipped. The bump terms they added were the same for every candidate, so
  chosen orientations do not change
- `--json` logs for multi-model runs could contain NUL bytes in chain and
  residue names: `MoverSnapshot.chainLabel`/`compId` returned slices into a
  by-value copy of the snapshot
- `batch --jsonl` kept only the last few lines: each line was written at file
  offset 0 and overwrote the previous ones; lines are now appended
- `_entity` written as key-value pairs (single-entity files) now sets residue
  entity types; such files previously had no polymer residues and no pKa
  sites (#264)
- Protonation override files now use a bounded whole-file read, avoiding Zig
  0.16 `EndOfStream` failures on ordinary short control files
- N-terminal PRO residues now actually receive NH2+ hydrogens (H2, H3) during hydrogen placement. Previously the dispatch required a backbone amide plan in the residue's standard plan list, but PRO has no such plan (its N is already bonded to CA and CD), so `placeNtermNH2Pro` was effectively dead code and N-terminal PRO received no backbone H at all. Uncovered during review of #252.
- N-terminal NH2 groups are now optimized as rigid-body pairs instead of as two independent single-H rotators (#253). A new `rotate_nh2` mover hint and `nh2_rotator` mover kind bundle H2/H3 so both atoms rotate together around the N-CA axis. Previously `placeNtermNH2Pro` and `placeNtermNH2Neutral` tagged each H with `.rotate`, causing `mover_gen` to emit two independent single-H rotators that sampled the same axis independently and destroyed the HNH angle set by the placer. The fix applies to both N-terminal PRO (all modes) and non-PRO neutral N-termini.
- `NtermResult` now carries a `missing_ref` counter that the placer tallies into `PlacementResult.n_skipped_missing_ref` (#254). Previously `placeNtermNH3`, `placeNtermNH2Pro`, `placeNtermNH2Neutral`, and `place3primeOH` silently returned `{ 0, 0 }` when required backbone atoms were missing at the target altloc, so the terminal-H shortfall was invisible in the run summary. The existing "N H skipped due to missing reference atoms" warning now covers the N-terminal path too.

## [0.3.0] - 2026-03-30

### Added

#### Performance Optimization
- Mover centroid early-exit: skip distant mover-vs-mover scoring using precomputed bounding spheres (#46)
- Multithreaded singleton and fine-search optimization via `std.Thread.Pool` (#46)
- SIMD Vec3 operations using `@Vector(4, T)` for dot/sub/add (#46)
- Schraudolph's fast exp approximation for contact scoring (#46)
- Structure-of-Arrays (SoA) layout for cache-efficient scoring (#46)
- Fixed-point integer arithmetic for mmCIF float formatting (#46)

#### Batch Processing
- `zreduce batch` subcommand for directory-level parallel processing (#47)
- File-level parallelism with atomic work-stealing and per-thread arena allocators (#47)
- Optional `--jsonl` aggregated JSONL log with mutex-protected streaming writer (#47)
- CCD dictionary loaded once and shared read-only across worker threads (#47)
- Configurable thread count (`-j N`) with auto-detect default (#47)

### Changed
- CLI refactored to subcommand structure: `zreduce run` (single file) and `zreduce batch` (directory) (#47)
- Processing pipeline extracted to `run.zig` as reusable `processFile()` (#47)
- `OptConfig` gains `n_threads` field (0 = auto-detect, batch sets to 1) (#47)
- Scoring logic extracted from `optimizer.zig` to `scoring.zig` (#46)

### Fixed
- CellList `GridTooLarge` handled with pairwise fallback instead of panic (#45)
- `writeFixedFloat` guarded against NaN/Inf with `std.fmt` fallback (#46)
- `scoreMover` buffer overflow check at runtime instead of elided `debug.assert` (#46)
- Thread pool initialization guard covers both singleton and fine-search phases (#46)

### Performance

| Structure | Before | After | Speedup |
|-----------|--------|-------|---------|
| AF-P0A9J6 (309 res) | 49s | 0.03s | 1600x |
| AF-P22523 (1486 res) | 5m17s | 0.39s | 810x |
| AF-P76347 (2339 res) | 3.1s | 1.0s | 3.2x |

Batch: 4370 E. coli proteome structures in 72s (10.5x CPU utilization).

## [0.2.0] - 2026-03-29

### Added

#### Placement Quality
- Duplicate hydrogen prevention with (name, altloc) duplicate key (#12)
- Parent atom metadata inheritance (altloc, occupancy, b_factor) for placed H (#12)
- Bond-graph-based neighbor inference replacing distance-only 1.9A cutoff (#14)
- Bond topology tables for all 20 standard amino acids (#14)
- Residue/atom-specific chemical typing (donor, acceptor, aromatic, charge) (#15)
- Terminal chemistry: C-terminal/OXT handling, N-terminal positive charge (#18)
- Conformer-aware placement with per-altloc atom lookup and blank fallback (#13)
- Insertion code (`pdbx_PDB_ins_code`) in residue boundary detection (#28)

#### Optimization Pipeline
- Mover generation from placed atoms using MoverHint metadata (#11)
- Standard residue rotator movers: OH/SH, methyl, NH3+ (#11)
- Amide flipper integration for ASN/GLN (#17)
- Histidine ring flipper integration with 6-state protonation model (#17)
- HD1 placement plan for HIS (previously only HE2) (#17)
- State-dependent donor/acceptor flags per His flip orientation (#26)
- Fine angular search after coarse optimization (#16)
- Iterative greedy for large cliques (replaces simple independent greedy) (#16)
- CellList spatial index for O(N) scoring, ~13x speedup on large structures (#35)
- CCD-based rotator generation for non-standard residues (nucleotides, glycans, ligands) (#36)
- `--no-opt` and `--no-flip` flags now functional

#### mmCIF Parsing
- Chain-break detection from `_pdbx_poly_seq_scheme` with gap-aware N/C-terminal chemistry
- `_pdbx_unobs_or_zero_occ_atoms` count for diagnostics
- Correct altloc output for added hydrogen atoms in mmCIF writer
- CIF value quoting for `[`, `]`, `{`, `}`, bare `?` and `.`

#### Quality & Diagnostics
- Model validation: sentinel position, NaN, Inf coordinate detection
- `--validate` CLI flag for detailed diagnostics
- Mover generation logging with skip count and per-skip warnings
- `isAbsentH` unified sentinel detection across codebase

### Fixed
- Flipper absent H atoms (sentinel position) excluded from output
- CIF quoting for SMILES containing `[` (e.g., `[Na+]`)
- Terminal annotation flag merge (OR) instead of type replacement
- Conformer iteration override of meta.altloc for shared backbone atoms

### Changed
- `applyChemistry` now detects chain-break terminals in addition to chain boundaries
- `executePlan` uses bond topology with distance-based fallback
- `generateMovers` accepts CCD dictionary for non-standard residue mover creation
- `scoreMover` uses CellList neighbor queries instead of all-atom iteration

## [0.1.0] - 2026-03-26

### Added

#### Phase 1: Foundation
- CIF tokenizer and parser with zero-copy design (Document/Block/Loop/Pair)
- CIF value helpers (null detection, float/int parsing with uncertainty suffix)
- Vec3 math library with rotation, dihedral angle, and angle computation
- Element table with comptime VDW radii and atom flags (45 atom types)
- Molecular model structs (Atom, Residue, Chain, Bond, Model)
- Spatial neighbor list using counting-sort cell grid (CellList)
- mmCIF atom_site extraction with chain/residue boundary tracking

#### Phase 2: Hydrogen Placement
- Type 1-6 hydrogen placement geometry functions
- Hardcoded placement plans for all 20 standard amino acids
- Streaming CCD component dictionary parser (~40K blocks without full Document tree)
- CCD-derived placement via bond topology hybridization analysis
- Unified placer: standard plans for 20 AA, CCD fallback for non-standard residues

#### Phase 3: Optimization Engine
- Dot sphere generation with concentric ring algorithm (matching original reduce)
- Dot-sphere bump/H-bond scorer with contact/clash/H-bond classification
- OH/SH rotation movers (12 orientations at 30 degree intervals)
- NH3+ rotation movers (3 orientations)
- Methyl rotation movers (3 orientations)
- Asn/Gln amide flip movers (2 orientations, penalty 0.0/0.5)
- His ring flip movers (6 orientations: 2 flip states x 3 protonation states)
- Interaction graph construction and connected component detection
- Clique optimizer: singleton, brute-force (up to 100K states), greedy fallback

#### Phase 4: Output + CLI
- mmCIF writer with atom_site loop and custom _zreduce_log category
- JSON log writer for optimization results
- Full CLI pipeline: parse args, read mmCIF, load CCD, place H, write output
- CLI flags: -o (output), -d (CCD dict), --json (log), --no-opt, --no-flip
- End-to-end integration tests (placement, multi-chain, bond lengths, round-trip)
