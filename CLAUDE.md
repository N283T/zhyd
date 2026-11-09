# CLAUDE.md — zhyd project instructions

## What is zhyd?

A Zig implementation of the [reduce](https://github.com/rlabduke/reduce) hydrogen placement tool for macromolecular structures (mmCIF format). Reads mmCIF, adds hydrogen atoms, optimizes orientations via clique-based scoring.

## Build & Test

```bash
zig build                           # Debug build
zig build -Doptimize=ReleaseFast    # Optimized build
zig build test --summary all        # Run all tests (~490)
```

## Run

```bash
# Single file
./zig-out/bin/zhyd run input.cif -o output.cif
./zig-out/bin/zhyd run input.cif -d components.cif -o output.cif  # with CCD
./zig-out/bin/zhyd run input.cif -d components.cif --sdf ligand.sdf -o output.cif  # with SDF
./zig-out/bin/zhyd run input.cif -o output.cif --validate         # with diagnostics

# Batch processing (parallel)
./zig-out/bin/zhyd batch input_dir/ -o output_dir/
./zig-out/bin/zhyd batch input_dir/ -d components.cif --jsonl log.jsonl
./zig-out/bin/zhyd batch input_dir/ -j 4  # limit to 4 threads
./zig-out/bin/zhyd batch input_dir/ --timing  # per-stage timing to stderr

# pKa reports (single file or directory batch)
./zig-out/bin/zhyd pka input.cif -o pka.json
./zig-out/bin/zhyd pka input_dir/ -o pka_dir/ -j 8 --jsonl pka.jsonl
```

## Project structure

```
src/
  main.zig              CLI entry point (subcommand dispatch: run/batch)
  run.zig               Single-file processing pipeline
  batch.zig             Batch processing (parallel, JSONL log)
  thread_pool.zig       Generic ThreadPool(Context, Result) (atomic chunk counter)
  timing.zig            Per-stage timing for --timing
  pka_run.zig           Per-file pKa pipeline (shared by single-file and batch)
  pka_batch.zig         `pka <dir>` batch mode (parallel, JSONL log)
  root.zig              Library re-exports
  validate.zig          Post-placement model validation
  math.zig              Vec3, rotation, dihedral
  element.zig           AtomType, VDW radii, AtomFlags, mergeFlags
  gzip.zig              Gzip I/O via Zig std.compress.flate
  integration_test.zig  End-to-end pipeline integration tests
  real_file_test.zig    End-to-end tests with real PDB structures
  cif.zig               CIF module re-exports
  cif/                  CIF parser (zero-copy tokenizer)
    char_table.zig      Character classification table
    parser.zig          CIF document parser
    tokenizer.zig       Zero-copy CIF tokenizer
    types.zig           CIF data types (Document, Block, Loop, Pair)
    value.zig           CIF value helpers (null, float, int parsing)
  mmcif.zig             _atom_site + _pdbx_poly_seq_scheme + _pdbx_unobs_or_zero_occ_atoms
  mmcif/
    conn.zig            _struct_conn and _pdbx_entity_branch_link parsing
    inline_comp.zig     Inline component dictionary and leaving atom flags
  ccd.zig               CCD component dictionary (streaming parser)
  ccd_binary.zig        Binary format for fast CCD load/save
  sdf.zig               SDF/MOL V2000 parser (non-CCD ligand topology)
  pdb.zig               PDB format parser (ATOM/HETATM records)
  model.zig             Model module re-exports
  model/                Atom, Residue, Chain, Bond, CellList
    atom.zig            Atom struct and helpers
    bond.zig            Bond struct
    chain.zig           Chain struct
    fixed_string.zig    Fixed-capacity space-padded PDB-style identifiers
    model.zig           Aggregate model container
    neighbor.zig        Spatial cell list for neighbor lookups
    residue.zig         Residue struct
  place.zig             Place module re-exports
  place/                Hydrogen placement
    placer.zig          Main placement logic (conformer-aware, altloc)
    placer_test.zig     Placement pipeline integration tests
    standard.zig        20 AA placement plans + MoverHint
    ccd_derive.zig      CCD-derived placement plan generation
    distance_derive.zig Distance-based bond inference fallback
    nucleotide.zig      DNA/RNA nucleotide placement plans
    modified.zig        Modified amino acid placement plans (MSE, SEP, etc.)
    topology.zig        Bond topology tables for 20 AAs
    chemistry.zig       Residue/atom-specific chemical annotations
    geometry.zig        Hydrogen placement geometry functions (Types 1-6)
    execute.zig         Plan execution and geometry dispatch
    lookup.zig          Atom lookup utilities for placement
    bond_policy.zig     Bond length policies (X-ray vs neutron)
    protonation.zig     Protonation state overrides
    terminal.zig        N-terminal and 3'-terminal H placement
    water.zig           Water hydrogen placement
  optimize/             Optimization engine
    optimize.zig        Optimize module re-exports
    optimizer.zig       Clique-based search (tabulated brute force / elimination) + fine angular search
    elimination.zig     Exact max-sum variable elimination for cliques beyond brute force
    scoring.zig         Mover scoring: CellList + SoA layout + mover neighbor lists
    scorer.zig          reduce-style dot scoring; pairDotScore integrates it per atom pair
    exclusion.zig       Atoms within 3 bonds of a mover atom (distance-inferred), not scored
    mover.zig           Mover struct, Orientation, isAbsentH
    rotator.zig         OH/SH, NH3+, NH2 rotators + fine orientations (methyls stay staggered)
    flipper.zig         Asn/Gln amide flip, His ring flip
    mover_gen.zig       Mover generation from placed atoms (standard + CCD)
    fix.zig             Mover state override from fix file
    clique.zig          Exact interaction graph (any orientation within scoring range) + components
    dot_sphere.zig      Concentric ring dot generation
  writer.zig            Writer module re-exports
  writer/               mmCIF + PDB + JSON output
    mmcif_writer.zig    mmCIF format output
    pdb_writer.zig      PDB format output
    json_writer.zig     JSON optimization log output
    format.zig          Value formatting and fixed-point float helpers
  test_data/            Test CIF/PDB fixtures
examples/
  data/                 Input structures (AF models, fold_test2)
  result/               zhyd output
```

## Key conventions

- **Atom names**: PDB 4-char space-padded (`[4]u8`). Use `nameSlice()` for trimmed view.
- **Altloc**: `' '` = blank/shared. Conformer-aware: `findAtomPos(mdl, res, name, target_altloc)`.
- **MoverHint**: Stored on Atom from PlacementPlan. Drives mover generation.
- **ParentMeta**: Inherits altloc/occupancy/b_factor from parent heavy atom to placed H.
- **Bond topology**: `topology.zig` for standard AAs, CCD for non-standard, SDF or distance-based inference as fallback.
- **Chemistry**: `chemistry.zig` for donor/acceptor/aromatic/charge annotations.
- **Terminal detection**: chain boundary + `is_chain_break_before` from `_pdbx_poly_seq_scheme`.
- **Sentinel**: `ABSENT_H_POS = (1000,1000,1000)` for flipper absent H. Use `isAbsentH()`.

## Testing

- `zig build test` runs all tests (~490)
- Test fixtures in `src/test_data/`: tiny.cif, tiny.pdb, multi_chain.cif, multi_chain.pdb, multi_model.cif, multi_model_null.cif, multi_model_with_h.cif, ala_with_h.cif, ala_altloc.cif, ala_stretched.cif, ala_cterm.cif, gap_chain.cif, ins_code.cif, asn.cif, his.cif, hetatm.pdb, disulfide.cif, disulfide_with_hbond.cif, branch_link.cif, entity_type.cif, inline_comp.cif, leaving_atom.cif, nterm_disorder.cif, unknown_ligand.cif
- Real-world test data in `examples/data/` (AF models, fold_test2 with DNA/RNA/glycan)

## Performance

- ReleaseFast build for benchmarks
- CellList spatial index + precomputed mover neighbor lists + SoA layout for scoring
- SIMD Vec3 + fast exp approximation
- Singleton and fine-search optimization run sequentially (optimizer-internal parallelism is a TODO since the Zig 0.16 migration)
- ~0.03s for 309-residue, ~1.0s for 2339-residue (Apple Silicon)
- Batch: 4370 E. coli structures in 72s (file-level parallelism)

## Open areas

- CCD dihedral estimation uses fixed heuristics (not computed from ideal coords)
- Distance-based bond inference cannot detect aromatic rings; bond order promotion uses valence heuristics only
- State-dependent chemistry updates for flip movers are per-orientation flags only
