# zhyd

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

zhyd is a Zig implementation of the [reduce](https://github.com/rlabduke/reduce)
hydrogen placement tool for macromolecular structures. It reads mmCIF/PDB files,
adds hydrogens, optimizes supported rotatable/flippable groups, and writes mmCIF
or PDB output.

Detailed design and feature documentation will live in separate docs.

## Requirements

- [Zig](https://ziglang.org/) 0.16+

## Build and test

```bash
zig build
zig build -Doptimize=ReleaseFast
zig build test --summary all
```

The optimized executable is written to `zig-out/bin/zhyd`.

### Optional PROPKA engine

`zhyd pka` can use PROPKA 3.5.1 through
[zpropka](https://github.com/N283T/zpropka), a Zig port that reproduces
PROPKA's output. zpropka is LGPL-2.1, so it is not part of the default
(MIT-only) build. Check it out next to this repository and enable it:

```bash
git clone https://github.com/N283T/zpropka ../zpropka
zig build -Dpropka=true -Doptimize=ReleaseFast
```

## Usage

```bash
# Single structure
zhyd run input.cif -o output.cif
zhyd run input.pdb -o output.pdb

# Compressed input and output
zhyd run input.cif.gz -o output.cif
zhyd batch input_dir/ -o output_dir/ --gz

# Use a CCD dictionary and optional ligand topology
zhyd run input.cif -d components.cif -o output.cif
zhyd run input.cif -d components.cif --sdf ligand.sdf -o output.cif

# Placement controls
zhyd run input.cif -o output.cif --no-opt
zhyd run input.cif -o output.cif --no-flip
zhyd run input.cif -o output.cif --water
zhyd run input.cif -o output.cif --strip-h

# Diagnostics and control files
zhyd run input.cif -o output.cif --validate
zhyd run input.cif -o output.cif --json log.json
zhyd run input.cif -o output.cif --state-json protonation-state.json
zhyd run input.cif -o output.cif --protonation protonation.txt
zhyd run input.cif -o output.cif --ph 7.0   # states from predicted pKa (experimental)
zhyd run input.cif --dump-movers movers.txt --no-opt
zhyd run input.cif -o output.cif --fix fix.txt

# Batch processing
zhyd batch input_dir/ -o output_dir/
zhyd batch input_dir/ -d components.cif --jsonl log.jsonl
zhyd batch input_dir/ -j 4
zhyd batch input_dir/ -o output_dir/ --timing   # per-stage timing to stderr

# Precompile a CCD dictionary
zhyd compile-dict components.cif -o components.zccd

# Experimental pKa report (does not change coordinates)
zhyd pka input.cif -o pka.json
zhyd pka input.cif --model-only -o model-pka-baseline.json
zhyd pka input_dir/ -o pka_dir/ -j 8 --jsonl pka.jsonl   # one <stem>.json per file
```

## Commands

```text
zhyd run [OPTIONS] <input.cif|input.pdb>
zhyd batch [OPTIONS] <input_dir>
zhyd pka [OPTIONS] <input.cif|input.pdb|input_dir>
zhyd compile-dict [OPTIONS] <input.cif>
```

Use `zhyd <command> --help` for the full option list.

Common `run` and `batch` options:

| Option | Description |
| --- | --- |
| `-d, --dict PATH` | CCD dictionary |
| `-s, --sdf PATH` | SDF/MOL ligand topology |
| `-o, --output PATH` | Output file or directory |
| `--json PATH` / `--jsonl PATH` | Single-file JSON log / batch JSONL log |
| `--state-json PATH` | Applied residue state/formal-charge sidecar (`run` only) |
| `--protonation PATH` | Residue protonation override file |
| `--ph VALUE` | Protonation states from predicted pKa at this pH (experimental; overrides win) |
| `--fix PATH` | Force mover states from a control file |
| `--no-opt` | Skip optimization |
| `--no-flip` | Disable Asn/Gln/His flips |
| `--water` | Add water hydrogens |
| `--strip-h` | Remove existing hydrogens before placement |
| `--bond-mode MODE` | `neutron` or `xray` bond lengths |
| `--isotope NAME` | `hydrogen`/`h` or `deuterium`/`d` |
| `--nterm MODE` | `auto`, `aggressive`, or `neutral` |
| `--model VALUE` | `all` or a model number |
| `--timing` | Per-stage timing breakdown and wall time on stderr (also `pka`) |

Batch-only options:

| Option | Description |
| --- | --- |
| `-j, --threads N` | Thread count |
| `--quiet` | Suppress progress output |
| `--gz` | Write gzip-compressed output (`.cif.gz`) |

`--timing` in batch mode sums each stage over the files (so stage sums can
exceed the wall time with several threads) and adds the file count and
throughput. Machine-readable `<CMD>_<STAGE>_TIME_MS:` lines follow the table.

Given a directory, `zhyd pka` writes `<output_dir>/<stem>.json` (default
`<input_dir>_pka/`) for every `.cif`, `.pdb` or `.ent` file, optionally
gzipped. Each report is byte-identical to `zhyd pka <input_dir>/<file> -o
<output_dir>/<stem>.json` with the same options. Parameter files are loaded
once and files run in parallel (`-j N`; `--engine propka` always runs on one
thread). `--jsonl PATH` logs one line per file (`input`, `output`, `status`,
`n_sites`, `seconds`, `error`). A failing file does not stop the batch, but the
exit status is non-zero; when two files share a stem (`x.cif`, `x.cif.gz`) the
first in name order is processed and the other is reported as failed.

## Control files

### Protonation overrides

One override per line:

```text
# chain:auth_seq[:ins_code] comp_id state
A:57 HIS HIE
A:102 ASP OD2
B:14 GLU DEPROTONATED
C:88 LYS NEUTRAL
D:5 CYS THIOLATE
```

Supported states:

- `HIS`: `AUTO`, `HID`, `HIE`, `HIP`
- `ASP`: `DEPROTONATED`, `OD1`, `OD2`, `PROTONATED` (optimizer picks the oxygen)
- `GLU`: `DEPROTONATED`, `OE1`, `OE2`, `PROTONATED` (optimizer picks the oxygen)
- `LYS`: `CHARGED`, `NEUTRAL`
- `CYS`: `THIOL`, `THIOLATE`

Applied state reports and experimental pKa reports use the versioned
`zhyd.protonation.v1` schema in `docs/schema/`.

`run --ph VALUE` applies the macrostates implied by the empirical engine's
pKa at that pH: HIS → HIP, ASP/GLU → protonated, LYS → neutral and CYS →
thiolate when the predicted pKa crosses the pH (TYR, ARG and termini keep
their defaults). The optimizer still picks tautomers, ring flips and the
protonated carboxyl oxygen; `--protonation` entries win. It is experimental:
see [Phase 4](docs/specs/zhyd-protonation-pka.md) for the rule and its
evaluation.

With `-Dpropka=true`, `pka` uses PROPKA 3.5.1 by default (`--engine propka`).
Its values match `propka3` run on the same structure to the two decimals
PROPKA prints. Crystal waters are dropped before PROPKA, as PDB2PQR does;
`--propka-keep-waters` keeps them for an exact reproduction of PROPKA run on
the deposited file. Without PROPKA, `pka` uses `--engine empirical`, an
experimental MIT reimplementation of the PROPKA 3 protein model, written from
[a specification](docs/specs/zhyd-pka-empirical-engine.md) and available in
every build; on protein-only PKAD-R structures 99.84 % of its sites are
within 0.02 of PROPKA (ligand and nucleic-acid groups are not modelled yet).
By default it runs in zhyd mode: corrections of PROPKA's model quirks
and a CYS/HIS recalibration (`zhyd-empirical-v1.3`, experimental and
unvalidated; on a PKAD-R-disjoint external set v1.1 had MAE 1.10 against
1.32 for PROPKA and 1.21 for the null model);
`--empirical-mode parity` reproduces the PROPKA 3 model instead;
`--empirical-fixes` selects corrections one by one and `--parameters`
overrides per-site model pKa and desolvation scale in zhyd mode.
`--engine builtin` selects the earlier experimental built-in model below
(also used by `--model-only`).

The built-in parameter set, `zhyd-pka-protein-v0.2`, was fitted on a
cluster-split PKAD-R training set. On held-out proteins it is slightly better
than intrinsic model pKa values alone (MAE 0.80 vs 0.85), but it does poorly on
deeply buried groups such as the SNase stress set (MAE 1.95, against 0.93 for
PROPKA). With either engine, disulfide, covalently modified and
metal-coordinated residues are listed as `excluded_sites` instead of being
predicted. See `docs/specs/zhyd-protonation-pka.md` and `scripts/pka/`.

### Mover fixes

One override per line:

```text
# chain:auth_seq[:ins_code] comp_id target value
A:57 ASN amide FLIP
A:88 HIS his HID_FLIP
B:14 SER OG 6
```

Targets and values:

- `amide`: `ORIGINAL`, `FLIP`
- `his`: `HIE`, `HID`, `HIE_FLIP`, `HID_FLIP`; for a pKa-derived HIP, `0`
  (original ring) or `1` (flipped)
- `carboxyl` (protonated Asp/Glu with `PROTONATED`): `0` (H on OD1/OE1) or
  `1` (H on OD2/OE2)
- rotators: center atom name plus coarse orientation index, for example `SER OG 6`

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for release history.

## Contributing

Contributions are welcome. Please keep changes focused, follow the existing code
style, and run `zig fmt .` plus `zig build test` before submitting a PR.

## License

zhyd is released under the MIT License. See [LICENSE](LICENSE) for details.

### Attribution

zhyd is an independent from-scratch Zig reimplementation. No source code from
the original [reduce](https://github.com/rlabduke/reduce) tool is used, but the
hydrogen placement algorithm, scoring heuristics, and bond topology tables it
pioneered are owed to J. Michael Word, Duke University / UCSF, and contributors.

Please cite the original work if you use zhyd in research:

> Word, et al. (1999) "Asparagine and glutamine: using hydrogen atom contacts in
> the choice of side-chain amide orientation." J. Mol. Biol. 285, 1735-1747.
> [doi:10.1006/jmbi.1998.2401](https://doi.org/10.1006/jmbi.1998.2401)

## References

- [Original reduce (C++)](https://github.com/rlabduke/reduce)
- [Word et al. (1999) J Mol Biol 285:1735-1747](https://doi.org/10.1006/jmbi.1998.2401)
