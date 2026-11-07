# zreduce

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

zreduce is a Zig implementation of the [reduce](https://github.com/rlabduke/reduce)
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

The optimized executable is written to `zig-out/bin/zreduce`.

### Optional PROPKA engine

`zreduce pka` can use PROPKA 3.5.1 through
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
zreduce run input.cif -o output.cif
zreduce run input.pdb -o output.pdb

# Compressed input and output
zreduce run input.cif.gz -o output.cif
zreduce batch input_dir/ -o output_dir/ --gz

# Use a CCD dictionary and optional ligand topology
zreduce run input.cif -d components.cif -o output.cif
zreduce run input.cif -d components.cif --sdf ligand.sdf -o output.cif

# Placement controls
zreduce run input.cif -o output.cif --no-opt
zreduce run input.cif -o output.cif --no-flip
zreduce run input.cif -o output.cif --water
zreduce run input.cif -o output.cif --strip-h

# Diagnostics and control files
zreduce run input.cif -o output.cif --validate
zreduce run input.cif -o output.cif --json log.json
zreduce run input.cif -o output.cif --state-json protonation-state.json
zreduce run input.cif -o output.cif --protonation protonation.txt
zreduce run input.cif --dump-movers movers.txt --no-opt
zreduce run input.cif -o output.cif --fix fix.txt

# Batch processing
zreduce batch input_dir/ -o output_dir/
zreduce batch input_dir/ -d components.cif --jsonl log.jsonl
zreduce batch input_dir/ -j 4

# Precompile a CCD dictionary
zreduce compile-dict components.cif -o components.zccd

# Experimental pKa report (does not change coordinates)
zreduce pka input.cif -o pka.json
zreduce pka input.cif --model-only -o model-pka-baseline.json
```

## Commands

```text
zreduce run [OPTIONS] <input.cif|input.pdb>
zreduce batch [OPTIONS] <input_dir>
zreduce pka [OPTIONS] <input.cif|input.pdb>
zreduce compile-dict [OPTIONS] <input.cif>
```

Use `zreduce <command> --help` for the full option list.

Common `run` and `batch` options:

| Option | Description |
| --- | --- |
| `-d, --dict PATH` | CCD dictionary |
| `-s, --sdf PATH` | SDF/MOL ligand topology |
| `-o, --output PATH` | Output file or directory |
| `--json PATH` / `--jsonl PATH` | Single-file JSON log / batch JSONL log |
| `--state-json PATH` | Applied residue state/formal-charge sidecar (`run` only) |
| `--protonation PATH` | Residue protonation override file |
| `--fix PATH` | Force mover states from a control file |
| `--no-opt` | Skip optimization |
| `--no-flip` | Disable Asn/Gln/His flips |
| `--water` | Add water hydrogens |
| `--strip-h` | Remove existing hydrogens before placement |
| `--bond-mode MODE` | `neutron` or `xray` bond lengths |
| `--isotope NAME` | `hydrogen`/`h` or `deuterium`/`d` |
| `--nterm MODE` | `auto`, `aggressive`, or `neutral` |
| `--model VALUE` | `all` or a model number |

Batch-only options:

| Option | Description |
| --- | --- |
| `-j, --threads N` | Thread count |
| `--quiet` | Suppress progress output |
| `--gz` | Write gzip-compressed output (`.cif.gz`) |

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
- `ASP`: `DEPROTONATED`, `OD1`, `OD2`
- `GLU`: `DEPROTONATED`, `OE1`, `OE2`
- `LYS`: `CHARGED`, `NEUTRAL`
- `CYS`: `THIOL`, `THIOLATE`

Applied state reports and experimental pKa reports use the versioned
`zhyd.protonation.v1` schema in `docs/schema/`. The pKa model is not yet
validated for production charge assignment and intentionally does not expose a
`--ph` coordinate-assignment option.

With `-Dpropka=true`, `pka` uses PROPKA 3.5.1 by default (`--engine propka`).
Its values match `propka3` run on the same structure to the two decimals
PROPKA prints. Crystal waters are dropped before PROPKA, as PDB2PQR does;
`--propka-keep-waters` keeps them for an exact reproduction of PROPKA run on
the deposited file. `--engine empirical` selects an experimental MIT
reimplementation of the PROPKA 3 protein model, written from
[a specification](docs/specs/zhyd-pka-empirical-engine.md) and available in
every build; on protein-only PKAD-R structures 99.84 % of its sites are
within 0.02 of PROPKA (ligand and nucleic-acid groups are not modelled yet).
`--empirical-mode zhyd` turns on the corrections of PROPKA's model quirks
adopted so far (`zhyd-empirical-v1`, experimental and unvalidated);
`--empirical-fixes` selects corrections one by one.
`--engine builtin` selects the experimental built-in model below.

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
- `his`: `HIE`, `HID`, `HIE_FLIP`, `HID_FLIP`
- rotators: center atom name plus coarse orientation index, for example `SER OG 6`

## Changelog

See [CHANGELOG.md](CHANGELOG.md) for release history.

## Contributing

Contributions are welcome. Please keep changes focused, follow the existing code
style, and run `zig fmt .` plus `zig build test` before submitting a PR.

## License

zreduce is released under the MIT License. See [LICENSE](LICENSE) for details.

### Attribution

zreduce is an independent from-scratch Zig reimplementation. No source code from
the original [reduce](https://github.com/rlabduke/reduce) tool is used, but the
hydrogen placement algorithm, scoring heuristics, and bond topology tables it
pioneered are owed to J. Michael Word, Duke University / UCSF, and contributors.

Please cite the original work if you use zreduce in research:

> Word, et al. (1999) "Asparagine and glutamine: using hydrogen atom contacts in
> the choice of side-chain amide orientation." J. Mol. Biol. 285, 1735-1747.
> [doi:10.1006/jmbi.1998.2401](https://doi.org/10.1006/jmbi.1998.2401)

## References

- [Original reduce (C++)](https://github.com/rlabduke/reduce)
- [Word et al. (1999) J Mol Biol 285:1735-1747](https://doi.org/10.1006/jmbi.1998.2401)
