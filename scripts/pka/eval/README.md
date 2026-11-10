# Literature evaluation set (2026-11)

`literature-2026-11.csv`: 47 experimental pKa values from recent papers,
collected 2026-11-09 to check the empirical engine outside PKAD-R and the
external set. All PDB entries except HCA II (4Y0J) and ubiquitin (1UBQ) are
in no RCSB 30 % sequence cluster of the fitting data; those two are marked
`known-protein-new-site` (sites not in the fitting data). Values were checked
against the papers. It is an evaluation set only: no parameter was selected
on it.

Columns are the canonical ones of `pkadr.py` plus:

- `lead_id`: source paper (below).
- `condition`: measurement condition. L01 is listed at 0 and 300 mM NaCl;
  score one of them (`--exclude-condition 0mM_NaCl`).
- `value_qualifier`: `point`, or `upper_bound` for "pKa below the value"
  (sites that did not titrate in the measured range).

| Lead | Protein (PDB) | Sites | Reference |
|---|---|---|---|
| L01 | MaSp1 / MaSp2 NTD (7WIO, 8GS7) | ASP, GLU | Oktaviani et al. 2023, doi:10.1021/acs.biomac.2c01344 |
| L02 | *B. thermoproteolyticus* ferredoxin (7YL8) | ASP | Wada et al. 2024, doi:10.7554/eLife.102506 |
| L13 | Human carbonic anhydrase II (4Y0J) | HIS | Fisher et al. 2025, doi:10.1002/cbic.202400930 |
| L20 | LCC-ICCG PETase (6THT, 8OTA) | HIS | Charlier et al. 2024, doi:10.1016/j.bpj.2024.04.026 |
| L23 | Human ubiquitin (1UBQ) | LYS | Pawloski et al. 2022, doi:10.1021/acs.biochem.2c00085 |
| L27 | Caspase-9 CARD (4RHW, chain E) | HIS | Rawal et al. 2026, doi:10.64898/2026.06.13.732084 (preprint) |
| L32 | *B. subtilis* RNase P protein (1A6F) | HIS | Mosley et al. 2015 |
| L35 | *E. coli* GDP-mannose mannosyl hydrolase (1RYA) | HIS | Legler et al. 2002, doi:10.1021/bi020362e |

`4RHW_E` is chain E of 4RHW written as its own file (the CARD filament
subunit measured); fetch 4RHW and extract chain E before running `zhyd pka`.

```sh
python3 scripts/pka/pkadr.py fetch scripts/pka/eval/literature-2026-11.csv --structures cif/
zhyd pka cif/ -o reports/
python3 scripts/pka/score_literature.py scripts/pka/eval/literature-2026-11.csv \
    zhyd=reports/ --exclude-condition 0mM_NaCl --verbose
```
