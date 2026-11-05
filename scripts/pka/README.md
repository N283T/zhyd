# pKa validation harness

These scripts keep correctness oracles and experimental data outside the MIT
implementation. They use only the Python standard library.

## PROPKA oracle

Install PROPKA separately under its own LGPL-2.1 terms, then run:

```sh
python3 scripts/pka/run_propka_oracle.py structures/*.pdb --output-dir oracle/
zreduce pka structures/1abc.pdb -o predictions/1abc.json
python3 scripts/pka/compare_propka.py --zhyd predictions/1abc.json --propka oracle/1abc.pka -o comparisons/1abc.csv
```

The harness invokes the executable as a subprocess and parses only its public
output. No PROPKA source or parameters are imported.

## Experimental evaluation (PKAD-R)

Do not commit redistributed PKAD-R data. Download the standalone JSON file
from the official PKAD-R site (no license text is published; the authors ask
that the paper be cited), then:

```sh
python3 scripts/pka/pkadr.py normalize PKAD-R-250211.json -o canonical.csv
python3 scripts/pka/pkadr.py fetch canonical.csv --structures structures/
python3 scripts/pka/pkadr.py predict canonical.csv --structures structures/ \
    --reports reports/ -o predictions.csv [--parameters params.json]
python3 scripts/pka/evaluate.py predictions.csv --output metrics.json
```

`normalize` keeps numeric, PDB-mapped entries for supported sites and records
what it excluded. `cluster_id` joins PDB ids linked by PKAD-R's
`Sequence_Identity_>_30%` and `Alternative_PDBs` columns or by protein name, so
related proteins never cross between train and test. SNase rows get
`holdout_family=snase` and are reported only as `stress:snase`.

`evaluate.py` reports MAE, RMSE and pKa-shift direction per split, residue and
burial stratum (per split, never mixed), plus a null model that predicts the
intrinsic model pKa.

## Calibration

```sh
python3 scripts/pka/calibrate.py canonical.csv --structures structures/ \
    --version my-fit -o fit.json [--objective mae|rmse]
```

Coordinate descent over a bounded grid reads only training-split rows. The
fitted set used as the built-in default, with data hashes and held-out
results, is in `params/zhyd-pka-protein-v0.2.json`.

Sources: the [PKAD-R site](https://compbio.clemson.edu/pkad-r/), the
[PROPKA repository](https://github.com/jensengroup/propka), and the primary
papers linked from `docs/specs/zhyd-protonation-pka.md`.
