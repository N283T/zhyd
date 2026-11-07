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

## Engine parity (zhyd reports)

Compare two `zreduce pka` reports on the same input, e.g. the zpropka oracle
against the empirical engine
([spec](../../docs/specs/zhyd-pka-empirical-engine.md)):

```sh
zreduce pka 1abc.cif --engine propka -o oracle/1abc.json
zreduce pka 1abc.cif --engine empirical -o candidate/1abc.json
python3 scripts/pka/compare_engines.py --oracle oracle/1abc.json --candidate candidate/1abc.json -o parity/1abc.csv
```

It prints the share of sites within `--tolerance` (default 0.02), the largest
difference in the predicted value and in each term, and sites present in only
one report; the exit status is 1 when the parity gate fails.

## Experimental evaluation (PKAD-R)

Do not commit redistributed PKAD-R data. Download the standalone JSON file
from the official PKAD-R site (no license text is published; the authors ask
that the paper be cited), then:

```sh
python3 scripts/pka/pkadr.py normalize PKAD-R-250211.json -o canonical.csv
python3 scripts/pka/pkadr.py fetch canonical.csv --structures structures/
python3 scripts/pka/pkadr.py predict canonical.csv --structures structures/ \
    --reports reports/ -o predictions.csv [--engine builtin --parameters params.json]
python3 scripts/pka/evaluate.py predictions.csv --output metrics.json
```

`predict --engine NAME` selects the zreduce pKa engine, and `--extra-arg ARG`
(repeatable) appends options to the `zreduce pka` command, e.g.
`--engine empirical --extra-arg=--empirical-fixes --extra-arg=Q2,Q7` for one
ablation run of the empirical engine.

`calibrate_empirical.py` fits one site kind of the empirical engine (zhyd
mode) on the train split: an exhaustive grid over the kind's model pKa and
desolvation scale, with leave-one-cluster-out CV deciding which parameters
are kept. `--cache` stores the per-row predictions of every grid point, so
re-running the selection costs nothing; `--fixed` holds kinds already
fitted. The record of the v1.1 fit is `params/zhyd-empirical-v1.1.json`.

`normalize` keeps numeric, PDB-mapped entries for supported sites and records
what it excluded. `cluster_id` joins PDB ids linked by PKAD-R's
`Sequence_Identity_>_30%` and `Alternative_PDBs` columns or by protein name, so
related proteins never cross between train and test. SNase rows get
`holdout_family=snase` and are reported only as `stress:snase`.

`evaluate.py` reports MAE, RMSE and pKa-shift direction per split, residue and
burial stratum (per split, never mixed), plus a null model that predicts the
intrinsic model pKa.

## External test set (disjoint from PKAD-R)

PKAD-R train and test have both been looked at during development, so a
separate set is kept for final checks. `external_set.py` builds it from
PKAD-3_V2 and the time-stamped blind set of
[KaML-ESM](https://github.com/JanaShenLab/KaML-ESM) (CC BY-NC 4.0; download
them yourself, do not commit them), dropping every PKAD-R structure and every
structure that shares an RCSB 30 % sequence cluster
([clusters-by-entity-30.txt](https://cdn.rcsb.org/resources/sequence/clusters/clusters-by-entity-30.txt))
with one:

```sh
python3 scripts/pka/external_set.py --pkad3 PKAD-3_V2.csv \
    --blind time_stamped/blind_test.csv --pkadr PKAD-R-250211.json \
    --clusters clusters-by-entity-30.txt -o external.csv
python3 scripts/pka/pkadr.py fetch external.csv --structures ext_structures/
python3 scripts/pka/pkadr.py predict external.csv --structures ext_structures/ \
    --reports ext_reports/ -o ext_predictions.csv
python3 scripts/pka/evaluate.py ext_predictions.csv --output ext_metrics.json
```

Rows carry `holdout_family=external`, so `evaluate.py` reports them only as
`stress:external`. Use it sparingly: every look makes it less independent.

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
