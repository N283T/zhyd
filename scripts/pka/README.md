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

## Experimental evaluation

Do not commit redistributed PKAD-R data without first checking its data-use
terms. Download the standalone file from the official PKAD-R site and normalize
it to CSV with these required columns:

`comp_id,experimental_pka,predicted_pka,model_pka,burial_fraction,cluster_id`

Then run:

```sh
python3 scripts/pka/evaluate.py canonical.csv --output metrics.json
```

`cluster_id` must represent a protein/sequence cluster, not an individual row.
The deterministic hash split prevents members of one cluster crossing between
train and test. The report includes MAE, RMSE, residue type, surface/buried
strata, and pKa-shift direction accuracy. Keep SNase mutants in a named,
holdout-only cluster family and report them separately as a stress test.

Sources: the [PKAD-R site](https://compbio.clemson.edu/pkad-r/), the
[PROPKA repository](https://github.com/jensengroup/propka), and the primary
papers linked from `docs/specs/zhyd-protonation-pka.md`.
