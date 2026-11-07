# zhyd protonation and pKa architecture

Status: three pKa engines. PROPKA 3.5.1 via the optional zpropka port
(reference method, reproduces `propka3`), the experimental MIT empirical
engine (zhyd mode `zhyd-empirical-v1.2`, the default engine without
PROPKA; parity mode `zhyd-empirical-v1-parity`), and the earlier experimental
built-in model (`zhyd-pka-protein-v0.2`, fitted, not accepted for
production). Schema v1.

## Scope and constraints

The long-term component prepares biomolecular structures by making three
outputs explicit: protonation macrostate, formal charge, and hydrogen
coordinates. The current executable and package remain named `zreduce`; `run`
and `batch` behavior remains compatible.

## pKa engines and licensing

| Engine | Where | License | Status |
|---|---|---|---|
| `propka` | [zpropka](https://github.com/N283T/zpropka), optional dependency (`zig build -Dpropka=true`) | LGPL-2.1 (derivative of PROPKA 3.5.1) | Reproduces `propka3` 3.5.1 output; default when built in |
| `empirical` | `src/pka_empirical.zig` | MIT | Experimental; default without PROPKA. Clean-room reimplementation of the PROPKA 3 protein model from a written spec. Parity mode: 99.84 % of protein sites within 0.02 of `propka`; zhyd mode (default): PKAD-R test MAE 0.71 |
| `builtin` | `src/pka.zig` | MIT | Experimental, independently written; superseded as default, used for `--model-only` |

zpropka is a faithful Zig port of PROPKA 3.5.1, kept in its own repository so
the default zreduce build stays MIT-only. zreduce talks to it only through
`zpropka.api.predictPdb`, from `src/pka_propka.zig`, which:

- writes heavy atoms as PDB text: single-character chain aliases for long
  auth IDs, standard residues as ATOM and modified residues as HETATM (PDB
  convention, which PROPKA's ligand handling expects), and `TER   ` after
  every polymer chain (PROPKA resets N-terminus detection there);
- drops crystal waters by default, as PDB2PQR does before PROPKA (PROPKA
  itself treats HOH as a ligand). On PKAD-R this changed nothing measurable
  (test MAE 0.677 vs 0.678); `--propka-keep-waters` restores PROPKA-on-file
  behavior;
- maps PROPKA groups back to residues by (chain alias, number, insertion
  code), moves disulfide / covalently bonded / metal-coordinated residues to
  `excluded_sites`, and maps determinants onto the schema terms:
  `burial` = desolvation, `hydrogen_bond` = side-chain + backbone H-bonds,
  `electrostatic` = Coulomb, `coupling` = 0 (PROPKA's terms are additive to
  the reported pKa);
- keeps PROPKA values outside -5..20 (PROPKA occasionally produces them for
  coupled cysteines), with low confidence and `pka_outside_physical_range`.

Verified against `propka3` 3.5.1 on 185 PKAD-R structures with waters kept:
14,291 sites, max |Δ| 0.005 (PROPKA prints two decimals).

Plan for the built-in engine: describe PROPKA's method as a specification
(equations, conditions, parameters cited to the papers, not code), write an
MIT implementation from that specification, and use zpropka as a test oracle
only. Improvements (metals, covalent cysteines, altloc ensembles,
recalibration) go into the MIT engine. The specification, clean-room
protocol, parity quirks and acceptance gates are in
[zhyd-pka-empirical-engine.md](zhyd-pka-empirical-engine.md);
`scripts/pka/compare_engines.py` checks parity against the oracle. Parity
mode is implemented as `--engine empirical`; its parity numbers and the
remaining deviations are in the spec's *Parity status* section. zhyd mode
(`--empirical-mode zhyd`) applies the adopted corrections and the CYS/TYR
recalibration; their rules, the PKAD-R ablation and the fit are in the
spec's *zhyd mode* and *Recalibration* sections.

The built-in engine is based on scientific concepts described by the primary
papers, especially continuous burial/desolvation, hydrogen bonds,
electrostatics, and coupled groups:

- Olsson et al., *JCTC* 7, 525–537 (2011),
  [doi:10.1021/ct100578z](https://doi.org/10.1021/ct100578z).
- Søndergaard et al., *JCTC* 7, 2284–2295 (2011),
  [doi:10.1021/ct200133y](https://doi.org/10.1021/ct200133y).
- Chen et al., PKAD-R (2025),
  [official database](https://compbio.clemson.edu/pkad-r/) and
  [doi:10.1142/S2737416525500164](https://doi.org/10.1142/S2737416525500164).
- The SNase pKa-cooperative challenge is treated as a difficult, distribution-
  shifted stress test rather than ordinary training data; see
  [doi:10.1002/prot.23110](https://doi.org/10.1002/prot.23110).

## Current architecture findings

The pre-existing pipeline already had useful separation but did not persist a
state model:

1. `place/protonation.zig` parsed residue overrides.
2. `place/chemistry.zig` mapped those overrides to donor/acceptor and signed
   atom flags.
3. `place/placer.zig` selected placement plans and added coordinates.
4. `optimize/flipper.zig` selected neutral HIS tautomer plus ring flip;
   rotators and amide flips were separate movers.
5. Writers emitted coordinates and mover choices, but not a stable residue-
   state/formal-charge record.

This phase preserves those modules and adds a parallel pKa/state boundary
(`pka.zig` plus `writer/protonation_json.zig`). It can later be extracted as a
`zhyd` module without moving placement internals first.

## Responsibility boundaries

| Responsibility | Input | Output | Must not decide |
|---|---|---|---|
| pKa prediction | heavy-atom structure, versioned parameters | pKa and named correction terms | tautomer or H coordinates |
| Macrostate selection | pKa plus an explicit pH/policy | protonated/deprotonated state and charge | H geometry |
| Tautomer selection | macrostate plus local geometry | HID/HIE or acid oxygen | total charge |
| H placement | explicit state and topology | H coordinates | pKa |
| Flip/rotamer/water optimization | placed alternatives and scorer | orientation | intrinsic pKa |

The experimental `pka` command reports a recommended state at a documented
reference pH of 7.0 but sets `applies_to_coordinates=false`. It deliberately has
no `--ph`: pH-driven coordinate assignment will not be exposed until prediction
accuracy and state application are validated together.

## Protein pKa model (parameters v0.2)

Supported sites are ASP, GLU, HIS, LYS, CYS, TYR, and observed N/C termini. Each
prediction is an additive, inspectable decomposition:

`pKa = model + burial + hydrogen_bond + electrostatic + coupling`.

- **Model pKa:** literature-scale priors centralized in `modelPka`.
- **Burial/desolvation:** continuous heavy-atom neighbor density within 10 Å.
  Burial raises acid pKa and lowers base pKa because it disfavors charge.
- **Hydrogen bonds:** distance-weighted polar contacts stabilize the charged
  acid/base form with a bounded shift.
- **Electrostatics:** a screened distance term currently treats ARG as fixed
  +1. Supported titratable sites are not double-counted here.
- **Coupling:** a deterministic mean-field iteration evaluates the expected
  charge of nearby titratable sites at the target site's transition pH.

All parameters are centralized in `pka.Parameters`. Bounds prevent a
malformed structure producing unbounded values. Low confidence is assigned to
deep burial or strong coupling.

- `zhyd-pka-protein-v0.1` used hand-set priors (burial 2.0, H-bond 0.35 per
  contact, cap 1.5).
- `zhyd-pka-protein-v0.2` (current) fits `burial_max_shift`, `hbond_strength`
  and `hbond_max_shift` on the PKAD-R training split only. The fit sets the
  burial shift to 0: the neighbor-count descriptor did not lower training
  error. Burial fraction is still computed for confidence and strata.
  Provenance: `scripts/pka/params/zhyd-pka-protein-v0.2.json`.
- `zreduce pka --parameters PATH` loads a JSON parameter set for research or
  recalibration. Unknown fields are rejected, and a non-default set must carry
  its own `version`, which the report echoes.

### Assumptions and fallbacks

- The first blank or `A` conformer is used; no ensemble averaging yet.
- Missing required site atoms omit the site and increment `omitted_sites`.
- Residues are classified before prediction, for both engines, in this order:
  a side-chain ligand atom ≤ 2.8 Å from a metal ion → `metal_coordinated`;
  CYS SG ≤ 2.5 Å from another CYS SG → `disulfide`; CYS SG ≤ 2.1 Å from any
  other non-solvent heavy atom, or flagged covalent by `_struct_conn` with no
  partner in the model → `covalently_bonded`. Metal comes first because
  `_struct_conn` also flags metalc bonds. These residues are listed in
  `excluded_sites`. `run --state-json` reports cystines as
  `disulfide_bonded`, other covalent CYS as `covalently_bonded` (source
  `covalent_bond`), and adds `metal_coordinated_state_not_modeled` to metal
  ligands.
- Explicit solvent (water entities, HOH/WAT/DOD/H2O) is the reference state of
  the model pKa values and is ignored by the burial and hydrogen-bond terms.
- Missing entity metadata falls back only for the 20 standard amino-acid names.
- C termini use OXT when present, otherwise backbone O.
- Neutral HIS prediction does not guess HID versus HIE; geometry selection is
  required and a warning is emitted.
- Existing manual overrides remain authoritative for coordinate placement and
  are recorded as `user_override` by `--state-json`.

### Non-goals for this phase

Ligand/nucleic-acid pKa, metal-bound protonation states, covalent ligand induction, constant-
pH ensembles, alternate-conformer ensembles, and exact coupled microstate
enumeration. Their extension boundary is a new `SiteKind`/site provider and a
versioned parameter family; none are silently approximated as protein residues.

## Stable interchange schema

`docs/schema/zhyd.protonation.v1.schema.json` defines two report kinds:

- `applied_protonation_state`, produced by `run --state-json`, records the state
  that actually controls coordinates, charge, and provenance and links the
  companion H-containing structure through `coordinate_file`. Its pKa is null.
  Disulfide-bonded cysteines (flagged by `_struct_conn` or SG–SG ≤ 2.5 Å) are
  reported as `disulfide_bonded` with source `covalent_bond`; a CYS override on
  such a residue is ignored with the `override_ignored_covalent_bond` warning.
  The `pka` command does not treat them as titratable sites.
- `pka_prediction`, produced by `pka`, records each correction term, predicted
  state/charge at reference pH, confidence, and warnings. It explicitly says
  the state does not apply to coordinates.

Downstream tools must key sites by model, chain identifiers, residue identifiers,
insertion code, and site type—not by a presentation string. New optional fields
may be added within schema v1; semantic changes require v2. Parameter revisions
do not change schema and must change `parameter_version`.

## Validation design

The harness under `scripts/pka/` keeps three evidence layers distinct:

1. **Regression:** all existing placement/optimization tests must pass unchanged.
2. **Oracle comparison:** invoke installed PROPKA as an external executable and
   compare residue-by-residue output. Oracle agreement diagnoses implementation
   differences but is not experimental truth.
3. **Independent experiment:** normalize PKAD-R measurements, split by protein or
   sequence cluster, and report MAE, RMSE, residue-stratified, surface/buried,
   and pKa-shift direction metrics. Never split individual residues randomly.

SNase mutants remain a named holdout stress test. Dataset source version,
download date, normalization code revision, PDB revision, parameter version,
split seed, cluster definition, and executable versions must accompany results.

## Development results (PKAD-R, 2026-11-05)

Data: PKAD-R `PKAD-R-250211.json` (sha256 `9ee3cfb8…`), 915 numeric,
PDB-mapped entries in 78 clusters / 185 PDB entries; 868 matched to a predicted
site (the rest are mainly CYS that are disulfide/modified in the deposited
structure). Split: SHA-256 of cluster id with seed `zhyd-pka-split-v1`, 80 %
train; SNase is a separate stress family.

| Split | n | null (model pKa) MAE / RMSE | v0.1 MAE / RMSE | v0.2 MAE / RMSE |
|---|---|---|---|---|
| train | 633 | 0.76 / 1.09 | 0.87 / 1.28 | 0.69 / 1.05 |
| test | 133 | 0.85 / 1.35 | 0.91 / 1.34 | 0.80 / 1.35 |
| SNase stress | 102 | 1.96 / 2.49 | 1.33 / 1.63 | 1.95 / 2.61 |

v0.1 and v0.2 columns include the solvent-exclusion fix.

Findings:

- v0.2 beats the null model on held-out test proteins in MAE, but only by a
  small margin; RMSE is unchanged. Gains come from surface residues (test
  surface MAE 0.51 vs 0.67 for v0.1).
- Buried groups are the main failure. With the burial shift removed, SNase
  buried residues get worse (MAE 3.28 vs 1.92 for v0.1). The v0.1 burial term
  helped buried groups but hurt surface ones because its descriptor saturates
  too early. A better desolvation descriptor (e.g. SASA- or depth-based, with
  explicit compensation by H-bonds and salt bridges) is the next modeling step.
- Active-site CYS are poorly predicted (test MAE ≈ 2.3). The PKAD-R CYS set is
  dominated by low-pKa enzyme thiols that need helix-dipole and backbone-amide
  terms the model does not have.
- These numbers are development estimates, not a blind benchmark: test and
  SNase metrics of an earlier fit were inspected during development, and the
  solvent fix was found while diagnosing them. Parameters were selected on
  train rows only.

Reproduce (data are not redistributed; download PKAD-R from the official site):

```sh
python3 scripts/pka/pkadr.py normalize PKAD-R-250211.json -o canonical.csv
python3 scripts/pka/pkadr.py fetch canonical.csv --structures structures/
python3 scripts/pka/pkadr.py predict canonical.csv --structures structures/ \
    --reports reports/ -o predictions.csv
python3 scripts/pka/evaluate.py predictions.csv --output metrics.json
python3 scripts/pka/calibrate.py canonical.csv --structures structures/ \
    --version my-fit -o fit.json   # train split only
```

### PROPKA engine on PKAD-R

Same split, PROPKA 3.5.1 via zpropka (waters dropped, exclusions applied):

| Split | n | PROPKA MAE / RMSE | built-in v0.2 MAE | null MAE |
|---|---|---|---|---|
| train | 550 | 0.87 / 1.40 | 0.69 | 0.79 |
| test | 129 | 1.13 / 4.50 | 0.80 | 0.80 |
| SNase stress | 102 | 0.93 / 1.21 | 1.95 | 1.96 |

PROPKA's test error is dominated by a few CYS (test CYS MAE 15, including a
52.96 prediction PROPKA itself makes); without CYS its test MAE is about
0.68. PROPKA is clearly better on buried groups (SNase).

## Phases and acceptance criteria

### Phase 0 — behavior inventory (complete)

- Existing `run`/`batch` CLI and 487-test baseline recorded.
- Override, chemistry, placement, flipper, and writer paths mapped.

### Phase 1 — explicit state contract (complete)

- Versioned schema includes state, formal charge, provenance, optional pKa,
  correction terms, confidence, and warnings.
- `run --state-json` reports default, geometry-selected, and user-override states.
- Existing run output is unchanged when the option is absent.

### Phase 2 — staged protein predictor (implemented, not scientifically accepted)

- Model-only baseline and independently implemented structure terms are unit
  tested and deterministic.
- `pka` emits valid JSON for mmCIF/PDB and does not mutate coordinates.
- Acceptance for production remains blocked until the evaluation thresholds are
  selected before fitting and met on cluster-held-out PKAD-R data.

### Phase 3 — reproducible validation (development run complete)

- External PROPKA runner/parser and canonical evaluator are committed.
- Leakage-resistant cluster split and required metrics are implemented.
- PKAD-R snapshot pinned by hash; clusters, null baseline, v0.1/v0.2 results
  and SNase stress results recorded above.
- Pending for acceptance: improved desolvation term, a fresh blind holdout,
  and thresholds chosen before that holdout is scored. Suggested gate: held-out
  MAE at least 0.1 below the null model and buried-site MAE below 1.5.

### Phase 4 — pH-driven state application (not started)

- Requires accepted Phase 2/3 thresholds, explicit user opt-in, state-to-
  placement integration tests for every supported kind, and charge/H-coordinate
  consistency validation. Only then may a `--ph` assignment interface be added.

### Phase 5 — microstates and new chemistries (not started)

- Exact/sampled coupled microstates with normalized probabilities and ensemble
  provenance; then ligand, nucleic-acid, and metal-aware providers under separate
  parameter versions and validation sets.
