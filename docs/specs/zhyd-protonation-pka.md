# zhyd protonation and pKa architecture

Status: implementation phase 1; experimental pKa parameters; schema v1.

## Scope and constraints

The long-term component prepares biomolecular structures by making three
outputs explicit: protonation macrostate, formal charge, and hydrogen
coordinates. The current executable and package remain named `zreduce`; `run`
and `batch` behavior remains compatible.

The implementation is independent and MIT-licensed. PROPKA is LGPL-2.1 and is
used only as an executable correctness oracle. No PROPKA code, control flow, or
parameter tables were copied or translated. The method is based on scientific
concepts described by the primary papers, especially continuous burial/
desolvation, hydrogen bonds, electrostatics, and coupled groups:

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

## Protein v0.1 model

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

All parameters are centralized under `zhyd-pka-protein-v0.1`. Bounds prevent a
malformed structure producing unbounded values. Low confidence is assigned to
deep burial or strong coupling. These are scientific priors, not fitted or
validated parameters; numerical pKa values must be treated as experimental.

### Assumptions and fallbacks

- The first blank or `A` conformer is used; no ensemble averaging yet.
- Missing required site atoms omit the site and increment `omitted_sites`.
- Missing entity metadata falls back only for the 20 standard amino-acid names.
- C termini use OXT when present, otherwise backbone O.
- Neutral HIS prediction does not guess HID versus HIE; geometry selection is
  required and a warning is emitted.
- Existing manual overrides remain authoritative for coordinate placement and
  are recorded as `user_override` by `--state-json`.

### Non-goals for v0.1

Ligand/nucleic-acid pKa, metal coordination, covalent ligand induction, constant-
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

### Phase 3 — reproducible validation (infrastructure complete; data run pending)

- External PROPKA runner/parser and canonical evaluator are committed.
- Leakage-resistant cluster split and required metrics are implemented.
- Pending: pin a permitted PKAD-R snapshot, define clusters, run baselines, and
  publish results; run SNase separately.

### Phase 4 — pH-driven state application (not started)

- Requires accepted Phase 2/3 thresholds, explicit user opt-in, state-to-
  placement integration tests for every supported kind, and charge/H-coordinate
  consistency validation. Only then may a `--ph` assignment interface be added.

### Phase 5 — microstates and new chemistries (not started)

- Exact/sampled coupled microstates with normalized probabilities and ensemble
  provenance; then ligand, nucleic-acid, and metal-aware providers under separate
  parameter versions and validation sets.
