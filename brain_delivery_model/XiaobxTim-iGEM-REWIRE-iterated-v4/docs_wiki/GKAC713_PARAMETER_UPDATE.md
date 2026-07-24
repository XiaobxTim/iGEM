# gkac713-informed v4 parameter update

This document records how the v4 model uses Han et al. 2022, Nucleic Acids
Research, DOI: `10.1093/nar/gkac713`, as a parameter-design reference.

The changes are implemented as literature-informed priors, not as final
calibrated constants. They should be updated when our own APOE editing,
expression and off-target experiments become available.

## Why this paper is relevant

The paper describes REWIRE, a gRNA-free programmable RNA base-editing platform.
For our project, the most relevant construct class is CU-REWIRE:

- RNA targeting is provided by engineered PUF domains.
- C-to-U editing is performed by APOBEC3A/A3A-family deaminase activity.
- HEK293T transfection assays commonly use a 48 h readout.
- Expanded 10-repeat PUF recognition improves RNA targeting specificity.
- C-to-U editing is strongest near a narrow editing window and is strongly
  sequence-context dependent.

This maps directly to our v4 model:

PUF design -> sequence/accessibility score -> effective `k_on` and `k_cat` ->
APOE112/APOE158 editing -> separate off-target burden classes.

## Implemented parameter/rule changes

### 1. 48 h as the HEK293T reference readout

The v4 config now uses:

- `simulation.default_t_end: 48.0`

The v4 design screen now compares:

- 24 h
- 48 h
- 72 h

The utility function includes a small penalty for choosing later-than-48 h
readouts. This prevents the model from recommending a late time point simply
because the cumulative ODE state had more time to increase.

Important interpretation note: the current v4 simulator still includes the
full delivery-expression chain, so its best design-screen time can still shift
toward 72 h when expression is delayed. That should be reported as a full-chain
prediction, not as a contradiction of the gkac713 HEK293T 48 h assay design.
For a pure HEK293T transfection-only digital twin, the upstream AAV/BBB modules
should be bypassed or replaced by a direct transfection input.

### 2. A3A is the default C-to-U editor

The default remains:

- `editing.editor_type: A3A`

This keeps the model aligned with CU-REWIRE instead of mixing APOBEC3A and
APOBEC1 assumptions.

### 3. Stronger APOE112 prior, lower APOE158 prior

The therapeutic objective in this model is APOE4-to-APOE3-like editing. That
means codon 112 is the preferred productive edit, while codon 158 is tracked as
a double-edit/APOE2-like risk channel unless the wet-lab strategy explicitly
chooses dual-site editing.

Updated defaults:

- `editing.k_cat_112: 1.0`
- `editing.k_cat_158: 0.15`

These are normalized effective turnover priors. They should be fitted with
APOE112/APOE158 amplicon time-course data.

### 4. Narrow editing-window rule

In `models/editing/sequence_to_kinetics.py`, context and distance now map to
`k_cat` as:

```text
kcat_scale = (0.04 + 0.96 * UC_context_score) * exp(-0.9 * |distance - 2|)
```

Interpretation:

- editing is strongest when the edited C is close to the preferred window
  around position +2 after the PUF binding site;
- UC-like context receives high catalytic weight;
- poor context or wrong spacing sharply lowers effective catalytic rate.

### 5. 10R PUF specificity rule

The v4 sequence-to-kinetics layer now calculates:

```text
puf_offtarget_scale = exp(-1.4 * max(PUF_repeats - 8, 0)) * exp(0.45 * mismatch_count)
```

Interpretation:

- 10R designs receive a much lower PUF-mediated off-target prior than 8R
  designs;
- mismatches in the intended design reduce confidence and increase aggregate
  off-target risk.

This modifier rescales:

- `editing.k_on_puf_off`
- `editing.k_cat_puf_off`

### 6. Separate deaminase background from PUF mismatch

Updated defaults:

- `editing.k_on_puf_off: 0.015`
- `editing.k_cat_puf_off: 0.006`
- `editing.k_deaminase_bg: 0.0003`

The design table controls active versus inactive editors using:

- `editor_activity_scale`

For inactive deaminase and PUF-only controls, this value is `0.0`, which sets
active editing and deaminase background to zero.

## Files changed

- `config/base_config.yaml`
- `config/parameter_provenance.yaml`
- `models/editing/sequence_to_kinetics.py`
- `scripts/run_iterated_analysis.py`
- `wetlab/templates/puf_design_candidates_template.csv`
- `tests/test_sequence_to_kinetics.py`

## What wet-lab data should update these priors

The most important update data are:

- APOE112 editing percentage at 24/48/72 h.
- APOE158 editing percentage at 24/48/72 h.
- Linked APOE112/APOE158 amplicon or long-read data to distinguish APOE3-like
  from APOE2-like products.
- Local bystander editing around the APOE target window.
- Top predicted PUF-mismatch off-target amplicon panel.
- RNA-seq for transcriptome-wide off-target validation.
- Mock, PUF-only and inactive-A3A controls.
- Editor mRNA qPCR and protein Western/fluorescence time courses.

## What remains literature-prior only

Until more data are available, the following should not be treated as measured
constants:

- exact `k_on`, `k_off` and `k_cat` separation;
- the numerical slope linking UC context to `k_cat`;
- the numerical slope linking 10R PUF to off-target reduction;
- deaminase-only background editing rate;
- RNA accessibility score calibration.

These are useful priors for design ranking, but the wiki should present them as
model assumptions that are updated through DBTL cycles.
