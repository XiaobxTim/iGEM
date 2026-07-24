# Cycle 6: v4 sequence-aware APOE editing digital twin

## Design

The previous iteration had a complete delivery-expression-editing chain, but
the editing layer was too coarse for APOE biology. It treated the desired RNA
editing product as one generic on-target state and collected off-target editing
into one aggregate state.

The v4 goal is to make the model answer an engineering design question:

Which PUF-deaminase design, expression level and sampling time should be tested
first to maximize APOE4-to-APOE3-like editing while controlling APOE2-like
double editing and mechanistically distinct off-target risks?

## Build

Added:

- APOE 112/158 multisite editing states and ODEs.
- Editor-type scaling for A3A, APOBEC1, engineered A3A and ProAPOBEC-like
  constructs.
- Three off-target mechanisms: local bystander, PUF-mediated mismatch and
  deaminase-only background.
- Sequence-to-kinetics mapping from PUF design features into `k_on` and `k_cat`
  modifiers.
- A candidate PUF design template for wet-lab planning.
- A v4 design screen that scans candidate design, expression level and sampling
  time.
- A utility score and Pareto-front selector for experiment recommendation.
- Parameter priors and sequence-to-kinetics rules updated from Han et al. 2022
  NAR gkac713 CU-REWIRE design results.

Key files:

- `models/editing/apoe_multisite_editing.py`
- `models/editing/sequence_to_kinetics.py`
- `wetlab/templates/puf_design_candidates_template.csv`
- `scripts/run_iterated_analysis.py`
- `docs/V4_OPTIMIZATION_SUMMARY_CN.md`
- `docs/GKAC713_PARAMETER_UPDATE.md`

## Test

New regression tests check that:

- the APOE multisite model runs through the full simulator;
- v4 metrics include APOE3-like, APOE2-like and separate off-target burdens;
- stronger PUF sequence features increase the expected kinetic modifiers.

The full iterated analysis now also writes:

- `reports/iterated_analysis/v4_design_screen.csv`
- `reports/iterated_analysis/v4_design_recommendation.md`

## Learn

This iteration changes the main story of the model. The AAV-BBB chain remains
available as a translational extension, but the core experimentally testable
model is now a HEK293T-compatible PUF-A3A digital twin:

sequence design -> kinetic parameters -> expression-editing ODE -> mechanistic
off-target classes -> Pareto-ranked next experiments.

The remaining highest-value work is to replace the current heuristic
sequence-to-kinetics coefficients with fitted coefficients from public REWIRE
data and then update them with the team's own wet-lab measurements.
