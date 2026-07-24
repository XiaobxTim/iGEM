# Cycle 4: REWIRE-style off-target panel integration

## Design

The previous model had a single aggregate `S_off` pool. That was useful for
ODE simplicity but weak for an RNA-editing project, because off-target burden
should eventually come from transcriptome-wide candidate sites.

## Build

Added:

- `models/editing/offtarget_panel.py`
- `wetlab/templates/offtarget_candidate_panel_template.csv`
- `reports/iterated_analysis/offtarget_panel_summary.csv`

The panel compresses candidate sites into:

- effective `S_off_init`
- panel-scaled `k_on_off`
- panel-scaled `k_cat_off`

## Test

Added `tests/test_offtarget_panel.py` to verify accessibility-weighted
aggregation.

## Learn

This gives us a clean bridge to the REWIRE pipeline. Once the team has a real
candidate table, we can replace the placeholder CSV without changing the ODE
core.
