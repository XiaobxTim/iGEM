# Cycle 2: Parameter provenance and uncertainty

## Design

Strong iGEM models explain where parameters come from and which are still
assumptions. Our previous combined model had many effective rates but no
machine-readable evidence table.

## Build

Added:

- `config/parameter_provenance.yaml`
- `optimization/iterated/uncertainty.py`
- Monte Carlo uncertainty scan in `scripts/run_iterated_analysis.py`

## Test

The analysis script samples uncertain parameters, runs the full ODE chain and
saves `uncertainty_samples.csv`.

## Learn

The model can now show fragility. If on-target editing varies mainly with
`editing.k_cat_on`, wet lab should prioritize editing time-course measurement.
If safety varies mainly with liver exposure, wet lab should prioritize vector
biodistribution assays.
