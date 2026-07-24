# Cycle 5: Wet-lab parameter fitting loop

## Design

The previous iteration compared predictions to assay observations but did not
update parameters. Strong iGEM modeling should close the loop: wet-lab data
should change the next model.

## Build

Added:

- `models/calibration/parameter_fit.py`
- `reports/iterated_analysis/parameter_fit_results.csv`
- `reports/iterated_analysis/best_fit_config.yaml`
- fitted trajectory plots under `reports/iterated_analysis/fitted/`

The fitter uses transparent random search over high-impact parameters. This is
simple enough to explain on the wiki and robust enough for early sparse data.

## Test

The full iterated analysis now runs baseline simulation, residual calculation,
off-target panel aggregation, random-search fitting and uncertainty analysis in
one command.

## Learn

The current observations are placeholders, so the fitted values should not be
treated as biology. The value of this cycle is the workflow: when real qPCR,
Western and amplicon data arrive, the same command will produce a fitted model.
