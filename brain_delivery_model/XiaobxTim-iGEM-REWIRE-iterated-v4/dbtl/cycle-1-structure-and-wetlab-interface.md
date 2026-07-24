# Cycle 1: Structure and wet-lab interface

## Design

Adopt a REWIRE-like project structure so the model can be used by wet lab and
wiki writers. The key missing interface was a way to compare ODE states with
real assays.

## Build

Added:

- `models/calibration/wetlab_bridge.py`
- `wetlab/templates/assay_observations_template.csv`
- `docs/IGEM_BENCHMARK_RESEARCH.md`
- `docs/MODEL_GAP_ANALYSIS.md`

## Test

The bridge maps model states to qPCR, Western blot, amplicon editing and
viability-style observations. It can compute residuals and a weighted RMSE once
real wet-lab data are available.

## Learn

The model now has a clean path for wet-lab calibration. The next bottleneck is
parameter credibility, so the next cycle should add provenance and uncertainty.
