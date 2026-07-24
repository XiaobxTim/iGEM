# Reproducibility checklist

- [x] Single command analysis entry point: `scripts/run_iterated_analysis.py`
- [x] Parameter priors are machine-readable: `config/parameter_provenance.yaml`
- [x] Wet-lab assay template is versioned: `wetlab/templates/assay_observations_template.csv`
- [x] Baseline outputs are saved under `reports/iterated_analysis/`
- [x] Uncertainty samples are saved as CSV
- [x] Off-target candidate panel is versioned as a CSV template
- [x] Parameter fitting writes all tried parameter sets and the best-fit config
- [x] DBTL records are stored in `dbtl/`
- [ ] Replace placeholder assay values with measured wet-lab data
- [ ] Add replicate-aware parameter fitting
- [ ] Add sequence-specific off-target panel from REWIRE pipeline outputs
