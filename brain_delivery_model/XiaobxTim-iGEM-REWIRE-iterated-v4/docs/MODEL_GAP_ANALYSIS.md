# Gap analysis against strong iGEM modeling projects

## Current strengths

- The project already has a full ODE chain from AAV administration to editing
  outcome.
- Module boundaries are interpretable: absorption, PBPK-like distribution, BBB,
  intracellular expression, PUF-APOBEC editing and dose optimization.
- The explicit PUF-APOBEC module tracks binding, unbinding, catalysis and
  off-target products.

## Main gaps before this iteration

| gap | why it matters | iteration added |
|---|---|---|
| No wet-lab data interface | The model could not be fitted or falsified by qPCR, Western or amplicon data | `models/calibration/wetlab_bridge.py` and assay CSV template |
| Parameters lacked provenance | Judges and wet-lab teammates could not tell which values were assumptions | `config/parameter_provenance.yaml` |
| No uncertainty analysis | A single trajectory hides fragile assumptions | Monte Carlo scan in `scripts/run_iterated_analysis.py` |
| No model-informed next experiments | The model did not tell wet lab what to measure next | `experimental_design_recommendations.md` output |
| No DBTL/decision/failure record | The model looked like a one-shot result rather than engineering | `dbtl/` logs |
| Weak wiki-ready explanation | Equations existed, but not the project story | `docs/` reports and updated README |
| Aggregate off-target pool | RNA editing risk should be grounded in candidate sites | REWIRE-style off-target candidate panel aggregation |
| No parameter update loop | Residuals alone do not close the DBTL loop | Random-search fitting over high-impact parameters |

## Remaining gaps for future real wet-lab integration

1. Replace placeholder assay observations with real measurements.
2. Replace placeholder off-target candidate panel with real REWIRE output.
3. Add replicate-level measurement noise rather than one SD per assay.
4. Promote top off-target sites to explicit ODE states if they dominate burden.
5. Separate brain cell populations if wet-lab data distinguishes neurons,
   astrocytes or endothelial cells.
6. Validate dose scaling experimentally before interpreting absolute dose.
