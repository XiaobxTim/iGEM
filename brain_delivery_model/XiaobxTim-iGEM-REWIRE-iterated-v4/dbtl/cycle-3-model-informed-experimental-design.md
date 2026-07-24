# Cycle 3: Model-informed experimental design

## Design

The model should actively suggest the next wet-lab experiments. This mirrors
good iGEM practice: models guide experiments, experiments update models.

## Build

The iterated analysis script ranks uncertain parameters by correlation with key
decision metrics and writes:

```text
reports/iterated_analysis/experimental_design_recommendations.md
```

## Test

Each recommendation maps a sensitive parameter to a real experiment, such as:

- brain vector qPCR for BBB delivery uncertainty
- PUF-APOBEC qPCR/Western for expression uncertainty
- amplicon sequencing time course for editing kinetics
- off-target panel sequencing for safety uncertainty

## Learn

The next wet-lab round should not measure everything equally. It should measure
the parameters that most change dose feasibility and specificity.
