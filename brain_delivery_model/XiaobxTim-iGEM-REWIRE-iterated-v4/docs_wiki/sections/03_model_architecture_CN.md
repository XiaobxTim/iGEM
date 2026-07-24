# Model Architecture

## Overall Structure

The model is an ODE-based multiscale simulation. Each module contributes a set
of state variables and fluxes, and the simulator integrates them over time.

![Full-chain state overview](../assets/full_chain_state_overview.png)

## State Vector Design

The global state vector allows all modules to run in one simulation:

- delivery states;
- blood/liver/brain distribution states;
- BBB transport states;
- intracellular vector/mRNA/protein states;
- APOE112/APOE158 editing states;
- off-target states;
- compatibility accumulators for old reporting code.

This architecture lets us keep earlier full-chain delivery modeling while
upgrading the RNA editing layer.

## Sequence-to-Kinetics Layer

The v4 model adds a new layer between sequence design and ODE parameters.

Candidate design features:

- PUF repeat number;
- PUF match score;
- mismatch count;
- RNA accessibility;
- editable C distance from PUF binding site;
- UC context score;
- editor type;
- editor activity scale.

These are converted into effective kinetic modifiers:

```text
k_on scale = exp(beta1 * PUF_score + beta2 * accessibility)
k_cat scale = (0.04 + 0.96 * UC_context) * exp(-0.9 * abs(distance - 2))
PUF off-target scale = exp(-1.4 * max(PUF_repeats - 8, 0)) * exp(0.45 * mismatch_count)
```

This means the model can compare biological designs, not only numerical
parameter sets.

## Off-Target Decomposition

Instead of one generic off-target number, v4 separates:

| Off-target class | Biological meaning | Wet-lab validation |
|---|---|---|
| Local bystander | PUF binds the correct target, but nearby C sites are edited | APOE amplicon sequencing around target window |
| PUF mismatch | PUF binds similar RNA sequences elsewhere | Top-K predicted off-target amplicon panel / RNA-seq |
| Deaminase background | APOBEC edits accessible RNA independent of correct PUF targeting | mock, PUF-only, inactive A3A, free-A3A controls |

This is important because each class requires a different engineering fix.

## Parameter Provenance

Parameter priors are documented in:

- `config/parameter_provenance.yaml`
- `docs/WIKI_PARAMETER_AND_REFERENCE_SOURCES_CN.md`
- `docs/GKAC713_PARAMETER_UPDATE.md`

The model clearly marks which values are literature-informed priors and which
need wet-lab calibration.

## Reproducibility

The whole v4 analysis can be regenerated with:

```bash
/opt/anaconda3/bin/python scripts/run_iterated_analysis.py --dose 1.0 --t_end 48 --dt 0.5 --fit_samples 20 --uncertainty_samples 20
```

The generated outputs are placed in:

```text
reports/iterated_analysis/
```
