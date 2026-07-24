# Parameter Sources, Rigor, and Validation Logic

## Literature-Grounded Priors

The model does not invent its core biological rules. We used reliable papers to
initialize the model:

- CU-REWIRE paper: PUF-APOBEC3A architecture, 48 h HEK293T reference, +2 editing
  window, UC context, 10R PUF specificity, deaminase background.
- PUF engineering papers: repeat-number effect and RNA recognition code.
- APOE biology papers: APOE112/APOE158 isoform definitions and APOE2-like LDLR
  binding risk.
- mRNA/protein turnover papers: broad degradation priors.
- AAV/PBPK papers: optional translational delivery priors.
- systems biology papers: identifiability, uncertainty, and model-guided
  experiment design.

Detailed source table:

```text
docs/WIKI_PARAMETER_AND_REFERENCE_SOURCES_CN.md
```

## Why This Is Rigorous

The model is rigorous because it avoids three common mistakes.

First, it does not collapse APOE biology into one on-target state. It separates
APOE3-like and APOE2-like products.

Second, it does not treat all off-targets as the same. It separates local
bystander, PUF mismatch, and deaminase background risks.

Third, it does not pretend that literature priors are final truths. Parameters
are documented with uncertainty ranges and are designed to be updated by wet-lab
data.

## Current Validation Strategy

The v4 code includes tests for:

- APOE multisite simulation;
- sequence-to-kinetics behavior;
- off-target panel integration;
- wet-lab bridge calculations.

The current analysis also produces:

- baseline simulation;
- fitted placeholder run;
- uncertainty samples;
- design screen;
- design recommendation.

## What Wet Lab Should Measure Next

| Priority | Data | Model parameter updated |
|---:|---|---|
| 1 | APOE112 editing at 24/48/72 h | `k_cat_112`, `k_on_112` |
| 2 | APOE158 editing at 24/48/72 h | `k_cat_158`, `k_on_158` |
| 3 | linked APOE112/APOE158 amplicon | APOE3-like vs APOE2-like split |
| 4 | nearby C bystander edits | `local_bystander_per_112/158` |
| 5 | top-K off-target amplicons | `k_on_puf_off`, `k_cat_puf_off` |
| 6 | mock / PUF-only / inactive A3A | `k_deaminase_bg` and background noise |
| 7 | qPCR and Western time courses | `k_tx`, `k_tl`, `k_deg_m`, `k_deg_p` |

## Honest Limitation Statement

Most kinetic constants are normalized effective parameters. Without rich
time-course data, the model may identify an effective editing strength such as
`k_cat / K_D` better than individual `k_on`, `k_off`, and `k_cat`. Future
iterations should add profile likelihood or parameter recovery before claiming
precise parameter estimates.
