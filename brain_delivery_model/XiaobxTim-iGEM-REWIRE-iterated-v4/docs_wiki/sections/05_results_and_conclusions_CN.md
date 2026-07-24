# Simulation Results and Conclusions

## Key Result 1: The Full Chain Can Be Simulated End to End

The model runs from delivery input through expression and RNA editing. This
makes it possible to ask whether a design fails because of delivery,
expression, binding, catalysis, or off-target risk.

![Full-chain state overview](../assets/full_chain_state_overview.png)

Conclusion: the model is not a single isolated editing equation. It connects
upstream exposure to downstream editing outcomes.

## Key Result 2: APOE3-like and APOE2-like Products Are Separated

The v4 RNA editing module tracks APOE112 and APOE158 separately.

![APOE multisite editing states](../assets/apoe_multisite_editing_states.png)

Conclusion: this makes the model biologically more accurate. It can distinguish
the intended APOE4-to-APOE3-like direction from APOE2-like double editing risk.

## Key Result 3: Off-Target Burden Is Mechanistically Interpretable

The model separates local bystander, PUF mismatch, and deaminase-only
background. The metrics plot summarizes target editing, specificity, and
off-target burden over time.

![APOE editing metrics](../assets/apoe_editing_metrics.png)

Conclusion: if off-target burden is high, we can tell which engineering lever
to change. For example:

- high local bystander -> move target C/window or change deaminase;
- high PUF mismatch -> redesign PUF target site or use 10R/ePUF;
- high deaminase background -> engineer deaminase or lower expression.

## Key Result 4: v4 Can Recommend Designs

The current design screen compares PUF-deaminase candidate, dose level, and
sampling time. With the current literature-informed priors, the top result in
`reports/iterated_analysis/v4_design_recommendation.md` is:

| Output | Current v4 recommendation |
|---|---|
| Design | `10R-proapobec` |
| Editor | `ProAPOBEC` |
| Expression level | `high` |
| Sampling time | `72 h` |
| Expected APOE3-like fraction | about `0.0666` |
| APOE2-like risk proxy | about `3.6e-06` |
| Off-target burden | about `0.00373` |

Important interpretation: 72 h appears because the current simulator still
contains the full delivery-expression chain. For a pure HEK293T transfection
model, 48 h remains the key literature-based assay point from CU-REWIRE.

## Key Result 5: The Model Gives a Testable Next Step

The model does not end with a curve. It proposes what data would most improve
the next cycle:

- APOE112 and APOE158 amplicon sequencing at 24/48/72 h;
- linked amplicon or long-read data to separate APOE3-like and APOE2-like;
- local bystander analysis;
- top-K off-target amplicon panel;
- mock, PUF-only, inactive A3A controls;
- editor mRNA and protein time courses.

## Main Conclusion for Wiki

The v4 model turns our project from trial-and-error design into a measurable
engineering loop. It uses literature-derived priors to make initial predictions,
then identifies the wet-lab data needed to replace those priors with our own
measurements.
