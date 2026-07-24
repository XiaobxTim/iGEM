# Model Overview

## What We Wanted to Solve

Our project aims to evaluate whether a programmable PUF-APOBEC RNA editor can
convert APOE4-like RNA toward a safer APOE3-like state while keeping off-target
editing under control. The key modeling question is not only "will editing
happen?", but "which design should we test first?"

The v4 model therefore works as a digital twin for design selection. It links:

1. delivery or transfection input;
2. intracellular editor expression;
3. PUF-guided RNA binding;
4. APOE112/APOE158 C-to-U editing;
5. local and transcriptome-wide off-target risk;
6. design ranking and next-experiment recommendation.

![Full model workflow](../diagrams/full_model_workflow.svg)

## Why the Model Matters

Without modeling, wet-lab design choices can become trial-and-error: choose one
PUF, choose one dose, choose one time point, then hope the result is useful.
Our model makes the choice explicit. It predicts what each candidate design is
expected to optimize and what risk it may introduce.

The model helps us decide:

- whether a design mainly drives APOE4-to-APOE3-like editing;
- whether it also creates unwanted APOE2-like double editing;
- whether off-target risk is local bystander, PUF mismatch, or deaminase
  background;
- whether weak editing is more likely caused by poor expression, poor binding,
  or poor catalysis;
- which time point and expression level should be tested next.

## Main v4 Upgrade Compared with Earlier Versions

Earlier versions treated on-target editing as one generic product. v4 separates
APOE112 and APOE158 because APOE biology depends on both sites:

- APOE4: Arg112 / Arg158
- APOE3-like: Cys112 / Arg158
- APOE2-like: Cys112 / Cys158

This prevents the model from incorrectly calling all successful editing
"APOE2-like" and lets us track therapeutic benefit and risk separately.

![APOE editing logic](../diagrams/apoe_editing_logic.svg)

## Website Takeaway

The modeling contribution is a decision-making system: it translates PUF design
features into kinetic parameters, simulates editing and risk, and recommends the
next experimental condition to test.
