# v4 Wiki Content Package

> This folder is prepared for the iGEM wiki team. It organizes the v4 modeling
> work into readable page sections, diagrams, figure assets, captions, and
> future-iteration notes.

## Recommended Wiki Page Order

Single assembled page:

- `wiki_model_page_full_CN.md`

Section-by-section source files:

1. `sections/01_model_overview_CN.md`
2. `sections/02_biological_process_CN.md`
3. `sections/03_model_architecture_CN.md`
4. `sections/04_parameters_and_rigor_CN.md`
5. `sections/05_results_and_conclusions_CN.md`
6. `sections/06_future_iterations_CN.md`

## Visual Assets

| Asset | Suggested use |
|---|---|
| `diagrams/full_model_workflow.svg` | Hero or first schematic: what the model does end to end. |
| `diagrams/apoe_editing_logic.svg` | Explain APOE112/APOE158 product logic. |
| `diagrams/dbtl_model_loop.svg` | Explain model-guided DBTL loop. |
| `assets/full_chain_state_overview.png` | Show full-chain simulation output. |
| `assets/module1_absorption_states.png` | Module 1 figure: absorption. |
| `assets/module2_distribution_states.png` | Module 2 figure: systemic distribution. |
| `assets/module3_bbb_transport_states.png` | Module 3 figure: BBB transport. |
| `assets/module4_expression_states.png` | Module 4 figure: intracellular expression. |
| `assets/apoe_multisite_editing_states.png` | Module 5 figure: APOE multisite editing states. |
| `assets/apoe_editing_fluxes.png` | Module 5 figure: editing fluxes. |
| `assets/apoe_editing_metrics.png` | Module 5 figure: editing/safety metrics. |

## One-Sentence Model Message

We built a sequence-aware, literature-grounded digital twin that connects
PUF-APOBEC design to expression, APOE112/APOE158 RNA editing, off-target risk,
and model-guided experimental selection.

## Notes for the Website Team

- The markdown files are written as direct wiki content, not internal notes.
- The SVG diagrams are lightweight and can be embedded directly.
- The PNG plots are generated from v4 simulation outputs.
- The model still uses literature-informed priors; the wiki should emphasize
  that future wet-lab data will update the model.
