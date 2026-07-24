# Figure Captions and Placement Guide

| File | Suggested title | Caption | Best page section |
|---|---|---|---|
| `diagrams/full_model_workflow.svg` | Full modeling workflow | The v4 model connects delivery or transfection input, editor expression, APOE RNA editing, off-target analysis, and design recommendation. | Overview / hero schematic |
| `diagrams/apoe_editing_logic.svg` | APOE112/APOE158 product logic | APOE4, APOE3-like, and APOE2-like states depend on codons 112 and 158. v4 tracks these products separately. | Biological logic |
| `diagrams/dbtl_model_loop.svg` | Model-guided DBTL loop | Each wet-lab result updates the model, and the updated model recommends the next design to test. | Rigor / future iteration |
| `assets/full_chain_state_overview.png` | Full-chain state simulation | Simulation output showing the full delivery-expression-editing chain. | Results |
| `assets/module1_absorption_states.png` | Module 1 absorption | Local dose absorption and drainage dynamics. | Biological process |
| `assets/module2_distribution_states.png` | Module 2 systemic distribution | Distribution among blood, liver, peripheral, and brain vascular compartments. | Biological process |
| `assets/module3_bbb_transport_states.png` | Module 3 BBB transport | Brain vascular, endothelial, endosomal, and interstitial transport states. | Biological process |
| `assets/module4_expression_states.png` | Module 4 expression | Intracellular vector, mRNA, and active editor protein dynamics. | Biological process |
| `assets/apoe_multisite_editing_states.png` | Module 5 APOE editing states | APOE4-like, APOE3-like, APOE2-like, mixed-edit, and off-target state dynamics. | Core editing model |
| `assets/apoe_editing_fluxes.png` | APOE editing fluxes | Binding, unbinding, catalysis, and off-target fluxes in the v4 editing module. | Methods / model details |
| `assets/apoe_editing_metrics.png` | Editing and safety metrics | Time-varying APOE editing fraction, off-target burden, specificity, and risk metrics. | Results |

## Notes for Website Implementation

- Use SVG diagrams for conceptual explanation.
- Use PNG simulation plots for evidence/results.
- Keep captions close to figures so readers know what each plot proves.
- Avoid presenting literature-informed priors as final measured constants.
