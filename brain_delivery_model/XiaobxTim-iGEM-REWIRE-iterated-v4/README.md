# brain_delivery_model

A runnable multiscale AAV delivery–expression–editing framework for brain-targeted therapeutic design.

## Iterated REWIRE-style version

This copy is the separated, iterated version of the combined AAV-PUF-APOBEC
model. It keeps the original full-chain ODE simulation and adds iGEM-style
engineering structure:

- wet-lab observation bridge for qPCR, Western blot and amplicon sequencing
- parameter provenance and uncertainty priors
- Monte Carlo uncertainty scan
- REWIRE-style off-target candidate panel aggregation
- transparent random-search parameter fitting against wet-lab observations
- model-informed next-experiment recommendations
- DBTL, decision and failure logs
- v4 APOE 112/158 multisite editing with separate APOE3-like, APOE2-like and mixed-edit outputs
- v4 mechanistic off-target split into local bystander, PUF mismatch and deaminase background classes
- v4 sequence-to-kinetics PUF design screen with Pareto-ranked experiment recommendations

Run the full iterated analysis:

```bash
/opt/anaconda3/bin/python scripts/run_iterated_analysis.py --dose 1.0 --t_end 72 --dt 0.2
```

Outputs are written to:

```text
reports/iterated_analysis/
```

Key new files:

```text
docs/V4_OPTIMIZATION_SUMMARY_CN.md
docs/GKAC713_PARAMETER_UPDATE.md
docs/WIKI_PARAMETER_AND_REFERENCE_SOURCES_CN.md
more_paper/literature_research_for_model_v4_CN.md
docs_wiki/wiki_model_page_full_CN.md
docs_wiki/README.md
models/editing/apoe_multisite_editing.py
models/editing/sequence_to_kinetics.py
models/calibration/wetlab_bridge.py
models/calibration/parameter_fit.py
models/editing/offtarget_panel.py
config/parameter_provenance.yaml
scripts/run_iterated_analysis.py
dbtl/
docs/
wetlab/templates/assay_observations_template.csv
wetlab/templates/offtarget_candidate_panel_template.csv
wetlab/templates/puf_design_candidates_template.csv
```

This project models the full chain from **administration** to **system-level decision metrics**:

1. **Module 1: Absorption**  
   Local depot → lymphatic uptake / direct blood entry
2. **Module 2: Systemic distribution**  
   Blood ↔ liver / peripheral tissues / brain vascular side
3. **Module 3: BBB transport**  
   Brain blood → endothelial surface → endocytosis → brain interstitial space
4. **Module 4: Brain intracellular uptake and expression**  
   Brain ISF → brain cell → nucleus → mRNA → `P_brain`
5. **Module 5: APOE multisite RNA-editing dynamics**  
   `P_brain` drives APOE112/APOE158 editing, APOE3-like/APOE2-like product formation,
   local bystander editing, PUF-mediated off-target editing and deaminase background editing
6. **Module 6: Dose optimization and therapeutic-window analysis**  
   Extracts efficacy/safety metrics and searches for feasible dose regions

---

## Current scope

This version implements a complete forward simulation chain for **Modules 1–5** and a post-processing / optimization layer for **Module 6**.

### Implemented components
- global state vector definition
- route-specific absorption configuration
- module-1 absorption RHS
- module-2 systemic distribution RHS
- module-3 BBB transport RHS
- module-4 intracellular uptake and expression RHS
- module-5 competitive on/off-target editing RHS
- module-6 metric extraction, therapeutic-window evaluation, and dose scan
- plotting utilities for all modules
- runnable `main.py` supporting:
  - single simulation mode
  - optimization mode (`--optimize`)

---

## Project structure

```text
brain_delivery_model/
├── config/
│   └── base_config.yaml
├── models/
│   ├── pbpk/
│   │   ├── lymphatic_absorption.py
│   │   └── organ_distribution.py
│   ├── bbb/
│   │   └── bbb_transport.py
│   ├── intracellular/
│   │   └── aav_intracellular.py
│   ├── editing/
│   │   └── competitive_editing.py
│   └── full_model/
│       ├── state_vector.py
│       ├── rhs_aav.py
│       └── simulator.py
├── optimization/
│   ├── metrics.py
│   ├── therapeutic_window.py
│   ├── min_effective_dose.py
│   └── report_module6.py
├── utils/
│   ├── config_loader.py
│   ├── io.py
│   └── plotting.py
├── main.py
└── outputs/
```

---

## Single simulation

Run one full forward simulation:

```bash
python main.py --route footpad --dose 1.0 --t_end 168 --dt 0.2
```

### Common arguments
- `--route`: `footpad`, `iv`, or `im`
- `--dose`: initial administered AAV amount
- `--t_end`: simulation horizon (hours)
- `--dt`: output sampling step
- `--output_dir`: output folder for figures/results

### Example
```bash
python main.py --route footpad --dose 1.0 --t_end 72 --dt 0.1
```

Outputs are saved under `outputs/run_full_model/` by default.

---

## Optimization mode (Module 6)

Run dose scan and therapeutic-window analysis:

```bash
python main.py --optimize --t_end 168 --dt 0.2
```

### Optional arguments
- `--scan_output_dir`: output directory for dose-scan tables and plots

### Example
```bash
python main.py --optimize --scan_output_dir outputs/run_module6 --t_end 168 --dt 0.2
```

This will:
- scan the configured dose grid from `config/base_config.yaml`
- extract summary metrics such as:
  - final on-target editing rate
  - final off-target burden
  - final specificity index
  - peak brain expression
  - liver AUC
  - blood Cmax
- evaluate feasibility under the configured thresholds
- report the minimum feasible dose if one exists

---

## Key outputs

### Single simulation outputs
Depending on enabled modules, the plotting utilities generate figures such as:
- `module1_states.png`
- `module1_fluxes.png`
- `module2_states.png`
- `module2_fluxes.png`
- `module3_states.png`
- `module3_fluxes.png`
- `module4_states.png`
- `module4_fluxes.png`
- `module5_states.png`
- `module5_fluxes.png`
- `module5_metrics.png`
- `module1_module2_states.png` / combined state plots

### Optimization outputs
Module 6 generates:
- `module6_dose_scan.csv`
- `dose_vs_on_target.png`
- `dose_vs_off_target.png`
- `dose_vs_specificity.png`
- `dose_vs_auc_liver.png`
- `dose_vs_cmax_blood.png`
- `dose_vs_pbrain_peak.png`
- `dose_feasibility_map.png`

---

## Notes on interpretation

- **Modules 1–5** are mechanistic / forward-simulation modules.
- **Module 6** is a decision layer built on top of Modules 1–5.
- Current parameters are suitable for exploration and debugging, but not yet equivalent to calibrated in vivo values.
- Many absorption, BBB, and intracellular efficiency parameters are expected to be **fitted later using wet-lab data**.
- Degradation-type parameters can often be initialized from literature-derived half-lives.

---

## Development roadmap

### Already implemented
- full forward chain from administration to editing
- optimization layer for dose scanning and therapeutic-window evaluation

### Planned / future refinement
- parameter calibration with wet-lab data
- sensitivity analysis
- tighter coupling between expression data and optimization workflow
- optional refinement of off-target pools and disease-effect readouts

---

## Quick start recommendation

For qualitative behavior checks:

```bash
python main.py --route footpad --dose 1.0 --t_end 168 --dt 0.2
```

For decision-layer testing:

```bash
python main.py --optimize --t_end 168 --dt 0.2
```
