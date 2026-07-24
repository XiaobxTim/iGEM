# Model Parameter and Reference Sources

> Wiki-ready version. This page explains which online/literature-derived
> parameters or rules are used in our v4 model, where they come from, how they
> are used, and what still needs wet-lab calibration.

## 1. How to Read This Page

Our model uses normalized ODE units. Many values are therefore **effective
priors**, not directly measured biochemical constants. We classify sources into
three levels:

- **Direct literature value/rule**: the paper reports a value, window, motif,
  or experimental condition that can be used almost directly.
- **Literature-informed prior**: the paper supports the range or direction of
  a parameter, but the exact number is adapted to our normalized model.
- **Model assumption / placeholder**: the parameter is necessary to run the
  model but must be replaced by our own wet-lab data.

The most important principle is honesty: these literature references make the
model scientifically grounded, but they do not replace our own APOE wet-lab
measurements.

---

## 2. Parameter Sources Used in the Current v4 Model

| Model part | Parameter(s) in v4 | Current v4 value or rule | What it means | Source type | Literature / URL | How we use it | Caveat / what wet lab should update |
|---|---|---:|---|---|---|---|---|
| Simulation time | `simulation.default_t_end` | `48.0 h` | Default HEK293T-style readout time | Direct experimental condition | Han et al., 2022, NAR. https://doi.org/10.1093/nar/gkac713 | CU-REWIRE assays commonly use 48 h after transfection, so v4 uses 48 h as the default cellular reference point. | Our simulator still contains AAV/BBB delivery modules; full-chain predictions can shift later, e.g. 72 h. HEK293T-only simulations should bypass delivery modules. |
| Editor identity | `editing.editor_type` | `A3A` | Default C-to-U deaminase identity | Direct system match | Han et al., 2022, NAR. https://doi.org/10.1093/nar/gkac713 | CU-REWIRE uses APOBEC3A/A3A for C-to-U editing, so v4 does not mix generic APOBEC, APOBEC1, and A3A assumptions. | If we change enzyme, `editor_type` must change to `APOBEC1`, `engineered_A3A`, or `ProAPOBEC`. |
| Inactive editor control | `editor_activity_scale` in `puf_design_candidates_template.csv` | active: `1.0`; inactive/PUF-only: `0.0` | Turns catalytic editing on/off for controls | Direct experimental logic | Han et al., 2022 reports inactive A3A E72A abolishes detectable CU-REWIRE editing. https://doi.org/10.1093/nar/gkac713 | Used to model inactive deaminase and PUF-only controls as binding/expression controls without catalytic activity. | Actual control plasmids still need wet-lab confirmation for expression and localization. |
| APOE codon logic | `S_APOE4`, `S_APOE3_like`, `S_APOE2_like`, `S_APOE158_only` | state definitions | Separates APOE4, APOE3-like, APOE2-like, and incomplete editing | Direct biological rule | Hatters et al., 2006; Huang & Mahley, 2014; Bu, 2009. URLs in reference list. | Editing codon 112 only is interpreted as APOE3-like; editing both 112 and 158 is APOE2-like. | RNA editing does not guarantee protein-level functional rescue; protein and lipid-binding assays are still needed later. |
| APOE2-like risk | `ldlr_binding_risk_proxy` | `APOE2_like_fraction` | Risk proxy for codon-158 editing | Direct biology, model proxy | Hatters et al., 2006; Huang & Mahley, 2014. | APOE2 has poor LDLR binding, so v4 penalizes unwanted double editing if the therapeutic goal is APOE4-to-APOE3-like. | Whether APOE2-like is harmful depends on project goal and disease context; it is a proxy, not a full lipid-metabolism model. |
| APOE112 catalysis | `editing.k_cat_112` | `1.0` | Effective catalytic turnover at target codon 112 | Literature-informed prior | Han et al., 2022 reports CU-REWIRE C-to-U editing can reach high editing in preferred contexts. https://doi.org/10.1093/nar/gkac713 | Set higher than APOE158 because APOE112 is the desired APOE4-to-APOE3-like edit in our current design. | Normalized value, not a real s^-1 or h^-1 biochemical measurement. Must fit with APOE112 amplicon time course. |
| APOE158 catalysis | `editing.k_cat_158` | `0.15` | Effective catalytic turnover at codon 158 | Literature-informed prior + therapeutic choice | Han et al., 2022 editing-window logic; APOE isoform biology from Hatters et al., 2006. | Lower default because 158 is treated as double-edit/APOE2-like risk unless dual-site editing is chosen. | If the team chooses APOE2-like as the goal, this prior and utility penalty must be changed. |
| PUF binding to APOE112 | `editing.k_on_112` | `1.0`, scaled by sequence-to-kinetics | Effective binding/capture of target APOE112 RNA | Literature-informed prior | PUF recognition code papers: Cheong & Hall, 2006; Campbell et al., 2014; Zhao et al., 2018. | Sequence features modify `k_on_112` through PUF score, repeat number, mismatch count, and accessibility. | Actual `k_on` cannot be separated from `k_cat` using only endpoint editing data; needs EMSA/SPR or rich time-course fitting. |
| PUF binding to APOE158 | `editing.k_on_158` | `0.35`, scaled by sequence-to-kinetics | Effective binding/capture of APOE158 window | Literature-informed prior | Same PUF recognition literature as above. | Lower than `k_on_112` because current candidate table prioritizes APOE112 targeting. | Must be recalculated from the real APOE158 target-site design. |
| Editor-RNA dissociation | `editing.k_off_112`, `editing.k_off_158` | `0.5`, `0.7` | Effective dissociation from target RNA | Model assumption anchored by PUF binding literature | Cheong & Hall, 2006; Zhao et al., 2018. | Needed for explicit `E + S <-> C -> E + product` kinetics. | These are not directly measured. If we cannot identify them, wiki should report effective `k_cat/KD` instead. |
| Editing window | `distance_112`, `distance_158`; `kcat_scale = (0.04 + 0.96 * UC_context_score) * exp(-0.9 * abs(distance - 2))` | best near `distance = 2` | Editing C position relative to PUF binding site | Direct rule converted to model function | Han et al., 2022, NAR. https://doi.org/10.1093/nar/gkac713 | Converts candidate PUF design into `k_cat` modifier. Designs placing the editable C near the preferred window receive higher catalytic rate. | The exact exponential slope is our normalized prior; it should be fitted using REWIRE public data and our APOE data. |
| Sequence context | `uc_context_112`, `uc_context_158`, `uc_context_scale` | 0-1 context score; A3A default `1.0` | UC-like context increases cytidine editing | Direct rule, normalized score | Han et al., 2022; Huang et al., 2020. | Used to increase `k_cat` for UC/UCN-like target sites and reduce it for poor contexts. | Context score is currently hand-curated; next version should compute it from real APOE sequence automatically. |
| RNA structure/accessibility | `accessibility_112`, `accessibility_158` in design table | 0-1 score | Exposed/single-stranded RNA is more editable/bindable | Literature-informed prior | Huang et al., 2020; RNA-structure prediction should be added with RNAfold or SHAPE-style data. | Used in sequence-to-kinetics to increase `k_on` when target RNA is predicted accessible. | Current scores are placeholders; should come from RNAfold, SHAPE-MaP, or target amplicon/structure analysis. |
| PUF repeat number | `puf_repeats` | 8 or 10 in template | Recognition length/specificity | Direct design rule | Zhao et al., 2018; Han et al., 2022. | 10R designs receive lower PUF-mismatch off-target prior than 8R designs. | More repeats are not always better; PUF9/10 often peak in affinity, and very long PUFs may have structural/expression burdens. |
| PUF mismatch off-target | `puf_offtarget_scale = exp(-1.4 * max(repeats - 8, 0)) * exp(0.45 * mismatch_count)` | computed per design | PUF-mediated binding to similar transcript sites | Literature-informed prior | Han et al., 2022 reports optimization reducing global off-target; Zhao et al., 2018 supports repeat-number specificity effects. | Rescales `editing.k_on_puf_off` and `editing.k_cat_puf_off`. | The exact coefficients are model priors, not measured constants. Need transcriptome search + RNA-seq/off-target amplicons. |
| PUF-mediated off-target binding/catalysis | `editing.k_on_puf_off`, `editing.k_cat_puf_off` | `0.015`, `0.006` before design scaling | Aggregate off-target caused by PUF binding to similar sequences | Literature-informed prior | Han et al., 2022; REWIRE public dataset GSE155734. | Forms a separate ODE branch `S_puf_off -> C_puf_off -> E_puf_off`. | Needs top-K predicted off-target panel and RNA-seq to calibrate. |
| Deaminase background editing | `editing.k_deaminase_bg` | `0.0003` before editor/design scaling | Off-target editing caused by A3A activity independent of correct PUF targeting | Direct mechanism, normalized prior | Han et al., 2022 reports APOBEC3A alone causes broad RNA off-target; Huang et al., 2020 also discusses A3A off-target. | Separate ODE branch `S_deaminase_bg -> E_deaminase_bg`; reduced/zeroed for inactive controls. | Exact rate depends on expression, localization, cell type, and deaminase variant. Must use mock, PUF-only, inactive-A3A, and free-A3A controls. |
| Local bystander editing | `local_bystander_per_112`, `local_bystander_per_158` | `0.04`, `0.06` | Nearby C edits after correct target engagement | Literature-informed prior | Han et al., 2022; Huang et al., 2020. | Separate accumulator `B_local_bystander` so local bystander is not mixed with transcriptome off-target. | Needs amplicon sequencing of nearby Cs around APOE112/158. |
| ProAPOBEC scaling | `proapobec_cat_scale`, `proapobec_background_scale` | `1.25`, `0.5` | Engineered deaminase expected to improve activity and/or reduce background | Literature-informed prior | Han et al., 2025, Nature Communications. https://doi.org/10.1038/s41467-025-64748-6 | Allows candidate `10R-proapobec` to be compared with A3A designs. | 2025 system is newer and variant-specific; our values are placeholders until we choose a real ProAPOBEC sequence. |
| APOBEC1 scaling | `apobec1_cat_scale`, `apobec1_background_scale` | `0.7`, `0.4` | Allows non-A3A editor comparison without mixing assumptions | Model assumption based on APOBEC-family differences | Huang et al., 2020; APOBEC/base-editor off-target literature. | Used only if `editor_type=APOBEC1` or APOBEC1-like. | We do not currently have APOBEC1-specific REWIRE data for APOE; should not be highlighted as final. |
| mRNA degradation | `intracellular.k_deg_m`, `editing.k_deg_apoe` | `0.08`, `0.1` | mRNA turnover | Literature-informed prior | Wu et al., 2020 reports HEK293T median mRNA half-life around 5.5 h. https://pmc.ncbi.nlm.nih.gov/articles/PMC7175251/ | Gives broad rate prior. For half-life `t1/2`, `k = ln(2)/t1/2`; 5.5 h corresponds to about 0.126 h^-1. | Editor mRNA and APOE mRNA may differ. Need qPCR time courses and mRNA stability assays. |
| Protein degradation | `intracellular.k_deg_p` | `0.03` | Editor protein turnover | Literature-informed prior | Schwanhausser et al., 2011; Boisvert et al., 2012; newly synthesized protein half-life studies. | Broad prior covering stable fusion proteins; `0.03 h^-1` corresponds to half-life about 23 h. | Fusion protein stability depends on linker, tags, localization, degradation signals. Need Western/fluorescence time course. |
| Transcription / translation | `intracellular.k_tx`, `intracellular.k_tl` | `0.4`, `0.9` | Effective expression from delivered vector/plasmid | Model assumption bounded by gene-expression literature | Schwanhausser et al., 2011; general mammalian expression dynamics. | Needed to connect delivery/plasmid input to active editor abundance. | Must be fitted using editor mRNA qPCR and protein Western/fluorescence. |
| AAV liver distribution | `distribution.k_blood_to_liver` | `0.08` | Liver as major systemic AAV sink | Literature-informed prior | AAV PBPK paper: https://doi.org/10.1016/j.xphs.2023.10.005 | Keeps liver exposure penalty in translational extension. | Not calibrated to our capsid/route/dose; use liver vector-genome qPCR if doing animal work. |
| AAV brain entry | `distribution.k_blood_to_brain`, BBB rates | `0.01`, `0.02`, `0.01`, etc. | Low but nonzero CNS entry / BBB transport | Literature-informed prior | Duque et al., 2009; AAV9/PBPK and BBB-capsid literature. | Keeps the AAV-BBB extension biologically plausible. | Highly capsid-, species-, age-, route-, and receptor-dependent. Not core for HEK293T model. |
| Therapeutic-window thresholds | `E_on_target_threshold`, `E_off_max`, `SI_min`, `AUC_liver_max`, `Cmax_blood_max` | `0.2`, `0.05`, `5.0`, `5.0`, `1.0` | Decision thresholds for dose/design screening | Model assumption | Inspired by iGEM model-design practice and safety/efficacy tradeoff logic. | Used to rank feasible regions and compare design options. | These are project-level choices, not biological constants. Team should revise after defining acceptable efficacy/safety standards. |
| Utility score weights | `utility()` in `scripts/run_iterated_analysis.py` | weights for APOE3-like benefit and off-target penalties | Converts outputs into experiment recommendation | Model-design assumption informed by optimal experimental design literature | Raue et al., 2009; Jeong et al., 2018; Banga & Balsa-Canto, 2008. | Makes the model guide which design/time/dose to test next. | Weights should be transparent and sensitivity-tested; different therapeutic goals imply different weights. |

---

## 3. Modeling Design References and How We Borrowed Them

| Design idea in our model | Reference in APA format | URL / DOI | How we borrowed it | Limitation |
|---|---|---|---|---|
| REWIRE as the core biological template | Han, W., Huang, W., Wei, T., Ye, Y., Mao, M., & Wang, Z. (2022). Programmable RNA base editing with a single gRNA-free enzyme. *Nucleic Acids Research, 50*(16), 9580-9595. | https://doi.org/10.1093/nar/gkac713 | This is the main reason our v4 core is PUF + APOBEC3A and why we model editing as sequence-guided C-to-U RNA editing without guide RNA. | The tested targets are not APOE, so editing rates are used as priors, not final APOE values. |
| Public benchmark data for REWIRE | Han et al. (2022) sequencing datasets. GEO accession GSE155734; National Omics Data Encyclopedia OEP001013. | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE155734 | Gives a future path to fit sequence-to-kinetics coefficients from public on-target/off-target data. | We have not yet downloaded/reprocessed these data into model coefficients. |
| APOBEC3A editing motif and RNA-structure dependence | Huang, X., Lv, J., Li, Y., Mao, S., Li, Z., Jing, Z., Sun, Y., Zhang, X., Shen, S., Wang, X., Di, M., Ge, J., Huang, X., Zuo, E., & Chi, T. (2020). Programmable C-to-U RNA editing using the human APOBEC3A deaminase. *The EMBO Journal, 39*(22), e104741. | https://doi.org/10.15252/embj.2020104741 | Supports the UC-context and accessibility terms in `sequence_to_kinetics.py`, and supports separating A3A background editing from PUF-mediated targeting. | CURE is Cas13-guided, not PUF-guided, so its spatial window differs from REWIRE. |
| ProAPOBEC / engineered deaminase extension | Han, W., Yuan, B., Fan, X., Li, W., Yuan, Y., Zhang, Y., Wang, S., Shan, S., Hafner, M., Wang, Z., & Qiu, Z. (2025). Effective in vivo RNA base editing via engineered cytidine deaminase APOBECs fused with PUF proteins. *Nature Communications, 16*, 9727. | https://doi.org/10.1038/s41467-025-64748-6 | Justifies `editor_type=ProAPOBEC`, `proapobec_cat_scale`, and lower background prior for engineered deaminase candidates. | Published after the original REWIRE work and variant-specific; do not overclaim until the exact ProAPOBEC is chosen. |
| PUF repeat-number effect | Zhao, Y.-Y., Mao, M.-W., Zhang, W.-J., Wang, J., Li, H.-T., Yang, Y., Wang, Z., & Wu, J.-W. (2018). Expanding RNA binding specificity and affinity of engineered PUF domains. *Nucleic Acids Research, 46*(9), 4771-4782. | https://doi.org/10.1093/nar/gky134 | Supports treating 8R/10R as design variables and using 10R as a specificity-favoring default. | Binding affinity is not identical to final editing efficiency. |
| PUF recognition code | Cheong, C.-G., & Hall, T. M. T. (2006). Engineering RNA sequence specificity of Pumilio repeats. *Proceedings of the National Academy of Sciences, 103*(37), 13635-13639. | https://doi.org/10.1073/pnas.0606294103 | Supports computing PUF-target matching from repeat-level recognition rather than treating PUF score as arbitrary. | Early engineering study; does not include deaminase fusion or transcriptome off-target. |
| PUF specificity code for endogenous RNA control | Campbell, Z. T., Valley, C. T., & Wickens, M. (2014). A protein-RNA specificity code enables targeted activation of an endogenous human transcript. *Nature Structural & Molecular Biology, 21*(8), 732-738. | https://doi.org/10.1038/nsmb.2847 | Supports the idea that engineered PUF domains can rationally target endogenous mRNA sequence motifs in human cells. | Translational activation system, not RNA base editing. |
| APOE isoform definitions | Hatters, D. M., Peters-Libeu, C. A., & Weisgraber, K. H. (2006). Apolipoprotein E structure: Insights into function. *Trends in Biochemical Sciences, 31*(8), 445-454. | https://doi.org/10.1016/j.tibs.2006.06.008 | Supports the v4 state split: APOE4 = Arg112/Arg158; APOE3-like = Cys112/Arg158; APOE2-like = Cys112/Cys158. | Protein-structure review, not RNA editing experiment. |
| APOE receptor-binding risk | Huang, Y., & Mahley, R. W. (2014). Apolipoprotein E: Structure and function in lipid metabolism, neurobiology, and Alzheimer's diseases. *Neurobiology of Disease, 72*, 3-12. | https://doi.org/10.1016/j.nbd.2014.08.025 | Supports using APOE2-like double editing as an LDLR-binding risk proxy when the goal is APOE4-to-APOE3-like. | LDLR binding is only one dimension of APOE biology. |
| APOE and AD/pathway context | Bu, G. (2009). Apolipoprotein E and its receptors in Alzheimer's disease: Pathways, pathogenesis and therapy. *Nature Reviews Neuroscience, 10*, 333-344. | https://doi.org/10.1038/nrn2620 | Supports the project-level relevance of APOE isoform modification in neurodegenerative disease context. | Broad review; not a parameter source for editing kinetics. |
| mRNA turnover prior | Wu, D., et al. (2020). Acrylonitrile-mediated nascent RNA sequencing for transcriptome-wide profiling of cellular RNA dynamics. *Advanced Science, 7*(8), 1902453. | https://pmc.ncbi.nlm.nih.gov/articles/PMC7175251/ | Supports broad mRNA degradation priors; HEK293T median mRNA half-life around 5.5 h corresponds to `k_deg` around 0.126 h^-1. | Median across many mRNAs; APOE/editor mRNA may differ. |
| Gene-expression kinetic model concept | Schwanhausser, B., Busse, D., Li, N., Dittmar, G., Schuchhardt, J., Wolf, J., Chen, W., & Selbach, M. (2011). Global quantification of mammalian gene expression control. *Nature, 473*, 337-342. | https://doi.org/10.1038/nature10098 | Supports modeling transcription, translation, mRNA degradation, and protein degradation as separate ODE layers. | Mouse fibroblast global dataset; not HEK293T PUF-A3A specific. |
| Protein turnover range | Boisvert, F.-M., Ahmad, Y., Gierlinski, M., Charriere, F., Lamont, D., Scott, M., Barton, G., & Lamond, A. I. (2012). A quantitative spatial proteomics analysis of proteome turnover in human cells. *Molecular & Cellular Proteomics, 11*(3), M111.011429. | https://doi.org/10.1074/mcp.M111.011429 | Supports broad `k_deg_p` prior because mammalian protein half-lives vary widely. | Not specific to PUF-APOBEC fusion proteins. |
| AAV PBPK extension | Zhang, Y., et al. (2024). Whole-body disposition and physiologically based pharmacokinetic modeling of adeno-associated viruses and the transgene product. *Journal of Pharmaceutical Sciences, 113*(1), 141-157. | https://doi.org/10.1016/j.xphs.2023.10.005 | Supports using liver as a major AAV sink and brain as lower but nonzero exposure in the translational extension. | Mouse AAV8/9 mAb-expression context; not our exact capsid or route. |
| AAV9 BBB concept | Duque, S., Joussemet, B., Riviere, C., Marais, T., Dubreil, L., Douar, A. M., Fyfe, J., Moullier, P., Colle, M. A., & Barkats, M. (2009). Intravenous administration of self-complementary AAV9 enables transgene delivery to adult motor neurons. *Molecular Therapy, 17*(7), 1187-1196. | https://doi.org/10.1038/mt.2009.71 | Supports the nonzero blood-to-brain/Bbb entry concept in our optional AAV-BBB module. | AAV9 CNS delivery is species-, age-, route-, and capsid-dependent. |
| AAV CNS vs liver tradeoff | Challis, R. C., et al. (2019). Systemic AAV vectors for widespread and targeted gene delivery in rodents. *Nature Protocols, 14*, 379-414. | https://doi.org/10.1038/s41596-018-0097-3 | Supports keeping liver/peripheral exposure as a safety penalty when discussing systemic AAV delivery. | Protocol/review-style source; not a parameter-calibration dataset. |
| ODE identifiability | Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmuller, U., & Timmer, J. (2009). Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. *Bioinformatics, 25*(15), 1923-1929. | https://doi.org/10.1093/bioinformatics/btp358 | Supports future profile-likelihood analysis and warns that `k_on`, `k_off`, and `k_cat` may not be separately identifiable. | Method paper; needs real time-course data before implementation is meaningful. |
| Model-informed experimental design | Banga, J. R., & Balsa-Canto, E. (2008). Parameter estimation and optimal experimental design. *Essays in Biochemistry, 45*, 195-209. | https://doi.org/10.1042/BSE0450195 | Supports using the model to choose next experiments rather than only fit past data. | General systems biology method; needs project-specific experiment costs and readouts. |
| Experimental design / model reduction | Jeong, J. E., Zhuang, Q., Transtrum, M. K., Zhou, E., & Qiu, P. (2018). Experimental design and model reduction in systems biology. *Quantitative Biology, 6*, 287-306. | https://doi.org/10.1007/s40484-018-0150-9 | Supports our DBTL cycle: identify uncertain parameters, recommend experiments, reduce model complexity if parameters cannot be identified. | Method reference; does not provide biological constants. |
| ODE validation and model selection | Hasdemir, D., Hoefsloot, H. C. J., & Smilde, A. K. (2015). Validation and selection of ODE based systems biology models: How to arrive at more reliable decisions. *BMC Systems Biology, 9*, 32. | https://doi.org/10.1186/s12918-015-0180-0 | Supports comparing alternative editing models and validating on held-out conditions. | Needs enough wet-lab conditions; not useful if all data are placeholders. |

---

## 4. What We Should Say on the Wiki

### Short text for the Model page

We used literature-derived priors to initialize our model before our own
wet-lab data were available. The most important source was the CU-REWIRE paper
by Han et al. (2022), which defined the PUF-APOBEC3A RNA-editing architecture,
the 48 h HEK293T assay setting, the preferred editing window, UC-context
preference, and the need to distinguish PUF-mediated off-targets from
deaminase-only background activity. PUF engineering papers were used to justify
sequence-aware binding scores and the 8R/10R design comparison. APOE structure
and receptor-binding literature was used to split the editing products into
APOE3-like and APOE2-like states. Gene-expression and AAV-distribution papers
were used only as broad priors for expression and translational delivery
extensions.

These parameters are not final biological constants. They are starting priors
that make the model mechanistically reasonable and testable. In later DBTL
cycles, qPCR, Western blot, APOE amplicon sequencing, local bystander analysis,
off-target amplicon panels, and RNA-seq will be used to fit or replace them.

### Recommended wording for parameter honesty

Because our wet-lab data are still incomplete, we report most kinetic constants
as normalized effective parameters. For example, `k_cat_112`, `k_on_112`, and
`k_off_112` should not be interpreted as separately measured biochemical
constants. With editing time-course data alone, these may only be identifiable
as an effective editing strength such as `k_cat / K_D`. Therefore, we keep
uncertainty ranges and use model predictions to rank designs, not to claim exact
absolute rates.

---

## 5. Parameters That Most Need Our Own Wet-Lab Data

| Priority | Parameter(s) | Required wet-lab data | Why it matters |
|---:|---|---|---|
| 1 | `k_cat_112`, `k_cat_158`, `k_on_112`, `k_on_158` | APOE112/APOE158 amplicon editing at 24/48/72 h under low/medium/high plasmid dose | Determines whether the model can predict real APOE editing. |
| 2 | `apoe3_like_fraction`, `apoe2_like_fraction` | Linked APOE112/APOE158 amplicon or long-read sequencing | Distinguishes desired APOE3-like product from double-edited APOE2-like product. |
| 3 | `local_bystander_per_112`, `local_bystander_per_158` | Amplicon sequencing of nearby Cs around APOE target windows | Required to know whether correct targeting still creates local collateral edits. |
| 4 | `k_on_puf_off`, `k_cat_puf_off` | Top-K predicted off-target amplicon panel and/or RNA-seq | Calibrates PUF-mediated off-target risk. |
| 5 | `k_deaminase_bg` | Mock, PUF-only, inactive A3A, and free-A3A controls | Separates targeting failure from intrinsic deaminase background. |
| 6 | `k_tx`, `k_tl`, `k_deg_m`, `k_deg_p` | Editor mRNA qPCR and protein Western/fluorescence time courses | Tells whether poor editing is caused by weak expression or weak catalysis. |
| 7 | AAV/BBB parameters | Blood/liver/brain vector-genome qPCR after delivery | Needed only if we move from HEK293T digital twin to in vivo delivery modeling. |

---

## 6. Full Reference List in APA Style

Banga, J. R., & Balsa-Canto, E. (2008). Parameter estimation and optimal
experimental design. *Essays in Biochemistry, 45*, 195-209.
https://doi.org/10.1042/BSE0450195

Boisvert, F.-M., Ahmad, Y., Gierlinski, M., Charriere, F., Lamont, D.,
Scott, M., Barton, G., & Lamond, A. I. (2012). A quantitative spatial
proteomics analysis of proteome turnover in human cells. *Molecular &
Cellular Proteomics, 11*(3), M111.011429.
https://doi.org/10.1074/mcp.M111.011429

Bu, G. (2009). Apolipoprotein E and its receptors in Alzheimer's disease:
Pathways, pathogenesis and therapy. *Nature Reviews Neuroscience, 10*, 333-344.
https://doi.org/10.1038/nrn2620

Campbell, Z. T., Valley, C. T., & Wickens, M. (2014). A protein-RNA
specificity code enables targeted activation of an endogenous human transcript.
*Nature Structural & Molecular Biology, 21*(8), 732-738.
https://doi.org/10.1038/nsmb.2847

Cheong, C.-G., & Hall, T. M. T. (2006). Engineering RNA sequence specificity of
Pumilio repeats. *Proceedings of the National Academy of Sciences, 103*(37),
13635-13639. https://doi.org/10.1073/pnas.0606294103

Duque, S., Joussemet, B., Riviere, C., Marais, T., Dubreil, L., Douar, A. M.,
Fyfe, J., Moullier, P., Colle, M. A., & Barkats, M. (2009). Intravenous
administration of self-complementary AAV9 enables transgene delivery to adult
motor neurons. *Molecular Therapy, 17*(7), 1187-1196.
https://doi.org/10.1038/mt.2009.71

Han, W., Huang, W., Wei, T., Ye, Y., Mao, M., & Wang, Z. (2022).
Programmable RNA base editing with a single gRNA-free enzyme. *Nucleic Acids
Research, 50*(16), 9580-9595. https://doi.org/10.1093/nar/gkac713

Han, W., Yuan, B., Fan, X., Li, W., Yuan, Y., Zhang, Y., Wang, S., Shan, S.,
Hafner, M., Wang, Z., & Qiu, Z. (2025). Effective in vivo RNA base editing via
engineered cytidine deaminase APOBECs fused with PUF proteins. *Nature
Communications, 16*, 9727. https://doi.org/10.1038/s41467-025-64748-6

Hasdemir, D., Hoefsloot, H. C. J., & Smilde, A. K. (2015). Validation and
selection of ODE based systems biology models: How to arrive at more reliable
decisions. *BMC Systems Biology, 9*, 32.
https://doi.org/10.1186/s12918-015-0180-0

Hatters, D. M., Peters-Libeu, C. A., & Weisgraber, K. H. (2006).
Apolipoprotein E structure: Insights into function. *Trends in Biochemical
Sciences, 31*(8), 445-454. https://doi.org/10.1016/j.tibs.2006.06.008

Huang, X., Lv, J., Li, Y., Mao, S., Li, Z., Jing, Z., Sun, Y., Zhang, X.,
Shen, S., Wang, X., Di, M., Ge, J., Huang, X., Zuo, E., & Chi, T. (2020).
Programmable C-to-U RNA editing using the human APOBEC3A deaminase. *The EMBO
Journal, 39*(22), e104741. https://doi.org/10.15252/embj.2020104741

Huang, Y., & Mahley, R. W. (2014). Apolipoprotein E: Structure and function in
lipid metabolism, neurobiology, and Alzheimer's diseases. *Neurobiology of
Disease, 72*, 3-12. https://doi.org/10.1016/j.nbd.2014.08.025

Jeong, J. E., Zhuang, Q., Transtrum, M. K., Zhou, E., & Qiu, P. (2018).
Experimental design and model reduction in systems biology. *Quantitative
Biology, 6*, 287-306. https://doi.org/10.1007/s40484-018-0150-9

Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmuller,
U., & Timmer, J. (2009). Structural and practical identifiability analysis of
partially observed dynamical models by exploiting the profile likelihood.
*Bioinformatics, 25*(15), 1923-1929.
https://doi.org/10.1093/bioinformatics/btp358

Schwanhausser, B., Busse, D., Li, N., Dittmar, G., Schuchhardt, J., Wolf, J.,
Chen, W., & Selbach, M. (2011). Global quantification of mammalian gene
expression control. *Nature, 473*, 337-342.
https://doi.org/10.1038/nature10098

Wu, D., et al. (2020). Acrylonitrile-mediated nascent RNA sequencing for
transcriptome-wide profiling of cellular RNA dynamics. *Advanced Science, 7*(8),
1902453. https://pmc.ncbi.nlm.nih.gov/articles/PMC7175251/

Zhao, Y.-Y., Mao, M.-W., Zhang, W.-J., Wang, J., Li, H.-T., Yang, Y.,
Wang, Z., & Wu, J.-W. (2018). Expanding RNA binding specificity and affinity
of engineered PUF domains. *Nucleic Acids Research, 46*(9), 4771-4782.
https://doi.org/10.1093/nar/gky134

Zhang, Y., et al. (2024). Whole-body disposition and physiologically based
pharmacokinetic modeling of adeno-associated viruses and the transgene product.
*Journal of Pharmaceutical Sciences, 113*(1), 141-157.
https://doi.org/10.1016/j.xphs.2023.10.005
