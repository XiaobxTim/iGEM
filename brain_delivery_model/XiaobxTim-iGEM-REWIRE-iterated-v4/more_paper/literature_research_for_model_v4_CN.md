# v4 模型文献调研整理

> 目的：这份文档把“论文 + 有价值的部分 + 为什么可信 / 怎么用到当前模型 + 缺点隐患”整理成方便阅读和后续 wiki 改写的形式。它不是单纯参考文献列表，而是解释每篇文献如何支持我们当前 v4 模型。

## 1. 总览：这些文献分别支持模型的哪一部分

| 模型问题 | 最相关文献 | 用到 v4 哪里 |
|---|---|---|
| PUF-APOBEC RNA editing 的核心机制 | Han et al., 2022, NAR, CU-REWIRE | `editor_type=A3A`、PUF 后 +2 编辑窗口、UC context、10R PUF、off-target 拆分 |
| engineered deaminase / ProAPOBEC 的扩展 | Han et al., 2025, Nature Communications | `editor_type=ProAPOBEC`、`proapobec_cat_scale`、`proapobec_background_scale` |
| PUF repeat number 和 RNA binding specificity | Zhao et al., 2018; Cheong & Hall, 2006; Campbell et al., 2014 | `puf_repeats`、`puf_score`、`mismatch_count`、sequence-to-kinetics |
| APOBEC3A motif / accessibility | Huang et al., 2020, EMBO Journal | `uc_context_112/158`、`accessibility_112/158`、`k_deaminase_bg` |
| APOE2/3/4 生物学定义和风险 | Hatters et al., 2006; Huang & Mahley, 2014; Bu, 2009 | APOE112/158 双位点状态、`ldlr_binding_risk_proxy` |
| mRNA / protein 表达动力学 | Wu et al., 2020; Schwanhausser et al., 2011; Boisvert et al., 2012 | `k_tx`、`k_tl`、`k_deg_m`、`k_deg_p` |
| AAV/BBB 递送扩展 | Duque et al., 2009; Zhang et al., 2024; Challis et al., 2019 | AAV-BBB translational extension、liver/brain exposure |
| 模型可信度和实验设计方法 | Raue et al., 2009; Banga & Balsa-Canto, 2008; Jeong et al., 2018; Hasdemir et al., 2015 | identifiability、parameter fitting、实验推荐、模型比较 |

---

## 2. 核心 RNA Editing 文献

### Paper 1: CU-REWIRE 原始框架

**APA**

Han, W., Huang, W., Wei, T., Ye, Y., Mao, M., & Wang, Z. (2022). Programmable RNA base editing with a single gRNA-free enzyme. *Nucleic Acids Research, 50*(16), 9580-9595. https://doi.org/10.1093/nar/gkac713

**有价值的部分**

- 证明 PUF + APOBEC3A 可以构成 gRNA-free C-to-U RNA editor。
- CU-REWIRE 使用 APOBEC3A/A3A，而不是泛泛的 APOBEC。
- HEK293T 转染实验常用 48 h 作为检测时间。
- C-to-U 编辑有窄编辑窗口，主要偏向 PUF binding site 后第 2 位。
- UC/UCN context 明显提高编辑效率，差 context 会降低编辑。
- 10-repeat PUF 相比短识别长度有更好 transcriptome specificity。
- free APOBEC3A 可产生 deaminase-only RNA off-target。
- inactive A3A E72A 可作为失活编辑对照。
- 公开数据可用于未来 benchmark：GEO GSE155734。

**为什么可信 / 怎么用到当前模型**

这是 REWIRE / CU-REWIRE 的原始核心论文，和我们项目的 PUF-A3A RNA editing 架构最直接对应。v4 中以下内容直接或间接受它指导：

- `editing.editor_type: A3A`
- `simulation.default_t_end: 48.0`
- `distance_112/158 -> kcat_scale`
- `uc_context_112/158 -> kcat_scale`
- `puf_repeats -> puf_offtarget_scale`
- `editor_activity_scale` 用于 inactive editor / PUF-only control
- off-target 拆为 local bystander、PUF mismatch、deaminase background

**缺点 / 隐患**

- 论文中很多结果来自 EGFP、CTNNB1、EZH2、SCN1A 等位点，不是 APOE。
- 论文中的编辑率不能直接当作 APOE112 的真实编辑率。
- 48 h 是 HEK293T 转染体系参考点，不等于 AAV-BBB 全链路的最优时间。
- 我们当前把论文规则转成指数函数和比例系数，这些系数仍需用公开数据和自有实验重新拟合。

---

### Paper 2: ProAPOBEC / 下一代 PUF-APOBEC

**APA**

Han, W., Yuan, B., Fan, X., Li, W., Yuan, Y., Zhang, Y., Wang, S., Shan, S., Hafner, M., Wang, Z., & Qiu, Z. (2025). Effective in vivo RNA base editing via engineered cytidine deaminase APOBECs fused with PUF proteins. *Nature Communications, 16*, 9727. https://doi.org/10.1038/s41467-025-64748-6

**有价值的部分**

- 展示 engineered APOBEC / ProAPOBEC 与 PUF 融合后可以提升 RNA editing 表现。
- 支持把 deaminase identity 作为模型变量，而不是只写“APOBEC”。
- 对 APOE 相关 C-to-U editing 有直接参考价值。
- 提示 engineered deaminase 可能同时影响 on-target activity 和 background off-target。

**为什么可信 / 怎么用到当前模型**

这是 REWIRE 体系后续升级方向，适合用来扩展 v4 的候选构型。v4 中：

- `editor_type=ProAPOBEC`
- `proapobec_cat_scale: 1.25`
- `proapobec_background_scale: 0.5`
- 候选表中有 `10R-proapobec`

这些不是最终实测值，而是允许模型比较 A3A 与 engineered deaminase 的文献先验。

**缺点 / 隐患**

- 2025 新工作较新，外部重复验证还少。
- ProAPOBEC 的性能高度 variant-specific，不能用一个固定倍数代表所有 engineered APOBEC。
- 如果我们没有选择具体 ProAPOBEC 序列，这部分只能作为未来方向。

---

### Paper 3: APOBEC3A C-to-U RNA editing motif

**APA**

Huang, X., Lv, J., Li, Y., Mao, S., Li, Z., Jing, Z., Sun, Y., Zhang, X., Shen, S., Wang, X., Di, M., Ge, J., Huang, X., Zuo, E., & Chi, T. (2020). Programmable C-to-U RNA editing using the human APOBEC3A deaminase. *The EMBO Journal, 39*(22), e104741. https://doi.org/10.15252/embj.2020104741

**有价值的部分**

- 支持 A3A 对 RNA C-to-U editing 的 UC context 偏好。
- 支持 RNA accessibility / local RNA structure 会影响编辑。
- 支持 APOBEC3A 可能产生 transcriptome-wide off-target，需要单独建模。

**为什么可信 / 怎么用到当前模型**

这篇不是 PUF-REWIRE，但它是 APOBEC3A RNA editing 的关键机制参考。v4 用它来支持：

- `uc_context_112`
- `uc_context_158`
- `accessibility_112`
- `accessibility_158`
- `k_deaminase_bg`

**缺点 / 隐患**

- 这篇的系统和 REWIRE 不完全一样，空间窗口不能直接照搬。
- 它支持 motif 和 background 方向，但不能给出我们 APOE 位点的最终参数。

---

## 3. PUF 设计与 RNA Binding 文献

### Paper 4: PUF repeat number 和 binding specificity

**APA**

Zhao, Y.-Y., Mao, M.-W., Zhang, W.-J., Wang, J., Li, H.-T., Yang, Y., Wang, Z., & Wu, J.-W. (2018). Expanding RNA binding specificity and affinity of engineered PUF domains. *Nucleic Acids Research, 46*(9), 4771-4782. https://doi.org/10.1093/nar/gky134

**有价值的部分**

- 比较 8R、9R、10R、12R、16R 等 engineered PUF。
- 9R / 10R 往往可以提升 specificity / affinity。
- 不是 repeat 越多越好，过长可能带来结构和功能风险。

**为什么可信 / 怎么用到当前模型**

这篇直接支持我们把 `puf_repeats` 作为设计变量。v4 中：

- 8R design 保留为 legacy / lower specificity 对照；
- 10R design 作为默认推荐方向；
- `puf_offtarget_scale` 会因 10R 而降低。

**缺点 / 隐患**

- PUF binding affinity 不等于最终 RNA editing efficiency。
- 需要结合 APOBEC catalytic window 和 RNA accessibility 一起判断。

---

### Paper 5: PUF recognition code 工程化

**APA**

Cheong, C.-G., & Hall, T. M. T. (2006). Engineering RNA sequence specificity of Pumilio repeats. *Proceedings of the National Academy of Sciences, 103*(37), 13635-13639. https://doi.org/10.1073/pnas.0606294103

**有价值的部分**

- 证明 PUF repeat 可以通过氨基酸识别代码设计 RNA base specificity。
- 支持“PUF target sequence 可以理性设计”这个建模前提。

**为什么可信 / 怎么用到当前模型**

v4 当前先用 `puf_score` 和 `mismatch_count` 代表 PUF-target match。下一步可以根据这类 recognition code 直接从 PUF repeat 序列计算 `puf_score`。

**缺点 / 隐患**

- 不是 RNA editing 论文。
- 不提供 APOBEC catalytic 参数。
- 需要我们知道真实 PUF repeat 氨基酸序列才能精确应用。

---

### Paper 6: PUF code 应用于内源 RNA

**APA**

Campbell, Z. T., Valley, C. T., & Wickens, M. (2014). A protein-RNA specificity code enables targeted activation of an endogenous human transcript. *Nature Structural & Molecular Biology, 21*(8), 732-738. https://doi.org/10.1038/nsmb.2847

**有价值的部分**

- 说明 engineered PUF 可以被用于人细胞中内源 RNA 的特异识别和调控。
- 支持 PUF 作为 programmable RNA-binding domain 的可设计性。

**为什么可信 / 怎么用到当前模型**

支撑 v4 的 sequence-to-kinetics 思路：候选 PUF 序列不是随便给参数，而是可以作为工程设计输入。

**缺点 / 隐患**

- 不是 base editing。
- 不能直接推断 editing rate 或 off-target editing。

---

## 4. APOE 生物学文献

### Paper 7: APOE isoform 结构定义

**APA**

Hatters, D. M., Peters-Libeu, C. A., & Weisgraber, K. H. (2006). Apolipoprotein E structure: Insights into function. *Trends in Biochemical Sciences, 31*(8), 445-454. https://doi.org/10.1016/j.tibs.2006.06.008

**有价值的部分**

- APOE4: Arg112 / Arg158
- APOE3: Cys112 / Arg158
- APOE2: Cys112 / Cys158
- 说明 112 和 158 两个位点共同决定 isoform。

**为什么可信 / 怎么用到当前模型**

v4 直接用这个生物学定义修正 v3 的问题：

- 只编辑 112 位不再叫 APOE2-like，而是 APOE3-like；
- 同时编辑 112 和 158 才叫 APOE2-like；
- 只编辑 158 作为 mixed / incomplete edit 追踪。

**缺点 / 隐患**

- 这是蛋白 isoform 结构文献，不是 RNA editing 速率文献。
- RNA 层面的编辑比例不一定等于蛋白层面的功能恢复比例。

---

### Paper 8: APOE receptor-binding 与风险解释

**APA**

Huang, Y., & Mahley, R. W. (2014). Apolipoprotein E: Structure and function in lipid metabolism, neurobiology, and Alzheimer's diseases. *Neurobiology of Disease, 72*, 3-12. https://doi.org/10.1016/j.nbd.2014.08.025

**有价值的部分**

- APOE isoform 不只是标签差异，会影响 lipid transport、receptor binding 和 Alzheimer disease 相关路径。
- APOE2-like 双位点状态可能带来 LDL receptor binding 解释风险。

**为什么可信 / 怎么用到当前模型**

v4 增加 `ldlr_binding_risk_proxy`，把 APOE2-like fraction 作为需要监控的风险，而不是简单认为“编辑越多越好”。

**缺点 / 隐患**

- APOE2-like 是否一定不适合治疗，需要结合项目目标和疾病背景判断。
- 我们目前只是风险 proxy，没有建完整脂质代谢网络。

---

### Paper 9: APOE 与 Alzheimer disease 治疗背景

**APA**

Bu, G. (2009). Apolipoprotein E and its receptors in Alzheimer's disease: Pathways, pathogenesis and therapy. *Nature Reviews Neuroscience, 10*, 333-344. https://doi.org/10.1038/nrn2620

**有价值的部分**

- 支持 APOE isoform 是 Alzheimer disease 相关的重要治疗方向。
- 强调 APOE/receptor/lipid metabolism 是系统性问题。

**为什么可信 / 怎么用到当前模型**

用于支撑项目故事：我们为什么关心 APOE4-to-APOE3-like RNA editing。也解释为什么模型需要区分 APOE3-like 和 APOE2-like，而不是只输出 on-target editing。

**缺点 / 隐患**

- 综述性质，不提供 kinetic 参数。
- 不能单独作为模型参数来源。

---

## 5. 表达、递送和建模方法文献

### Paper 10: mRNA half-life

**APA**

Wu, D., et al. (2020). Acrylonitrile-mediated nascent RNA sequencing for transcriptome-wide profiling of cellular RNA dynamics. *Advanced Science, 7*(8), 1902453. https://pmc.ncbi.nlm.nih.gov/articles/PMC7175251/

**有价值的部分**

- HEK293T mRNA half-life 在数小时量级。
- median half-life 可用于估计 `k_deg = ln(2) / t_half`。

**为什么可信 / 怎么用到当前模型**

v4 中 `intracellular.k_deg_m` 和 `editing.k_deg_apoe` 使用 broad prior，而不是随便固定。5.5 h half-life 对应约 `0.126 h^-1`，和 v4 的 `0.08-0.1 h^-1` 在同一数量级。

**缺点 / 隐患**

- median mRNA half-life 不等于 APOE mRNA 或 PUF-A3A mRNA 的 half-life。
- 仍需 qPCR time course 校准。

---

### Paper 11: Gene expression ODE 结构

**APA**

Schwanhausser, B., Busse, D., Li, N., Dittmar, G., Schuchhardt, J., Wolf, J., Chen, W., & Selbach, M. (2011). Global quantification of mammalian gene expression control. *Nature, 473*, 337-342. https://doi.org/10.1038/nature10098

**有价值的部分**

- 支持把 mRNA abundance、translation、protein degradation 分开建模。

**为什么可信 / 怎么用到当前模型**

v4 module 4 用：

- `k_tx`
- `k_deg_m`
- `k_tl`
- `k_deg_p`

这对应 D -> M -> P 的表达动力学，而不是直接假设编辑器蛋白瞬间出现。

**缺点 / 隐患**

- 不是 PUF-A3A fusion protein。
- 参数仍需我们自己的 qPCR/Western fitting。

---

### Paper 12: AAV9 BBB / CNS delivery

**APA**

Duque, S., Joussemet, B., Riviere, C., Marais, T., Dubreil, L., Douar, A. M., Fyfe, J., Moullier, P., Colle, M. A., & Barkats, M. (2009). Intravenous administration of self-complementary AAV9 enables transgene delivery to adult motor neurons. *Molecular Therapy, 17*(7), 1187-1196. https://doi.org/10.1038/mt.2009.71

**有价值的部分**

- 支持 AAV9 类载体可以有非零 CNS delivery。

**为什么可信 / 怎么用到当前模型**

v4 保留 AAV-BBB 作为 translational extension。`distribution.k_blood_to_brain` 和 BBB transport 参数因此不设为 0。

**缺点 / 隐患**

- AAV9 delivery 强烈依赖动物、年龄、capsid、给药方式。
- 当前 HEK293T 核心模型不应该过度强调 AAV-BBB。

---

### Paper 13: AAV PBPK / liver-brain 分布

**APA**

Zhang, Y., et al. (2024). Whole-body disposition and physiologically based pharmacokinetic modeling of adeno-associated viruses and the transgene product. *Journal of Pharmaceutical Sciences, 113*(1), 141-157. https://doi.org/10.1016/j.xphs.2023.10.005

**有价值的部分**

- AAV liver exposure 通常远高于 brain exposure。
- PBPK 框架可支持多组织分布建模。

**为什么可信 / 怎么用到当前模型**

v4 的 AAV 扩展中保留 liver AUC 和 brain exposure，让模型不仅看 editing，还看递送安全性。

**缺点 / 隐患**

- PBPK 数据不是我们具体 capsid / dose / route。
- 不能直接用于 HEK293T 实验。

---

### Paper 14: Identifiability / profile likelihood

**APA**

Raue, A., Kreutz, C., Maiwald, T., Bachmann, J., Schilling, M., Klingmuller, U., & Timmer, J. (2009). Structural and practical identifiability analysis of partially observed dynamical models by exploiting the profile likelihood. *Bioinformatics, 25*(15), 1923-1929. https://doi.org/10.1093/bioinformatics/btp358

**有价值的部分**

- ODE 模型中参数可能不可识别。
- profile likelihood 可以判断参数是否能被数据唯一约束。

**为什么可信 / 怎么用到当前模型**

提醒我们不要声称 `k_on`、`k_off`、`k_cat` 都是精确值。未来可以加入 profile likelihood，判断哪些参数只能作为 `k_cat/KD` 这样的组合被识别。

**缺点 / 隐患**

- 方法论文，不提供生物参数。
- 需要真实 time-course 数据后才有意义。

---

### Paper 15: Optimal experimental design

**APA**

Banga, J. R., & Balsa-Canto, E. (2008). Parameter estimation and optimal experimental design. *Essays in Biochemistry, 45*, 195-209. https://doi.org/10.1042/BSE0450195

**有价值的部分**

- 模型不应该只拟合已有数据，也应该指导下一轮实验。

**为什么可信 / 怎么用到当前模型**

v4 的 design screen 和 recommendation 就是这个思想：让模型推荐设计、表达水平、时间点，而不是等实验做完才解释。

**缺点 / 隐患**

- 我们当前 utility weights 仍是团队目标假设。
- 后续需要根据实验成本、时间和风险重新设定权重。

---

## 6. 可以直接告诉评委的总结

我们不是简单把论文里的参数照搬进模型，而是把可靠文献转化为“可被实验更新的先验”。CU-REWIRE 论文告诉我们 PUF-A3A 的编辑窗口、UC context、10R PUF specificity 和 deaminase background；APOE 文献告诉我们必须区分 APOE3-like 和 APOE2-like；PUF engineering 文献告诉我们 PUF sequence 可以进入 kinetic 参数；系统生物学文献告诉我们需要 uncertainty、identifiability 和 model-guided experiment。

因此，v4 模型的意义是：

1. 用文献规则建立合理初始参数；
2. 用 ODE 模拟从表达、结合、催化到脱靶的过程；
3. 用 design screen 选择下一轮最值得做的实验；
4. 用湿实验数据不断替换文献先验，形成 DBTL 闭环。
