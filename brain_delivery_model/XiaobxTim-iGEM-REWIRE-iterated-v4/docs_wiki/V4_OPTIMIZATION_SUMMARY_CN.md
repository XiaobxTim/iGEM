# v4 优化方案总结：从“能跑的模拟”到“能指导实验的数字孪生”

## 一句话版

v3 已经能把 AAV 递送、脑内表达和 RNA 编辑串成一条完整链路；v4 的目标是让模型更像一个实验设计助手：它不只是画曲线，而是帮我们判断哪个 PUF-deaminase 构型、哪个表达水平、哪个采样时间更值得先做实验。

## 优化方案到底是什么意思

可以把 v3 理解成一条“快递路线”：

药物从给药部位出发，进入血液，过 BBB，到脑细胞里表达编辑器，最后产生一个 on-target 和一个 off-target 结果。

优化方案指出：这条路线已经有了，但最后的“编辑结果”说得太粗。原来模型像是在说“包裹送到了，产生了好结果或坏结果”；v4 要进一步问：

- 送到的是 APOE 的哪个位点？
- 只改 112 位，还是 112 和 158 位都改了？
- 这个产物应该叫 APOE3-like，还是 APOE2-like？
- 脱靶到底是目标附近误改、PUF 找错了相似序列，还是 deaminase 自己造成背景编辑？
- 一条候选 PUF 序列为什么会更快结合、更容易催化，还是更容易脱靶？
- 如果只能先做几组实验，应该优先测哪几个设计？

所以 v4 的核心不是“多加几个方程显得复杂”，而是把模型变成一个可解释的设计选择器。

## v4 已经完成的主要升级

### 1. APOE 112/158 双位点建模

v3 中只有一个笼统的 `S_on -> E_on`，容易把产物直接写成 APOE2-like。

v4 改成了 APOE4 的两个关键位点：

- APOE4: 112=Arg, 158=Arg
- APOE3-like: 112 被编辑，158 仍为 Arg
- APOE2-like: 112 和 158 都被编辑
- APOE158-only / mixed-edit: 只发生 158 位编辑或不完整编辑

这让模型能回答：我们当前设计更像是在做 APOE4 -> APOE3，还是在推动双位点 APOE2-like。

### 2. 区分 editor 类型

v4 在参数里加入 `editor_type`，可以区分：

- `A3A`
- `APOBEC1`
- `engineered_A3A`
- `ProAPOBEC`

不同 editor 会影响催化效率、序列偏好和背景脱靶风险，避免把不同 deaminase 的文献参数混在一起解释。

### 3. 把 off-target 拆成三类

v3 把脱靶压缩成一个 `S_off`。

v4 拆成：

- `B_local_bystander`: PUF 结合对了，但附近其他 C 被顺手编辑；
- `S_puf_off/C_puf_off/E_puf_off`: PUF 结合到相似 RNA 位点造成的脱靶；
- `S_deaminase_bg/E_deaminase_bg`: deaminase 背景活性造成的非靶向编辑。

这三类对应不同实验验证方式，也对应不同工程优化策略。

### 4. 加入 sequence-to-kinetics

v4 新增候选 PUF 设计表。模型会把设计特征转换成动力学参数：

- PUF repeat 数；
- PUF score；
- mismatch 数；
- RNA accessibility；
- 编辑 C 与 PUF 结合区距离；
- UC context；
- editor 类型。
- editor 活性比例，用来表示 inactive deaminase 或 PUF-only 对照。

这些特征会影响 `k_on_112`、`k_cat_112`、`k_on_158`、`k_cat_158`、local bystander 和 background editing。

也就是说，v4 不再只是手填参数，而是让“候选序列设计”真正进入 ODE 模型。

### 5. 加入设计筛选和 Pareto 推荐

v4 的分析脚本会扫描：

- 多个 PUF 候选构型；
- 低 / 中 / 高表达水平；
- 24 / 48 / 72 / 96 h 采样时间。

然后综合：

- APOE3-like 目标编辑；
- APOE2-like / LDL receptor binding 风险代理；
- local bystander；
- PUF mismatch off-target；
- deaminase background；
- brain expression burden；
- overall specificity。

最终输出设计推荐表和 Pareto 前沿，不只给一条曲线。

## v4 主要文件地址

### 核心状态和仿真

- `models/full_model/state_vector.py`
- `models/full_model/simulator.py`

### RNA editing 核心模型

- `models/editing/apoe_multisite_editing.py`
- `models/editing/module5.py`
- `models/editing/sequence_to_kinetics.py`
- `models/editing/offtarget_panel.py`

### 参数和候选设计

- `config/base_config.yaml`
- `config/parameter_provenance.yaml`
- `wetlab/templates/puf_design_candidates_template.csv`
- `wetlab/templates/assay_observations_template.csv`

### 湿实验连接和拟合

- `models/calibration/wetlab_bridge.py`
- `models/calibration/parameter_fit.py`

### 自动分析和输出

- `scripts/run_iterated_analysis.py`
- `reports/iterated_analysis/v4_design_screen.csv`
- `reports/iterated_analysis/v4_design_recommendation.md`

### 测试

- `tests/test_apoe_multisite_model.py`
- `tests/test_sequence_to_kinetics.py`

## 仍需真实实验或公开数据继续补强的地方

v4 已经搭好结构，但下面这些属于下一轮可以继续强化的内容：

- 用 REWIRE 公开数据训练 sequence-to-kinetics 系数；
- 用 ProAPOBEC 等外部数据做冻结测试；
- 做 parameter recovery；
- 做 profile likelihood / identifiability；
- 做多机制模型比较；
- 把你们自己的 qPCR、Western blot、amplicon-seq 数据填入模板后重新拟合。

这部分不是没做，而是 v4 已经把入口留好了，等数据或公开表格整理好后可以直接接入。
