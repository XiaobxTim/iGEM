# iGEM AAV-PUF-APOBEC Combined Model 说明

## 1. GitHub 仓库前半段在做什么

这个仓库原本是一个从给药到脑内编辑结果的多尺度 AAV 模型。它把整个过程拆成 6 个模块：

1. **Module 1: 吸收**
   AAV 从局部给药 depot 出发，一部分进入淋巴，一部分直接进入血液，还有一部分局部损失。

2. **Module 2: 全身分布**
   血液中的 AAV 在血液、肝脏、外周组织和脑血管侧之间转运，并考虑血液、肝脏、外周组织、脑血管侧清除。

3. **Module 3: BBB 转运**
   脑血管侧 AAV 先结合/进入脑内皮细胞表面，再进入内吞 compartment，最后有一部分成功跨 BBB 到脑间质液 ISF。

4. **Module 4: 脑细胞摄取与表达**
   AAV 从脑 ISF 进入脑细胞，再进入细胞核，经过转录得到 `mRNA_brain`，再翻译得到 `P_brain`。在整合模型里，`P_brain` 被解释为脑细胞内可用的 PUF-APOBEC 活性编辑器。

5. **Module 5: RNA 编辑**
   原仓库里这里是一个简化的竞争编辑模型：用 `P_brain` 直接驱动 on-target 和 off-target 编辑速率。

6. **Module 6: 剂量优化**
   在不同剂量下重复运行模型，提取 on-target 编辑率、off-target burden、specificity index、肝脏 AUC、血液 Cmax 等指标，判断治疗窗口。

## 2. 本地 notebook 后半段在做什么

`puf-apobec.ipynb` 主要建立了 PUF-APOBEC 对 APOE4 RNA 的酶促编辑模型：

```text
E + S <-> ES
ES -> E + P
```

其中：

- `E` 是游离活性 PUF-APOBEC 编辑器。
- `S` 是未编辑、可接近的 APOE4 mRNA。
- `ES` 是编辑器和底物形成的复合物。
- `P` 是编辑后的 APOE2-like RNA。

notebook 还定义了一个人为输入函数：

```text
u(t) = u_max * (1 - exp(-k_exp * t)) * exp(-k_decay * t)
```

它表示 AAV 表达导致的编辑器生成速率。然后通过 ODE 计算编辑效率：

```text
editing_fraction = P / (S + P)
```

notebook 最后还做了参数扫描，例如扫描 `u_max` 和 `k_cat`，画出最终编辑效率 heatmap，用于判断表达强度和催化效率对编辑结果的影响。

## 3. 两部分如何 combine

整合时最关键的接口是：

```text
GitHub 前半段的 P_brain(t)  ->  notebook 后半段的 E(t)
```

也就是说，原 notebook 里人为假设的 `u(t)` 不再单独使用，因为 GitHub 模型已经显式模拟了：

```text
AAV 给药 -> 吸收 -> 全身分布 -> BBB -> 脑细胞摄取 -> 转录 -> 翻译 -> P_brain
```

所以整合后的模型直接把 `P_brain` 当作游离活性编辑器，让它参与 RNA 结合、解离和催化。

## 4. 新的 Module 5 方程

整合模型保留了 on-target 和 off-target 两条反应链：

```text
P_brain + S_on  <->  ES_on   ->  P_brain + E_on
P_brain + S_off <->  ES_off  ->  P_brain + E_off
```

状态变量含义：

- `P_brain`: 游离活性 PUF-APOBEC 编辑器，由 Module 4 产生。
- `S_on`: 未编辑 APOE4 target RNA。
- `S_off`: 可接近的 off-target RNA。
- `ES_on`: on-target 编辑器-底物复合物。
- `ES_off`: off-target 编辑器-底物复合物。
- `E_on`: on-target 编辑产物，也就是希望得到的 APOE2-like RNA。
- `E_off`: off-target 编辑产物，是安全性风险指标。

核心 ODE：

```text
dP_brain/dt += - bind_on - bind_off + unbind_on + unbind_off + edit_on + edit_off

dS_on/dt    = k_prod_on  - k_deg_on*S_on  - bind_on  + unbind_on
dS_off/dt   = k_prod_off - k_deg_off*S_off - bind_off + unbind_off

dES_on/dt   = bind_on  - unbind_on  - edit_on
dES_off/dt  = bind_off - unbind_off - edit_off

dE_on/dt    = edit_on  - k_loss_on*E_on
dE_off/dt   = edit_off - k_loss_off*E_off
```

其中：

```text
bind_on   = k_on_on  * P_brain * S_on
unbind_on = k_off_on * ES_on
edit_on   = k_cat_on * ES_on

bind_off   = k_on_off  * P_brain * S_off
unbind_off = k_off_off * ES_off
edit_off   = k_cat_off * ES_off
```

## 5. 改动位置

- `models/editing/puf_apobec_mechanistic.py`
  新增显式 PUF-APOBEC 结合、解离、催化模型，注释较详细。

- `models/editing/module5.py`
  新增编辑模块分发层。配置为 `mechanistic_puf_apobec` 时使用新模型；配置为 `competitive` 时仍可使用原始简化模型。

- `models/full_model/state_vector.py`
  在全局状态向量里加入 `ES_on` 和 `ES_off`。

- `models/full_model/rhs_aav.py`
  把 Module 5 调用改成按配置选择模型。

- `config/base_config.yaml`
  默认启用 `model: mechanistic_puf_apobec`，并加入 `k_on_on`、`k_off_on`、`k_on_off`、`k_off_off`、RNA 生成和降解参数。

- `utils/plotting.py`
  绘图中加入 `P_brain`、`ES_on`、`ES_off` 和 notebook 风格的 editing fraction。

- `optimization/metrics.py`
  剂量优化使用统一的 Module 5 指标函数。

## 6. 如何运行

在仓库目录运行：

```bash
python main.py --route footpad --dose 1.0 --t_end 72 --dt 0.2 --output_dir outputs/run_combined_puf_apobec
```

如果你使用的是本机 Anaconda 环境，可以运行：

```bash
/opt/anaconda3/bin/python main.py --route footpad --dose 1.0 --t_end 72 --dt 0.2 --output_dir outputs/run_combined_puf_apobec
```

运行后会生成：

- `module1_states.png`
- `module2_states.png`
- `module3_states.png`
- `module4_states.png`
- `module5_states.png`
- `module5_fluxes.png`
- `module5_metrics.png`
- `combined_states.png`
- `simulation_results.npz`

## 7. 当前验证结果示例

使用 `dose=1.0`、`t_end=72 h`、`dt=0.2 h` 的一次测试中：

- `P_brain` final = 0.00320
- `E_on` final = 0.01343
- `E_off` final = 0.00129
- final on-target editing fraction = 0.01337
- final off-target editing fraction = 0.000258
- final specificity index = 10.41

这些数值主要用于检查模型链条是否连通，不应直接解释为真实实验预测。后续应根据 wet-lab 数据校准 AAV 转运、表达、结合和催化参数。
