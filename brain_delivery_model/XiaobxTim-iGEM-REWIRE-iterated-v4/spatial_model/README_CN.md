# PhysiCell–ParaView 三维脑递送模型

这个扩展把原有的数值 ODE 模型变成一个可运行的三维反应–扩散/多细胞模型：

- Python Module 1–2 产生脑血侧 AAV 边界曲线；
- PhysiCell 1.14.2 计算 BBB 转运和 1480 个固定细胞；
- BioFVM 计算组织内 AAV 扩散、衰减和细胞摄取；
- 每个神经元和星形胶质细胞都运行完整的 Module 4–5；
- ParaView 6.0.1 数据包含 VTI、VTP、VTM 和 PVD 时间序列。

模型使用 400 × 300 × 300 μm 组织域、10 μm 网格、280 个内皮细胞、800 个神经元和 400 个星形胶质细胞。三组剂量为 0.3、1.0 和 3.0，细胞位置均使用固定随机种子 42。

## 一次运行

安装固定版本 PhysiCell（已经存在时不会重复下载）：

```bash
conda run -n xbx_env python -m spatial_model.setup_physicell
```

运行默认 72 小时模拟：

```bash
conda run -n xbx_env python -m spatial_model.run_pipeline
```

短测试可使用：

```bash
conda run -n xbx_env python -m spatial_model.run_pipeline \
  --output-dir outputs/physicell_brain_delivery_short \
  --duration-hours 1 --boundary-dt-min 1 --output-interval-min 10
```

输出目录采用防覆盖和原子发布：目标目录已存在时程序会停止，不会覆盖旧结果。

## 在 ParaView 中打开

直接打开：

```text
outputs/physicell_brain_delivery/comparison/comparison.pvd
```

这是三剂量并排动画，三个组织中心位于 x = −450、0、450 μm。单独剂量的入口是：

```text
dose_0p3/paraview/simulation.pvd
dose_1p0/paraview/simulation.pvd
dose_3p0/paraview/simulation.pvd
```

推荐显示方式：

1. 打开 PVD 后点击 `Apply`，用 `Extract Block` 选择 `microenvironment` 或 `cells`。
2. 对微环境使用 `Volume` 或 `Slice`，按 `log10_extracellular_AAV` 着色。
3. 对细胞使用 `Glyph`，Glyph 类型选 `Sphere`，按 `radius_um` 缩放。
4. 细胞按 `cell_type` 着色：0 为内皮，1 为神经元，2 为星形胶质细胞。
5. 可切换 `editing_fraction`、`off_target_burden`、`editor_protein`、`distance_to_vessel_um` 等数组。

## 在 PhysiCell Studio 中查看/修改

每个剂量目录的 `PhysiCell_settings.xml` 是完整的 Studio 兼容配置。Studio 可以读取空间域、BioFVM 底物、细胞定义和时间设置；自定义 BBB/细胞内动力学源码位于 `spatial_model/physicell_project/` 和 `spatial_model/cpp_core/`。修改 Studio 参数后，需要重新运行项目构建与流水线，不能只依靠 Studio 的默认示例可执行文件替代自定义 C++ 模块。

## 输出结构

```text
physicell_brain_delivery/
├── comparison/comparison.pvd       # 推荐打开：三剂量并排动画
├── dose_0p3/
│   ├── PhysiCell_settings.xml      # Studio 兼容配置
│   ├── inputs/                     # PBPK 边界、细胞和参数
│   ├── raw/                        # PhysiCell MultiCellDS 原始结果
│   ├── paraview/simulation.pvd     # 单剂量 ParaView 时间序列
│   └── mass_balance.csv
├── dose_1p0/
├── dose_3p0/
├── dose_comparison_metrics.csv
├── scientific_checks.json
└── run_metadata.json
```

VTI 场数组包括 `extracellular_AAV`、其对数、血管/组织掩膜、距血管距离和 BBB 释放率；VTP 包含细胞位置、类型、半径、18 个动力学状态、编辑率与脱靶负担。

## 解释边界

空间与时间使用 μm 和 min；AAV 数量仍是相对于中剂量脑血峰值归一化的量。当前细胞摄取和组织参数是可配置的演示先验，尚未进行体内标定，因此结果适合机制展示、参数比较和实验设计，不应解释为临床剂量预测。
