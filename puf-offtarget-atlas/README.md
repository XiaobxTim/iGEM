# PUF-OffTarget Atlas

PUF-OffTarget Atlas 用于在 GENCODE 转录本序列中搜索**潜在的** PUF 识别位点，并综合序列相似性、转录本注释、可选的 GTEx 风格组织表达数据、可选的 RNAplfold 结构可及性以及谨慎表述的潜在后果，对候选位点进行优先级排序。

本工具不会声称预测到的结合或编辑事件真实发生。所有综合风险分数都是 **heuristic prioritization score（启发式优先级评分）**，不是经过实验校准的脱靶概率。

## 生物学假设

- 输入为 PUF 所识别的 RNA 序列，方向为 5′→3′，长度必须为 8–12 nt。
- 程序内部将 `U` 统一转换为 `T`。结果展示时仍使用 RNA 字母 `U`。
- GENCODE transcript FASTA 已经是转录本 5′→3′方向，因此默认只搜索输入 motif 本身。
- MVP 使用固定长度、仅替换型的 Hamming distance，允许 0–3 个 mismatch，不处理插入或缺失。
- reverse complement 搜索必须显式启用，并且每条结果都会记录搜索方向。
- `binding_only` 模式只描述潜在结合和占位影响，不代表发生碱基编辑。
- APOBEC C→U 和 ADAR A→I（测序及翻译读取为 G）在用户提供的编辑窗口内逐个碱基模拟；MVP 不组合模拟多个碱基同时编辑。

## 安装

需要 Python 3.11 或更高版本。

在现有 `xbx_env` 环境中安装：

```bash
cd /Users/benxiang/Desktop/iGEM/puf-offtarget-atlas
conda activate xbx_env
python -m pip install -e ".[dev]"
pufscan --help
```

也可以创建包含 ViennaRNA 的独立 Conda 环境：

```bash
conda env create -f environment.yml
conda activate puf-offtarget-atlas
```

RNAplfold 是可选依赖。在 macOS 或 Linux 上可通过 Bioconda 安装：

```bash
conda install -c bioconda viennarna
```

如果系统中没有 RNAplfold，扫描仍会继续：结构相关字段为 `NA`，报告中给出警告，综合评分会删除结构组件的权重并对剩余权重重新归一化。程序不会用随机值填补缺失结构数据。

## 快速开始

### 启动 Streamlit 页面

```bash
cd /Users/benxiang/Desktop/iGEM/puf-offtarget-atlas
conda activate xbx_env
streamlit run app/streamlit_app.py
```

浏览器默认访问 `http://localhost:8501`。页面预填了合成测试数据，可输入 `AACGUCUAUA` 并点击 **Run analysis** 完成一次端到端分析。

### 使用合成数据运行 CLI

```bash
pufscan scan \
  --query AACGUCUAUA \
  --max-mismatches 1 \
  --gencode-fasta tests/data/synthetic.fa \
  --gencode-gtf tests/data/synthetic.gtf \
  --gtex-expression tests/data/expression.tsv \
  --target-tissue Liver \
  --mode binding_only \
  --no-structure \
  --output results
```

每次运行会生成独立的时间戳目录，例如：

```text
results/AACGUCUAUA_20260721_120000/
```

示例配置也可以通过 YAML 启动：

```bash
pufscan scan --config examples/example_config.yaml --no-structure
```

显式提供的 CLI 参数会覆盖 YAML 中的 query 和路径配置。

## GENCODE Release 50 数据准备

默认使用 GENCODE Human Release 50 和 GRCh38.p14。下载命令为：

```bash
pufscan download-gencode \
  --release 50 \
  --output data/gencode_v50
```

下载器会获取 transcript FASTA、标准 comprehensive GTF，以及与全部 transcript 区域匹配的 all-regions comprehensive GTF。随后执行：

```bash
pufscan prepare-gencode \
  --fasta data/gencode_v50/gencode.v50.transcripts.fa.gz \
  --gtf data/gencode_v50/gencode.v50.chr_patch_hapl_scaff.annotation.gtf.gz \
  --output data/gencode_v50/prepared
```

下载过程使用 `.part` 临时文件和 HTTP Range 支持续传，已有文件会跳过，并记录 URL、文件大小和 SHA-256。预处理会生成：

- `transcripts.fa` 及其 `.fai` 索引；
- GTF 注释 Parquet；
- transcript coordinate mapping 数据；
- 输入文件校验清单。

使用 all-regions GTF 能尽可能覆盖 transcript FASTA 中的替代区域、patch 和 haplotype 转录本。只分析主染色体时也可改用 `gencode.v50.annotation.gtf.gz`。

## GTEx 表达数据准备

表达模块接受宽表格式的 TSV、CSV、gzip 或 Parquet 文件。第一列为 `gene_id`，其余列为不同组织的 TPM：

```text
gene_id       Liver    Brain_Cortex    Whole_Blood
ENSG...       12.3     0.8             4.1
```

准备命令：

```bash
pufscan prepare-gtex \
  --input GTEx_gene_median_tpm.tsv.gz \
  --output data/gtex/gene_median_tpm.parquet
```

该模块不依赖固定的 GTEx 下载地址。处理时会保留原始 Ensembl ID，同时去掉版本号用于跨数据库匹配；重复的稳定 ID 会被拒绝。表达数据库中缺失的基因保持为 `NA`，不会错误填成 0。

指定多个目标组织时，程序使用所选组织中的最大 TPM 作为 `expression_reference`；未指定目标组织时，使用全部组织中的最大 TPM。用于归一化的 P99 从完整预处理表达矩阵计算。

## 使用真实数据扫描

### Binding-only 模式

```bash
pufscan scan \
  --query AACGUCUAUA \
  --max-mismatches 2 \
  --gencode-fasta data/gencode_v50/prepared/transcripts.fa \
  --gencode-gtf data/gencode_v50/prepared/annotation.parquet \
  --gtex-expression data/gtex/gene_median_tpm.parquet \
  --target-tissue Liver \
  --mode binding_only \
  --structure \
  --structure-top-n 2000 \
  --output results
```

Binding-only 结果只表示潜在 PUF 结合候选，可用于讨论以下可能性：RNA 稳定性、翻译、RNA 加工或剪接受到影响，以及与内源 RNA-binding protein 发生竞争。所有这些解释都需要实验验证。

### Editor-fusion 模式

APOBEC 示例：

```bash
pufscan scan \
  --query AACGUCUAUA \
  --max-mismatches 2 \
  --gencode-fasta data/gencode_v50/prepared/transcripts.fa \
  --gencode-gtf data/gencode_v50/prepared/annotation.parquet \
  --gtex-expression data/gtex/gene_median_tpm.parquet \
  --target-tissue Liver \
  --mode editor_fusion \
  --editor APOBEC_C2U \
  --editing-window=-15:10 \
  --structure \
  --output results
```

ADAR 示例只需将编辑器改为：

```bash
--editor ADAR_A2I --editing-window=-15:10
```

选择 editor-fusion 模式时必须提供编辑窗口，否则程序会报错。窗口坐标相对于 motif 起点，以 0 为基准，并且两端都包含；超出转录本边界的部分会自动截断。输出中的候选碱基称为 **potential editable base**，不代表真实编辑位点。

### 可选 reverse complement 搜索

```bash
pufscan scan ... --search-reverse-complement
```

这不是默认行为。启用后请根据结果中的 `search_orientation` 区分 motif 和 reverse complement 命中。

## 自定义位置权重

创建一个 CSV 文件，每个 motif 位置一行：

```csv
position,weight
1,1.0
2,1.2
3,1.5
```

扫描时添加：

```bash
--position-weight-file weights.csv
```

默认序列评分公式为：

```text
sequence_score = 1 - sum(weight_i * mismatch_i) / sum(weight_i)
```

代码中预留了 `SequenceScoringModel` 和 `SubstitutionMatrixScorer` 接口，用于接入完整的 position-specific substitution matrix。当前项目不包含经过实验验证的通用 PUF recognition matrix，因此统一 mismatch 模型是一个重要局限。

## 综合风险优先级评分

评分组件包括：

- `sequence_score`：序列匹配评分；
- `accessibility_score`：RNA 局部单链可及性；
- `expression_score`：组织表达归一化评分；
- `consequence_score`：基于候选区域或模拟编辑后果的排序参数。

所有可用组件会被限制在 `[0.001, 1]`。缺失组件及其权重会被删除，剩余权重重新归一化，然后使用加权几何平均：

```text
risk_score = 100 × product(component_i ** normalized_weight_i)
```

默认权重为：序列 0.50、结构可及性 0.20、表达 0.20、后果 0.10。`configs/default.yaml` 中的 consequence score 只是初始排序参数，不代表真实生物效应强度。

风险优先级标签为：

| risk score | risk priority |
|---:|---|
| `[0, 25)` | Low priority |
| `[25, 50)` | Moderate priority |
| `[50, 75)` | High priority |
| `[75, 100]` | Very high priority |

这些标签必须解释为候选验证顺序，不能解释为 confirmed off-target 或真实脱靶概率。

## 输出文件与坐标约定

每次运行的结果目录包含：

```text
run_metadata.json
all_transcript_hits.parquet
all_transcript_hits.tsv.gz
unique_genomic_loci.tsv
potential_editing_events.tsv
top_hits.tsv
candidates.bed
candidates.bed12
report.html
summary.json
```

主要结果包括序列比对、mismatch、完整和去版本号 Ensembl ID、基因与转录本元数据、GENCODE tag、基因组 blocks、外显子连接与 UTR/CDS 区域、原始及归一化组织表达、RNAplfold 概率与 opening energy、potential editing events、后果证据等级、缺失特征、警告和综合 risk priority。

坐标约定：

- 程序内部坐标为 0-based、左闭右开；
- 结果表中的 transcript 坐标和 `genomic_blocks` 为 1-based、两端包含；
- BED 和 BED12 为标准的 0-based、左闭右开坐标；
- 跨 exon junction 的候选会保留多个 genomic blocks，并输出相应 BED12 记录。

## HTML 报告

`report.html` 是独立 HTML 文件，包含运行摘要、mismatch 分布、风险分数分布、Top genes、组织表达热图、转录本区域示意图、RNA 可及性图以及证据和局限说明。

也可以针对已有结果重新生成报告：

```bash
pufscan report --run-directory results/AACGUCUAUA_YYYYMMDD_HHMMSS
```

## Streamlit 功能

```bash
streamlit run app/streamlit_app.py
```

页面支持：

- 配置 query、mismatch、GENCODE、GTEx、目标组织、模式、编辑窗口和结构分析；
- 按 gene、transcript、region、mismatch 和 risk score 筛选或排序；
- 查看所选候选的 alignment、表达、结构可及性和潜在编辑事件；
- 下载 CSV、TSV、Parquet 和 HTML 报告；
- 在 GTEx 或 RNAplfold 缺失时继续运行，仅将对应字段标记为不可用。

## 环境检查

使用合成数据检查当前环境：

```bash
pufscan doctor \
  --gencode-fasta tests/data/synthetic.fa \
  --gencode-gtf tests/data/synthetic.gtf \
  --gtex-expression tests/data/expression.tsv \
  --output results
```

`doctor` 会检查 Python 版本、GENCODE 文件、FASTA/GTF 索引、GTEx 文件、RNAplfold 是否可用，以及输出目录权限。

## 结果解释与证据等级

报告将潜在后果分成三个等级：

- **Level 1**：由转录本坐标和标准遗传密码直接计算，例如模拟编辑后的 codon 变化；
- **Level 2**：由组织表达和预测结构可及性支持的候选风险排序；
- **Level 3**：基于可靠基因功能元数据推测的潜在系统影响，必须实验验证。

固定限制说明：

- Sequence similarity does not prove PUF binding.
- Predicted accessibility does not represent in-cell RNA structure.
- Bulk tissue expression does not guarantee co-localization with the PUF editor.
- Predicted editing consequences require experimental validation.
- The integrated risk score is a prioritization score, not a calibrated probability.

其他已知局限：

- MVP 不处理插入和缺失；
- 不组合模拟多个碱基同时编辑；
- 尚未整合细胞类型表达、亚细胞定位、eCLIP、miRNA/RBP overlap 或 RNA modification 数据；
- RNAplfold 预测不能代替细胞内真实 RNA 结构测量；
- bulk tissue TPM 不能证明 PUF editor 与目标 RNA 在同一细胞或亚细胞区域共定位；
- alternative/haplotype loci 需要 all-regions GTF，并且不一定能映射到主染色体位点。

## 使用真实 RNA-seq 数据校准

若要把当前启发式排序扩展为校准概率模型，应收集具有生物学重复的处理组和对照组 RNA-seq，并使用经过验证的流程检测编辑事件。建议：

1. 在 transcript position 层面建立实验阳性和阴性标签；
2. 避免从同一基因或高度相似转录本同时抽样到训练集与测试集；
3. 加入表达、可及性、局部序列上下文、定位和实验批次等特征；
4. 按基因或染色体划分训练集和测试集，降低数据泄漏；
5. 在独立数据上进行概率校准、外部验证和不确定性评估；
6. 不要把当前 heuristic risk score 当作监督学习的真实标签。

## 批量比较设计

项目为多个候选 PUF 设计预留了命令：

```bash
pufscan compare-designs \
  --queries candidate_pufs.csv \
  --gencode-fasta data/gencode_v50/prepared/transcripts.fa \
  --gencode-gtf data/gencode_v50/prepared/annotation.parquet \
  --output comparison_results
```

该接口用于汇总不同设计的 exact、1-mismatch、2-mismatch 命中、高优先级基因和 CDS 命中、目标组织暴露、最大风险、Top-10 平均风险以及 specificity gap。具体可用参数请以 `pufscan compare-designs --help` 为准。

## 配置与扩展

默认配置位于 `configs/default.yaml`，示例配置位于 `examples/example_config.yaml`。路径、结构参数、风险权重和后果排序参数均可通过配置文件管理。

代码已为以下扩展保留接口：PUF position-specific recognition matrix、插入/缺失近似匹配、亚细胞定位、单细胞或空间转录组、RNA modification、RBP/miRNA/eCLIP 注释、实验标签、机器学习概率校准和 RNA language model embedding。这些功能不属于当前 MVP，不能从现有输出中假定已经实现。

## 可复现性

每次运行都会记录：

- 软件版本；
- GENCODE release；
- 解析后的运行参数；
- 输入文件路径；
- 生成日期和运行时间；
- 搜索速度、扫描转录本数、总碱基数和候选数量；
- 可获得时的结构计算时间和内存信息。

下载和预处理清单记录文件大小及校验值，以便复核数据版本。

## 测试与代码质量

```bash
pytest -q
ruff check .
mypy src
```

项目包含小型合成 FASTA、GTF 和表达矩阵，因此测试不需要下载完整人类转录组，也不会伪造真实分析结果。

## 引用与数据来源

使用本工具开展分析时，请引用所使用的数据版本和软件：

- [GENCODE Human](https://www.gencodegenes.org/human/)
- [GTEx Portal](https://gtexportal.org/home/)
- [ViennaRNA Package / RNAplfold](https://www.tbi.univie.ac.at/RNA/)

报告或论文中还应记录具体 GENCODE release、GTEx 数据版本、ViennaRNA 版本、PUF query、允许的 mismatch、编辑器与编辑窗口，以及本项目的软件版本。
