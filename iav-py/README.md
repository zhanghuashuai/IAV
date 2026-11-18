# iav-py

甲型流感病毒生信分析的 **Python-first** 项目骨架。

核心思想：**只用 Python 写脚本和组织流程**，把 IRMA、MAFFT、IQ-TREE、TreeTime、PAML、Nextclade、Salmon 这些成熟工具当作“计算内核”调用。这样既保留了行业标准的准确性，又只需要维护一套 Python 代码。

---

## 1. 设计原则

| 层次 | 做法 |
|---|---|
| 编排 | Snakemake（Python DSL）或 `scripts/*.py` 直接运行 |
| 解析 | Biopython、pandas、PyYAML、pysam |
| 分析 | 自写 Python + 调用专业 CLI |
| 统计 | pydeseq2、gseapy、scipy、statsmodels |
| 可视化 | matplotlib、seaborn |
| 报告 | Markdown / Jupyter / nbconvert |

不追求“零外部二进制”。多序列比对、ML 建树、IRMA 组装、PAML/HyPhy 这些步骤，用 Python 重写既不划算也不可靠。

---

## 2. 目录结构

```text
iav-py/
├── Snakefile
├── config/
│   ├── config.yaml          # 主配置
│   └── markers.yaml         # 耐药/宿主适应标记表
├── envs/
│   ├── flubio.yml           # 病毒基因组主环境
│   └── flurnaseq.yml        # 转录组环境
├── scripts/
│   ├── utils.py
│   ├── qc.py                # cutadapt
│   ├── consensus.py         # IRMA 或 minimap2+pysam
│   ├── nextclade_typing.py  # Nextclade 分型
│   ├── markers.py           # 功能标记扫描
│   ├── phylo.py             # MAFFT + trimAl + IQ-TREE
│   ├── timetree.py          # TreeTime
│   ├── reassortment.py      # 8 节段拓扑比较
│   ├── selection.py         # PAML codeml
│   ├── plot_tree.py         # 快速画树
│   ├── quant.py             # Salmon
│   ├── deseq2.py            # pydeseq2
│   ├── enrich.py            # gseapy
│   └── report.py            # Markdown 报告
├── notebooks/
│   ├── 01_qc_and_markers.ipynb
│   └── 02_phylo.ipynb
├── data/
│   ├── raw/
│   ├── metadata/
│   └── phylo/
├── ref/
├── results/
├── figures/
├── reports/
└── tests/
```

---

## 3. 安装

推荐 WSL2 + Ubuntu + Miniforge：

```bash
wsl --install -d Ubuntu-22.04
```

在 WSL 中：

```bash
mamba env create -f envs/flubio.yml
mamba activate flubio

# 需要做宿主转录组时
mamba env create -f envs/flurnaseq.yml
mamba activate flurnaseq
```

检查关键工具：

```bash
IRMA 2>&1 | head -n 1
nextclade --version
mafft --version
iqtree2 --version
treetime --version
codeml 2>&1 | head -n 1
salmon --version
```

---

## 4. 准备数据

### 4.1 样本 FASTQ

编辑 `config/config.yaml`：

```yaml
samples:
  sampleA:
    r1: data/raw/sampleA_R1.fastq.gz
    r2: data/raw/sampleA_R2.fastq.gz
```

### 4.2 8 节段参考（仅 map 模式需要）

使用 `assembly.assembler: map` 时，`ref/segments.fa` 的 header 必须是节段名：

```text
>PB2
...
>PB1
...
>PA
...
>HA
...
>NP
...
>NA
...
>MP
...
>NS
...
```

也支持 `>PB2|A/xxx`、`>PB2_xxx`。

### 4.3 系统发育数据

```text
data/phylo/HA.fasta          # 多条 HA 序列
data/metadata/metadata.tsv   # 至少包含 name 和 date 两列
```

### 4.4 转录组数据（可选）

```text
ref/combined.fa              # 宿主转录本 + 病毒 8 段
data/metadata/rnaseq_metadata.tsv
```

`rnaseq_metadata.tsv` 第一列是样本名，至少有一个分组列（如 `condition`）。

---

## 5. 快速开始

```bash
snakemake -n                 # 先看 DAG
snakemake -j 8               # 运行核心流程
```

核心流程：

```text
FASTQ
  -> cutadapt
  -> IRMA / minimap2+pysam consensus
  -> Nextclade 分型
  -> 功能标记扫描
  -> reports/report.md
```

输出：

```text
results/qc/                  # 修剪后 FASTQ + JSON
results/consensus/           # 合并的 8 段 consensus + stats
results/irma/<sample>/       # IRMA 原始输出
results/nextclade/<sample>.tsv
results/markers/<sample>.tsv
reports/report.md
```

### 切换组装模式

```yaml
assembly:
  assembler: map       # 改为 map 后使用 minimap2 + pysam
  min_depth: 10
  consensus_threshold: 0.5
```

`irma` 模式需要 `irma` 在 PATH 中；`map` 模式需要 `minimap2`、`samtools` 和 `pysam`。

---

## 6. 可选模块

### 6.1 系统发育树

```bash
snakemake -j 8 results/phylo/HA.treefile
```

输出：

```text
results/phylo/HA.aln.fasta
results/phylo/HA.trim.fasta
results/phylo/HA.treefile
results/phylo/HA.iqtree
```

换节段时修改 `config.yaml` 的 `phylogeny.segment` 和 `input_fasta`，或者直接运行：

```bash
python scripts/phylo.py \
  --fasta data/phylo/NA.fasta \
  --segment NA \
  --outdir results/phylo
```

### 6.2 时间树

```bash
snakemake -j 8 results/phylo/HA.timetree.nwk
```

输出：

```text
results/phylo/HA.timetree.nwk
results/phylo/HA.clock.tsv
results/phylo/molecular_clock.txt
```

### 6.3 密码子比对与选择压力

```bash
snakemake -j 8 results/phylo/HA.codon.trim.fasta
snakemake -j 8 results/selection/HA.selection.tsv
```

`selection.py` 用 Biopython 调用 PAML codeml，做 M1a/M2a 和 M7/M8 比较，并尝试解析 BEB 位点。

### 6.4 重配筛查

先对 8 个节段分别建树，并确认 `config.yaml` 的 `reassortment.trees` 路径存在：

```bash
snakemake -j 8 results/reassortment/rf_matrix.tsv
```

输出：

```text
results/reassortment/rf_matrix.tsv
results/reassortment/candidates.tsv
results/reassortment/taxa_summary.tsv
```

方法是 clade 集合的 Jaccard 距离，属于筛查；最终判断需要结合 bootstrap 支持度和流行病学证据。

### 6.5 快速画树

```bash
snakemake -j 8 figures/HA_tree.png
```

出版级美化建议用 ggtree / iTOL / toytree。

### 6.6 宿主转录组

> Snakefile 的转录组规则复用 `samples` 列表；如果转录组样本与病毒样本不同，建议直接用 `scripts/quant.py` / `deseq2.py` / `enrich.py`，或另建一份配置。

```bash
# 1) 定量
snakemake -j 8 results/rnaseq/sampleA/quant.sf

# 2) 差异表达
snakemake -j 8 results/deseq2/deseq2_results.tsv

# 3) 富集
snakemake -j 8 results/enrich/enrichment.tsv
```

或在 `flurnaseq` 环境中直接运行：

```bash
python scripts/quant.py \
  --r1 data/raw/S1_R1.fastq.gz --r2 data/raw/S1_R2.fastq.gz \
  --index ref/salmon_idx --transcriptome ref/combined.fa \
  --outdir results/rnaseq/S1

python scripts/deseq2.py \
  --quant-dir results/rnaseq \
  --metadata data/metadata/rnaseq_metadata.tsv \
  --outdir results/deseq2 \
  --design '~condition' \
  --contrast condition,H5N1,mock \
  --tx2gene ref/tx2gene.tsv

python scripts/enrich.py \
  --deseq results/deseq2/deseq2_results.tsv \
  --outdir results/enrich \
  --gmt ref/h.all.v2024.1.Hs.symbols.gmt
```

---

## 7. Python 与外部工具对照

| 环节 | 本项目脚本 | 计算内核 | 纯 Python 替代 |
|---|---|---|---|
| QC | `qc.py` | cutadapt | cutadapt 本身就是 Python |
| 组装 | `consensus.py` | IRMA 或 minimap2 | pysam/mappy 可做 consensus，但 IRMA 更稳 |
| 分型 | `nextclade_typing.py` | Nextclade | augur clades + 自定义规则 |
| 比对 | `phylo.py` | MAFFT | 无可靠替代 |
| 建树 | `phylo.py` | IQ-TREE | 无可靠 ML 替代 |
| 时间树 | `timetree.py` | TreeTime | TreeTime 本身就是 Python |
| 选择压力 | `selection.py` | PAML codeml | Bio.Phylo.PAML 已是 Python 接口 |
| 重配 | `reassortment.py` | 自写 Python + dendropy | 纯 Python |
| 定量 | `quant.py` | Salmon | HTSeq/pysam 计数（较慢） |
| 差异表达 | `deseq2.py` | pydeseq2 | 纯 Python |
| 富集 | `enrich.py` | gseapy | 纯 Python |
| 画图 | `plot_tree.py` | matplotlib + Biopython | 纯 Python |

---

## 8. 报告

```bash
snakemake -j 8 reports/report.md
```

生成 `reports/report.md`，汇总：

- 组装/覆盖度；
- Nextclade 分型；
- 功能标记；
- 分子钟；
- 重配候选；
- 正选择结果；
- `figures/` 下的图片。

---

## 9. 注意事项

1. **编号体系**：HA 的 H1/H3/H5 numbering、NA 的 N1/N2 numbering 不同；`markers.yaml` 只是模板，使用前逐条核对。
2. **Nextclade 数据集名称**会更新，先运行 `nextclade dataset list | grep -i flu`。
3. **PAML 在 Windows 原生环境容易失败**，建议 WSL2/Linux。
4. **PAML 结果解析**依赖 Biopython 版本；若解析失败，直接查看 `results/selection/paml/rst`。
5. **重配筛查**是拓扑距离，不是因果证据；需要支持度和元数据。
6. **公共 RNA-seq** 有批次效应，不要直接把不同 GSE 混在一起做 DEG。
7. **GISAID 数据**受数据使用协议约束，不得二次分发。
8. **计算预测不是实验证据**；宿主适应、耐药、抗原逃逸最终需要实验验证。

---

## 10. 下一步

- 把 `config/markers.yaml` 扩展成你关心的亚型/药物组合；
- 用 `notebooks/` 做探索性分析和图；
- 把 `Snakefile` 拆成 `rules/*.smk`，接入更多样本；
- 加 CI（GitHub Actions）做脚本语法检查和单元测试。
