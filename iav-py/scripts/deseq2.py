#!/usr/bin/env python3
"""用 pydeseq2 做差异表达，输出结果表与 PCA/volcano/heatmap。"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import seaborn as sns  # noqa: E402

from utils import ensure_dir, setup_logging  # noqa: E402


def load_counts(quant_dir: str, tx2gene: str | None) -> pd.DataFrame:
    files = sorted(Path(quant_dir).glob("*/quant.sf"))
    if not files:
        sys.exit(f"[ERROR] {quant_dir} 下没有 */quant.sf")
    series = {}
    for f in files:
        df = pd.read_csv(f, sep="\t", index_col=0)
        if "NumReads" not in df.columns:
            sys.exit(f"[ERROR] 缺少 NumReads 列：{f}")
        series[f.parent.name] = df["NumReads"]
    counts = pd.DataFrame(series).fillna(0).T
    counts.index.name = "sample"

    if tx2gene:
        mapping = pd.read_csv(
            tx2gene, sep="\t", header=None, names=["tx", "gene"]
        )
        mapping = mapping.drop_duplicates("tx").set_index("tx")["gene"]
        gene = counts.columns.map(mapping)
        keep = pd.notna(gene)
        counts = counts.loc[:, keep].T.groupby(gene[keep]).sum().T
    return counts.astype(int)


def plot_pca(
    counts: pd.DataFrame, meta: pd.DataFrame, factor: str, out: Path
) -> None:
    x = np.log2(counts.values.astype(float) + 1)
    x = (x - x.mean(axis=0)) / (x.std(axis=0) + 1e-9)
    u, s, _vt = np.linalg.svd(x, full_matrices=False)
    pcs = u[:, :2] * s[:2]
    var = (s ** 2) / np.sum(s ** 2)

    fig, ax = plt.subplots(figsize=(5.5, 5))
    groups = meta[factor].astype(str).values
    for group in sorted(set(groups)):
        mask = groups == group
        ax.scatter(pcs[mask, 0], pcs[mask, 1], label=group, s=60)
    for i, sample in enumerate(meta.index):
        ax.annotate(str(sample), (pcs[i, 0], pcs[i, 1]),
                    fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax.set_xlabel(f"PC1 ({var[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({var[1] * 100:.1f}%)")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_volcano(res: pd.DataFrame, alpha: float, lfc: float, out: Path) -> None:
    df = res.dropna(subset=["padj", "log2FoldChange"]).copy()
    df["neglog10p"] = -np.log10(df["padj"].clip(lower=1e-300))
    sig = (df["padj"] < alpha) & (df["log2FoldChange"].abs() > lfc)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(df.loc[~sig, "log2FoldChange"], df.loc[~sig, "neglog10p"],
               s=8, color="grey", alpha=0.5)
    ax.scatter(df.loc[sig, "log2FoldChange"], df.loc[sig, "neglog10p"],
               s=10, color="#c0392b", alpha=0.8)
    ax.axvline(lfc, ls="--", lw=0.8, color="black")
    ax.axvline(-lfc, ls="--", lw=0.8, color="black")
    ax.axhline(-np.log10(alpha), ls="--", lw=0.8, color="black")
    ax.set_xlabel("log2 fold change")
    ax.set_ylabel("-log10 adjusted p")
    fig.tight_layout()
    fig.savefig(out, dpi=200)
    plt.close(fig)


def plot_heatmap(
    counts: pd.DataFrame, res: pd.DataFrame, out: Path, top_n: int = 50
) -> None:
    top = res.dropna(subset=["padj"]).sort_values("padj").head(top_n).index
    if len(top) < 2:
        return
    mat = np.log2(counts[top].astype(float) + 1).T
    z = mat.sub(mat.mean(axis=1), axis=0).div(mat.std(axis=1) + 1e-9, axis=0)
    grid = sns.clustermap(
        z, cmap="RdBu_r", center=0, figsize=(8, 10), yticklabels=True
    )
    grid.fig.savefig(out, dpi=200)
    plt.close(grid.fig)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--quant-dir", required=True,
                   help="包含 */quant.sf 的目录")
    p.add_argument("--metadata", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--design", default="~condition")
    p.add_argument("--contrast", required=True,
                   help="如 condition,H5N1,mock")
    p.add_argument("--tx2gene", default="")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=1.0)
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args()

    setup_logging()
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.ds import DeseqStats
    except ImportError:
        sys.exit(
            "[ERROR] 未安装 pydeseq2。请使用 envs/flurnaseq.yml，"
            "或 pip install pydeseq2。"
        )

    outdir = ensure_dir(args.outdir)
    counts = load_counts(args.quant_dir, args.tx2gene or None)
    meta = pd.read_csv(args.metadata, sep="\t", index_col=0)

    common = counts.index.intersection(meta.index)
    if len(common) < len(counts.index):
        missing = sorted(set(counts.index) - set(meta.index))
        print(f"[WARN] 以下样本没有元数据，将被丢弃：{missing}")
    counts = counts.loc[common]
    meta = meta.loc[common]

    contrast = [x.strip() for x in args.contrast.split(",")]
    if len(contrast) != 3:
        sys.exit("[ERROR] --contrast 需要三部分：factor,test,reference")

    try:
        dds = DeseqDataSet(
            counts=counts, metadata=meta, design=args.design,
            n_cpus=args.threads,
        )
    except TypeError:
        dds = DeseqDataSet(counts=counts, metadata=meta, design=args.design)
    dds.deseq2()

    try:
        stat = DeseqStats(dds, contrast=contrast, n_cpus=args.threads)
    except TypeError:
        stat = DeseqStats(dds, contrast=contrast)
    stat.summary()
    res = stat.results_df
    res.to_csv(outdir / "deseq2_results.tsv", sep="\t")

    plot_pca(counts, meta, contrast[0], outdir / "pca.png")
    plot_volcano(res, args.alpha, args.lfc, outdir / "volcano.png")
    plot_heatmap(counts, res, outdir / "heatmap.png")

    print(f"results: {outdir / 'deseq2_results.tsv'}")
    print(f"figures: {outdir}")


if __name__ == "__main__":
    main()
