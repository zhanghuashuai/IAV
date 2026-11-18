#!/usr/bin/env python3
"""差异基因富集分析（gseapy）。

两种模式：
  --gmt FILE       本地 GSEA prerank（推荐，可离线）
  无 --gmt          Enrichr 在线富集（需要网络）
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

from utils import ensure_dir, setup_logging  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--deseq", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--padj", type=float, default=0.05)
    p.add_argument("--lfc", type=float, default=1.0)
    p.add_argument("--id-map", default="",
                   help="两列 TSV：transcript_id<TAB>symbol")
    p.add_argument("--gmt", default="")
    p.add_argument("--gene-sets", default="KEGG_2021_Human")
    p.add_argument("--organism", default="Human")
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--min-size", type=int, default=5)
    p.add_argument("--max-size", type=int, default=500)
    args = p.parse_args()

    setup_logging()
    try:
        import gseapy
    except ImportError:
        sys.exit(
            "[ERROR] 未安装 gseapy。请使用 envs/flurnaseq.yml，"
            "或 pip install gseapy。"
        )

    res = pd.read_csv(args.deseq, sep="\t", index_col=0)
    for col in ("padj", "log2FoldChange"):
        if col not in res.columns:
            sys.exit(f"[ERROR] 结果表缺少 {col} 列。")

    if args.id_map:
        mapping = pd.read_csv(
            args.id_map, sep="\t", header=None, names=["id", "symbol"]
        )
        mapping = mapping.drop_duplicates("id").set_index("id")["symbol"]
        res.index = [mapping.get(str(i), i) for i in res.index]

    sig = res[
        (res["padj"] < args.padj)
        & (res["log2FoldChange"].abs() >= args.lfc)
    ]
    genes = [
        str(g) for g in sig.index
        if str(g) and str(g).lower() != "nan"
    ]
    if not genes:
        sys.exit("[ERROR] 没有满足阈值的显著基因。")

    outdir = ensure_dir(args.outdir)

    if args.gmt:
        ranked = res["log2FoldChange"].dropna().sort_values(ascending=False)
        rnk = pd.DataFrame(
            {"gene": ranked.index.astype(str), "score": ranked.values}
        )
        pre = gseapy.prerank(
            rnk=rnk,
            gene_sets=args.gmt,
            outdir=str(outdir),
            permutation_num=args.permutations,
            seed=42,
            min_size=args.min_size,
            max_size=args.max_size,
        )
        df = pre.res2d
        score_col = "NES"
    else:
        enr = gseapy.enrichr(
            gene_list=genes,
            gene_sets=args.gene_sets,
            outdir=str(outdir),
            organism=args.organism,
        )
        df = enr.results
        score_col = "Combined Score"

    if df is None or df.empty:
        sys.exit("[ERROR] 富集分析没有返回结果。")

    df.to_csv(outdir / "enrichment.tsv", sep="\t", index=False)

    if score_col in df.columns:
        plot_df = df.copy()
        plot_df[score_col] = pd.to_numeric(
            plot_df[score_col], errors="coerce"
        )
        plot_df = (
            plot_df.dropna(subset=[score_col])
            .sort_values(score_col)
            .tail(15)
        )
        if not plot_df.empty and "Term" in plot_df.columns:
            fig, ax = plt.subplots(figsize=(8, 6))
            ax.barh(plot_df["Term"], plot_df[score_col],
                    color="#2c7fb8")
            ax.set_xlabel(score_col)
            ax.set_title("Top enriched terms")
            fig.tight_layout()
            fig.savefig(outdir / "enrichment_bar.png", dpi=200)
            plt.close(fig)

    print(f"results: {outdir / 'enrichment.tsv'}")


if __name__ == "__main__":
    main()
