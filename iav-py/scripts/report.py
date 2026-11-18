#!/usr/bin/env python3
"""汇总 results/ 下的表格与 figures/ 下的图片，生成 Markdown 报告。"""
from __future__ import annotations

import argparse
import os
from datetime import datetime
from pathlib import Path

import pandas as pd

from utils import ensure_dir, setup_logging


def df_to_md(df: pd.DataFrame, max_rows: int = 15) -> str:
    if df is None or df.empty:
        return "_无数据_\n"
    df = df.head(max_rows)
    cols = [str(c) for c in df.columns]
    lines = [
        "| " + " | ".join(cols) + " |",
        "|" + "|".join(["---"] * len(cols)) + "|",
    ]
    for _, row in df.iterrows():
        lines.append(
            "| " + " | ".join(str(v) for v in row.values) + " |"
        )
    return "\n".join(lines) + "\n"


def read_pattern(results: Path, pattern: str) -> pd.DataFrame:
    frames = []
    for path in sorted(results.glob(pattern)):
        try:
            df = pd.read_csv(path, sep="\t")
        except Exception:  # noqa: BLE001
            continue
        if not df.empty:
            frames.append(df)
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--results", default="results")
    p.add_argument("--figures", default="figures")
    p.add_argument("--out", required=True)
    args = p.parse_args()

    setup_logging()
    results = Path(args.results)
    figures = Path(args.figures)

    stats = read_pattern(results, "consensus/*.stats.tsv")
    nextclade = read_pattern(results, "nextclade/*.tsv")
    markers = read_pattern(results, "markers/*.tsv")
    clock = read_pattern(results, "phylo/*.clock.tsv")
    candidates = read_pattern(results, "reassortment/candidates.tsv")
    selection = read_pattern(results, "selection/*.selection.tsv")

    found_markers = (
        markers[markers["found"] == "yes"]
        if "found" in markers.columns else pd.DataFrame()
    )

    lines = [
        "# 甲型流感病毒生信分析报告",
        "",
        f"生成时间：{datetime.now():%Y-%m-%d %H:%M}",
        "",
        "> 本报告由 `scripts/report.py` 自动汇总。",
        "> 所有位点与功能解释需结合文献、结构和实验验证。",
        "",
        "## 1. 组装与 QC",
        "",
        df_to_md(stats),
        "",
        "## 2. Nextclade 分型",
        "",
        df_to_md(nextclade, max_rows=20),
        "",
        "## 3. 功能标记",
        "",
        "### 3.1 检出的标记",
        "",
        df_to_md(found_markers, max_rows=30),
        "",
        "### 3.2 全部标记矩阵",
        "",
        df_to_md(markers, max_rows=30),
        "",
        "## 4. 系统发育与分子钟",
        "",
        df_to_md(clock),
        "",
        "## 5. 重配筛查",
        "",
        df_to_md(candidates),
        "",
        "## 6. 选择压力",
        "",
        df_to_md(selection),
        "",
        "## 7. 图",
        "",
    ]

    figure_files = sorted(figures.glob("*.png")) if figures.exists() else []
    if figure_files:
        for fig in figure_files:
            rel = Path(os.path.relpath(fig, out.parent))
            lines += [f"### {fig.name}", "", f"![{fig.stem}]({rel.as_posix()})", ""]
    else:
        lines += ["_暂无图片_", ""]

    lines += [
        "## 8. 方法与局限",
        "",
        "- 组装：IRMA 或 minimap2+pysam；深部变异阈值需 duplicates/UMI 支持。",
        "- 分型：Nextclade 数据集版本会影响 clade 与突变注释。",
        "- 标记：编号体系随亚型不同，必须逐条核对。",
        "- 重配：拓扑距离只是筛查，需要支持度和流行病学证据。",
        "- 选择压力：dN/dS 信号不等于功能适应。",
        "",
    ]

    out = Path(args.out)
    ensure_dir(out.parent)
    out.write_text("\n".join(lines), encoding="utf-8")
    print(f"report: {out}")


if __name__ == "__main__":
    main()
