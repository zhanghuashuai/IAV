#!/usr/bin/env python3
"""用 PAML codeml（通过 Biopython 调用）做位点模型选择压力分析。

输出：
  <prefix>.selection.tsv   模型比较（M1a vs M2a、M7 vs M8）
  <prefix>.beb_sites.tsv   BEB 正选择位点（best effort）
  paml/mlc、paml/rst       原始输出
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

from Bio import Phylo
from Bio.Phylo.PAML import codeml
from scipy.stats import chi2

from utils import ensure_dir, require_tool, setup_logging, write_tsv


def model_result(results: dict, model: int) -> dict:
    ns = results.get("NSsites", {})
    if model in ns:
        return ns[model] or {}
    if str(model) in ns:
        return ns[str(model)] or {}
    return {}


def get_lnl(results: dict, model: int) -> float | None:
    data = model_result(results, model)
    value = data.get("lnL")
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def get_omega(results: dict, model: int) -> float | None:
    data = model_result(results, model)
    params = data.get("parameters", {})
    if not isinstance(params, dict):
        return None
    value = params.get("omega")
    if isinstance(value, (list, tuple)):
        nums = []
        for item in value:
            try:
                nums.append(float(item))
            except (TypeError, ValueError):
                pass
        return max(nums) if nums else None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def parse_beb(rst_path: Path) -> list[dict]:
    rows: list[dict] = []
    if not rst_path.exists():
        return rows
    text = rst_path.read_text(encoding="utf-8", errors="replace")
    marker = "Bayes Empirical Bayes (BEB) analysis"
    if marker not in text:
        return rows
    sections = text.split(marker)[1:]
    for idx, section in enumerate(sections, start=1):
        for line in section.splitlines():
            m = re.match(r"\s*(\d+)\s+([A-Z])\s+([0-9.]+)(\*{0,2})", line)
            if not m:
                continue
            rows.append(
                {
                    "beb_section": idx,
                    "site": int(m.group(1)),
                    "aa": m.group(2),
                    "prob": float(m.group(3)),
                    "significance": m.group(4),
                }
            )
    return rows


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--aln", required=True, help="密码子比对 FASTA")
    p.add_argument("--tree", required=True)
    p.add_argument("--outdir", default="results/selection")
    p.add_argument("--prefix", required=True)
    args = p.parse_args()

    setup_logging()
    require_tool("codeml")

    aln = Path(args.aln)
    tree = Path(args.tree)
    if not aln.exists():
        sys.exit(f"[ERROR] 找不到比对：{aln}")
    if not tree.exists():
        sys.exit(f"[ERROR] 找不到树：{tree}")

    outdir = ensure_dir(args.outdir)
    workdir = ensure_dir(outdir / "paml")
    mlc = workdir / "mlc"
    rst = workdir / "rst"

    # IQ-TREE 的 .treefile 可能带 bootstrap 内部节点标签；
    # PAML 对某些节点标签兼容性不好，先清洗成纯 Newick。
    phy = Phylo.read(str(tree), "newick")
    for clade in phy.get_nonterminals():
        clade.name = None
    clean_tree = workdir / "tree.nwk"
    Phylo.write(phy, str(clean_tree), "newick")

    cml = codeml.Codeml()
    cml.alignment = str(aln)
    cml.tree = str(clean_tree)
    cml.out_file = str(mlc)
    cml.working_dir = str(workdir)
    cml.set_options(
        seqtype=1,
        CodonFreq=2,
        model=0,
        NSsites=[0, 1, 2, 7, 8],
        fix_omega=0,
        omega=1,
        cleandata=1,
        RateAncestor=1,
    )
    results = cml.run(verbose=True)
    if not results:
        sys.exit("[ERROR] codeml 没有返回可解析结果，请检查 paml/ 目录。")

    comparisons = [
        ("M1a_vs_M2a", 1, 2),
        ("M7_vs_M8", 7, 8),
    ]
    rows = []
    for name, null_model, alt_model in comparisons:
        lnl_null = get_lnl(results, null_model)
        lnl_alt = get_lnl(results, alt_model)
        if lnl_null is None or lnl_alt is None:
            continue
        lrt = max(0.0, 2.0 * (lnl_alt - lnl_null))
        p_value = float(chi2.sf(lrt, df=2))
        rows.append(
            {
                "comparison": name,
                "null_model": null_model,
                "alt_model": alt_model,
                "lnL_null": f"{lnl_null:.4f}",
                "lnL_alt": f"{lnl_alt:.4f}",
                "LRT": f"{lrt:.4f}",
                "df": 2,
                "p_value": f"{p_value:.4g}",
                "omega_alt": (
                    f"{get_omega(results, alt_model):.3f}"
                    if get_omega(results, alt_model) is not None else ""
                ),
            }
        )

    if not rows:
        sys.exit("[ERROR] 未能提取 M1a/M2a/M7/M8 的 lnL，请检查 PAML 输出。")

    write_tsv(rows, outdir / f"{args.prefix}.selection.tsv")

    beb = parse_beb(rst)
    if beb:
        write_tsv(beb, outdir / f"{args.prefix}.beb_sites.tsv")

    print(f"selection: {outdir / (args.prefix + '.selection.tsv')}")
    if beb:
        print(f"BEB sites: {outdir / (args.prefix + '.beb_sites.tsv')}")


if __name__ == "__main__":
    main()
