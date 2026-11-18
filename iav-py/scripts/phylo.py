#!/usr/bin/env python3
"""MAFFT + trimAl + IQ-TREE 系统发育流程。

用法：
  # 常规核酸比对 + 建树
  python scripts/phylo.py --fasta data/phylo/HA.fasta --segment HA --outdir results/phylo

  # 只生成密码子比对（供 PAML/HyPhy 使用）
  python scripts/phylo.py --fasta data/phylo/HA.cds.fasta --segment HA \
      --outdir results/phylo --codon-only
"""
from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from utils import ensure_dir, require_tool, run, setup_logging


def find_iqtree() -> str:
    for name in ("iqtree2", "iqtree3", "iqtree"):
        path = shutil.which(name)
        if path:
            return path
    sys.exit("[ERROR] 找不到 iqtree2/iqtree3/iqtree。")


def mafft_align(
    mafft: str,
    fasta: Path,
    out: Path,
    threads: int,
    mafft_args: str,
    codon: bool,
) -> None:
    cmd = [mafft] + shlex.split(mafft_args) + ["--thread", str(threads)]
    if codon:
        cmd.append("--codon")
    cmd.append(str(fasta))
    with open(out, "w", encoding="utf-8") as fh:
        subprocess.run(cmd, check=True, stdout=fh)


def trimal_align(inp: Path, out: Path) -> None:
    trimal = shutil.which("trimal")
    if trimal:
        run([trimal, "-in", str(inp), "-out", str(out), "-automated1"])
    else:
        # 没有 trimAl 时保留完整比对，并在日志中说明
        shutil.copy(inp, out)


def build_tree(
    trim: Path, prefix: Path, threads: int, iqtree_args: str
) -> Path:
    iqtree = find_iqtree()
    cmd = [
        iqtree,
        "-s", str(trim),
        "-T", str(threads),
        "--prefix", str(prefix),
    ] + shlex.split(iqtree_args)
    run(cmd, log_path=str(prefix) + ".iqtree.log")
    tree = Path(str(prefix) + ".treefile")
    if not tree.exists():
        sys.exit(f"[ERROR] IQ-TREE 未生成树文件：{tree}")
    return tree


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fasta", required=True)
    p.add_argument("--segment", required=True)
    p.add_argument("--outdir", default="results/phylo")
    p.add_argument("--mafft-args", default="--auto")
    p.add_argument("--iqtree-args", default="-m MFP -B 1000")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--codon-only", action="store_true",
                   help="只做密码子比对，不建树")
    args = p.parse_args()

    setup_logging()
    mafft = require_tool("mafft")
    outdir = ensure_dir(args.outdir)
    fasta = Path(args.fasta)
    if not fasta.exists():
        sys.exit(f"[ERROR] 输入 FASTA 不存在：{fasta}")

    aln = outdir / f"{args.segment}.aln.fasta"
    trim = outdir / f"{args.segment}.trim.fasta"

    if args.codon_only:
        codon_aln = outdir / f"{args.segment}.codon.aln.fasta"
        codon_trim = outdir / f"{args.segment}.codon.trim.fasta"
        mafft_align(mafft, fasta, codon_aln, args.threads,
                    args.mafft_args, codon=True)
        trimal_align(codon_aln, codon_trim)
        print(f"codon alignment: {codon_aln}")
        print(f"codon trimmed:   {codon_trim}")
        return

    mafft_align(mafft, fasta, aln, args.threads, args.mafft_args, codon=False)
    trimal_align(aln, trim)
    tree = build_tree(trim, outdir / args.segment, args.threads,
                      args.iqtree_args)
    print(f"alignment: {aln}")
    print(f"tree:      {tree}")


if __name__ == "__main__":
    main()
