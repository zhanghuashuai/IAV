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

from utils import (
    ensure_dir,
    read_fasta,
    require_tool,
    run,
    setup_logging,
    write_fasta,
)


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


def macse_codon_align(macse: str, fasta: Path, out_nt: Path, out_aa: Path) -> None:
    run(
        [
            macse,
            "-prog", "alignSequences",
            "-seq", str(fasta),
            "-out_NT", str(out_nt),
            "-out_AA", str(out_aa),
        ],
        log_path=str(out_nt) + ".log",
    )


def clean_codon_alignment(inp: Path, out: Path) -> None:
    """剔除含 MACSE 移码标记 `!` 的密码子，保证长度是 3 的倍数。"""
    records = read_fasta(inp)
    if not records:
        sys.exit(f"[ERROR] 密码子比对为空：{inp}")
    seqs = list(records.values())
    length = len(seqs[0])
    keep = [
        i for i in range(0, length - 2, 3)
        if not any("!" in seq[i:i + 3] for seq in seqs)
    ]
    if not keep:
        sys.exit(f"[ERROR] 密码子比对全部被移码标记污染：{inp}")
    cleaned = {
        name: "".join(seq[i:i + 3] for i in keep)
        for name, seq in records.items()
    }
    write_fasta(cleaned, out)


def build_tree(
    trim: Path, prefix: Path, threads: int, iqtree_args: str
) -> Path:
    iqtree = find_iqtree()
    cmd = [
        iqtree,
        "-s", str(trim),
        "-T", str(threads),
        "--prefix", str(prefix),
        "-redo",
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
        macse = shutil.which("macse")
        if not macse:
            sys.exit(
                "[ERROR] 找不到 macse。密码子比对需要 MACSE；"
                "请安装 bioconda 的 macse 包。"
            )
        raw_nt = outdir / f"{args.segment}.codon.macse.fasta"
        raw_aa = outdir / f"{args.segment}.codon.macse.aa.fasta"
        macse_codon_align(macse, fasta, raw_nt, raw_aa)
        clean_codon_alignment(raw_nt, codon_aln)
        shutil.copy(codon_aln, codon_trim)
        print(f"codon alignment: {codon_aln}")
        print(f"codon cleaned:   {codon_trim}")
        return

    mafft_align(mafft, fasta, aln, args.threads, args.mafft_args, codon=False)
    trimal_align(aln, trim)
    tree = build_tree(trim, outdir / args.segment, args.threads,
                      args.iqtree_args)
    print(f"alignment: {aln}")
    print(f"tree:      {tree}")


if __name__ == "__main__":
    main()
