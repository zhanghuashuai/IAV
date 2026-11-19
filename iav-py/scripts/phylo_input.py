#!/usr/bin/env python3
"""汇总系统发育输入 FASTA。

来源（按优先级去重）：
  1. 用户提供的多序列 FASTA（--user-fasta）
  2. results/consensus/*.fasta 中目标节段的 consensus
  3. 参考 FASTA（--reference-fasta）中目标节段的一条序列

不足 --min-sequences 时仍写出文件，但以非零码退出并给出明确提示。
"""
from __future__ import annotations

import argparse
import glob
import re
import sys
from pathlib import Path

from utils import read_fasta, setup_logging, write_fasta

KNOWN_SEGMENTS = {"PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"}
SEGMENT_BY_NUMBER = {
    "1": "PB2", "2": "PB1", "3": "PA", "4": "HA",
    "5": "NP", "6": "NA", "7": "MP", "8": "NS",
}
SEGMENT_NAMES = {
    "PB2": ("PB2", "polymerase PB2"),
    "PB1": ("PB1", "polymerase PB1"),
    "PA": ("PA", "polymerase PA"),
    "HA": ("HA", "hemagglutinin"),
    "NP": ("NP", "nucleocapsid"),
    "NA": ("NA", "neuraminidase"),
    "MP": ("MP", "matrix protein"),
    "NS": ("NS", "nonstructural", "nuclear export"),
}


def segment_of(header: str, default: str) -> str:
    # 先找裸节段名（支持 sample|HA、A_HA_H1、(HA) 等）
    tokens = [
        re.sub(r"[^A-Za-z0-9]", "", token).upper()
        for token in re.split(r"[|_\-\s.]+", header)
        if token
    ]
    for token in tokens:
        if token in KNOWN_SEGMENTS:
            return token
    # GenBank 风格：segment 4 hemagglutinin (HA)
    match = re.search(r"\bsegment\s+([1-8])\b", header, re.IGNORECASE)
    if match:
        return SEGMENT_BY_NUMBER[match.group(1)]
    upper = header.upper()
    for segment, names in SEGMENT_NAMES.items():
        for name in names:
            if name.upper() in upper:
                return segment
    return default.upper()


def record_name(header: str) -> str:
    if "|" in header:
        return header.split("|", 1)[0].strip()
    return header.split()[0]


def add_records(
    out: dict[str, str],
    path: Path,
    segment: str,
    rename: bool = False,
    limit: int | None = None,
) -> int:
    if not path.exists():
        return 0
    kept = 0
    for header, seq in read_fasta(path).items():
        if segment_of(header, segment) != segment:
            continue
        name = record_name(header) if rename else header
        if name in out:
            continue
        out[name] = seq
        kept += 1
        if limit is not None and kept >= limit:
            break
    return kept


def add_reference(
    out: dict[str, str], path: Path, segment: str, name: str | None = None
) -> int:
    if not path.exists():
        return 0
    for header, seq in read_fasta(path).items():
        if segment_of(header, segment) == segment:
            out.setdefault(name or f"reference_{segment}", seq)
            return 1
    return 0


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--segment", required=True)
    p.add_argument("--user-fasta", default="")
    p.add_argument("--consensus-glob", default="results/consensus/*.fasta")
    p.add_argument("--reference-fasta", default="")
    p.add_argument("--dataset-reference", default="",
                   help="Nextclade 数据集参考，作为外群")
    p.add_argument("--out", required=True)
    p.add_argument("--min-sequences", type=int, default=3)
    args = p.parse_args()

    setup_logging()
    segment = args.segment.upper()
    out: dict[str, str] = {}
    sources: list[str] = []

    if args.user_fasta:
        n = add_records(out, Path(args.user_fasta), segment)
        if n:
            sources.append(f"user={args.user_fasta}:{n}")

    n_consensus = 0
    for path in sorted(glob.glob(args.consensus_glob)):
        n_consensus += add_records(out, Path(path), segment, rename=True)
    if n_consensus:
        sources.append(f"consensus:{n_consensus}")

    if args.reference_fasta:
        n = add_reference(out, Path(args.reference_fasta), segment)
        if n:
            sources.append(f"reference={args.reference_fasta}:{n}")

    if args.dataset_reference:
        n = add_reference(
            out,
            Path(args.dataset_reference),
            segment,
            name=f"nextclade_reference_{segment}",
        )
        if n:
            sources.append(f"dataset_reference={args.dataset_reference}:{n}")

    write_fasta(out, args.out)
    total = len(out)
    print(
        f"{segment}: {total} sequences -> {args.out} "
        f"({', '.join(sources) or 'no source'})"
    )

    if total < args.min_sequences:
        sys.exit(
            f"[ERROR] {segment} 只有 {total} 条序列，至少需要 "
            f"{args.min_sequences} 条才能建树。\n"
            f"  已尝试：{', '.join(sources) or '无可用输入'}。\n"
            f"  请把参考序列放到 data/phylo/{segment}.fasta，"
            "或修改 config.yaml 的 phylogeny.input_fasta；"
            "也可以增加 samples 后重跑 consensus。\n"
            f"  已写出当前可用序列：{args.out}（可用于检查，但不足以建树）。"
        )


if __name__ == "__main__":
    main()
