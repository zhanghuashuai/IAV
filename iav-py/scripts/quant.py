#!/usr/bin/env python3
"""Salmon 转录组定量（宿主 + 病毒合并索引）。

如果索引不存在，会用 --transcriptome 自动构建。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from utils import ensure_dir, require_tool, run, setup_logging


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--r1", required=True)
    p.add_argument("--r2", required=True)
    p.add_argument("--index", required=True)
    p.add_argument("--transcriptome", default="")
    p.add_argument("--outdir", required=True)
    p.add_argument("--libtype", default="A")
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args()

    setup_logging()
    require_tool("salmon")
    ensure_dir(args.outdir)

    index = Path(args.index)
    if not index.exists():
        if not args.transcriptome:
            sys.exit(
                "[ERROR] Salmon 索引不存在，且没有提供 --transcriptome。"
            )
        if not Path(args.transcriptome).exists():
            sys.exit(f"[ERROR] 找不到转录组：{args.transcriptome}")
        run(
            [
                "salmon", "index",
                "-t", args.transcriptome,
                "-i", str(index),
                "-k", "31",
                "-p", str(args.threads),
            ]
        )

    run(
        [
            "salmon", "quant",
            "-i", str(index),
            "-l", args.libtype,
            "-1", args.r1,
            "-2", args.r2,
            "-p", str(args.threads),
            "-o", args.outdir,
        ]
    )

    quant = Path(args.outdir) / "quant.sf"
    if not quant.exists():
        sys.exit(f"[ERROR] Salmon 未生成 quant.sf：{quant}")
    print(f"quant: {quant}")


if __name__ == "__main__":
    main()
