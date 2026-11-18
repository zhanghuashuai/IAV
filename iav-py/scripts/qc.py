#!/usr/bin/env python3
"""双端 FASTQ 质控与接头修剪（cutadapt），输出 JSON 摘要。"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from utils import count_fastq_reads, ensure_dir, require_tool, run, setup_logging


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--r1", required=True)
    p.add_argument("--r2", required=True)
    p.add_argument("--out-r1", required=True)
    p.add_argument("--out-r2", required=True)
    p.add_argument("--report-json", required=True)
    p.add_argument("--adapter", default="")
    p.add_argument("--adapter-r2", default="")
    p.add_argument("--quality", type=int, default=20)
    p.add_argument("--min-length", type=int, default=50)
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args()

    setup_logging()
    cutadapt = require_tool("cutadapt")
    ensure_dir(Path(args.out_r1).parent)
    ensure_dir(Path(args.out_r2).parent)

    cmd = [
        cutadapt,
        "-j", str(args.threads),
        "-q", str(args.quality),
        "-m", str(args.min_length),
        "-o", args.out_r1,
        "-p", args.out_r2,
    ]
    if args.adapter:
        cmd += ["-a", args.adapter]
    if args.adapter_r2:
        cmd += ["-A", args.adapter_r2]
    cmd += [args.r1, args.r2]

    before_r1 = count_fastq_reads(args.r1)
    run(cmd, log_path=str(Path(args.report_json).with_suffix(".cutadapt.log")))
    after_r1 = count_fastq_reads(args.out_r1)
    after_r2 = count_fastq_reads(args.out_r2)

    report = {
        "r1": args.r1,
        "r2": args.r2,
        "out_r1": args.out_r1,
        "out_r2": args.out_r2,
        "reads_before_r1": before_r1,
        "reads_after_r1": after_r1,
        "reads_after_r2": after_r2,
        "retained_fraction": (
            round(after_r1 / before_r1, 4) if before_r1 else None
        ),
        "params": {
            "quality": args.quality,
            "min_length": args.min_length,
            "adapter": args.adapter,
            "adapter_r2": args.adapter_r2,
        },
        "command": cmd,
    }
    Path(args.report_json).write_text(
        json.dumps(report, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
