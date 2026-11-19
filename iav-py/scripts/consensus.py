#!/usr/bin/env python3
"""从双端 reads 生成 8 节段 consensus。

两种模式：
  irma : 调用 CDC IRMA（推荐，甲流专用）
  map  : minimap2 + samtools + pysam，自定义深度/一致性阈值
"""
from __future__ import annotations

import argparse
import logging
import re
import shutil
import subprocess
import sys
from collections import Counter
from pathlib import Path

from utils import (
    ensure_dir,
    read_fasta,
    run,
    run_pipeline,
    setup_logging,
    write_fasta,
    write_tsv,
)

SEGMENT_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def normalize_reference(
    ref_path: Path, segments: list[str], out_path: Path
) -> dict[str, str]:
    """读取多节段参考，按节段名重写 header 为 PB2/PB1/...。"""
    records = read_fasta(ref_path)
    mapping: dict[str, str] = {}
    for seg in segments:
        for header, seq in records.items():
            token = header.split()[0].upper()
            if (
                token == seg
                or token.startswith(seg + "|")
                or token.startswith(seg + "_")
            ):
                mapping[seg] = seq
                break
    missing = [s for s in segments if s not in mapping]
    if missing:
        sys.exit(
            "[ERROR] 参考 FASTA 缺少节段："
            + ", ".join(missing)
            + "。请让 header 为 PB2/PB1/...，或 PB2|xxx、PB2_xxx。"
        )
    write_fasta({s: mapping[s] for s in segments}, out_path)
    return mapping


def infer_segment(path: Path, segments: list[str]) -> str | None:
    stem = path.stem
    tokens = [t.upper() for t in re.split(r"[_\-.|]", stem)]
    for seg in sorted(segments, key=len, reverse=True):
        if stem.upper() == seg or seg in tokens:
            return seg
    return None


def mean_depth_from_bam(bam_path: Path, length: int) -> float | None:
    try:
        import pysam
    except ImportError:
        return None
    try:
        total = 0
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            for col in bam.pileup(truncate=True):
                total += col.nsegments
        return total / length if length else None
    except Exception:
        return None


def run_irma(
    args: argparse.Namespace, workdir: Path, segments: list[str]
) -> tuple[dict[str, str], list[dict]]:
    irma = shutil.which("IRMA") or shutil.which("irma")
    if not irma:
        sys.exit("[ERROR] 找不到 IRMA。请安装 bioconda 的 irma 包。")

    if workdir.exists():
        shutil.rmtree(workdir)
    ensure_dir(workdir.parent)

    cmd = [irma, args.irma_config, args.r1, args.r2, str(workdir)]
    log_path = str(workdir) + ".irma.log"
    try:
        run(cmd, log_path=log_path)
    except subprocess.CalledProcessError:
        sys.exit(
            f"[ERROR] IRMA 运行失败，详见日志：{log_path}\n"
            "  常见原因：--irma-config 名称不对；Illumina 双端 FASTQ 用 FLU"
            "（快速版 FLU-fast）。"
        )

    # IRMA 输出布局随版本/运行方式变化：
    #   新版单样本：workdir/A_HA_H1.fasta
    #   经典布局：  workdir/consensus/*.fasta
    files: list[Path] = []
    for d in (workdir / "consensus", workdir):
        if not d.is_dir():
            continue
        candidates = sorted(d.glob("*.fasta")) + sorted(d.glob("*.fa"))
        if candidates:
            files = candidates
            break
    found: dict[str, Path] = {}
    for f in files:
        seg = infer_segment(f, segments)
        if seg and seg not in found:
            found[seg] = f

    records: dict[str, str] = {}
    rows: list[dict] = []
    for seg in segments:
        if seg not in found:
            rows.append(
                {
                    "sample": args.sample,
                    "segment": seg,
                    "status": "missing",
                    "length": 0,
                    "n_count": 0,
                    "mean_depth": "",
                    "source": "IRMA",
                }
            )
            continue
        recs = read_fasta(found[seg])
        if not recs:
            rows.append(
                {
                    "sample": args.sample,
                    "segment": seg,
                    "status": "empty",
                    "length": 0,
                    "n_count": 0,
                    "mean_depth": "",
                    "source": "IRMA",
                }
            )
            continue
        seq = next(iter(recs.values()))
        records[f"{args.sample}|{seg}"] = seq
        bam_candidates = [
            p for p in workdir.glob(f"**/*{seg}*.bam") if p.is_file()
        ]
        depth = (
            mean_depth_from_bam(bam_candidates[0], len(seq))
            if bam_candidates
            else None
        )
        rows.append(
            {
                "sample": args.sample,
                "segment": seg,
                "status": "ok",
                "length": len(seq),
                "n_count": seq.count("N"),
                "mean_depth": f"{depth:.1f}" if depth is not None else "",
                "source": "IRMA",
            }
        )
    return records, rows


def call_consensus(
    bam_path: Path,
    ref_name: str,
    ref_seq: str,
    min_depth: int,
    threshold: float,
) -> tuple[str, list[int]]:
    import pysam

    calls = ["N"] * len(ref_seq)
    depths = [0] * len(ref_seq)
    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        for col in bam.pileup(
            ref_name,
            truncate=True,
            min_base_quality=20,
            min_mapping_quality=20,
            stepper="samtools",
        ):
            pos = col.reference_pos
            if pos >= len(ref_seq):
                continue
            bases = []
            for b in col.get_query_sequences(add_indels=False):
                if not b or b in "*#":
                    continue
                bases.append(b.upper())
            depth = len(bases)
            depths[pos] = depth
            if depth >= min_depth:
                base, count = Counter(bases).most_common(1)[0]
                if count / depth >= threshold:
                    calls[pos] = base
    return "".join(calls), depths


def run_map(
    args: argparse.Namespace, segments: list[str], workdir: Path
) -> tuple[dict[str, str], list[dict]]:
    minimap2 = shutil.which("minimap2")
    samtools = shutil.which("samtools")
    if not minimap2:
        sys.exit("[ERROR] 找不到 minimap2。")
    if not samtools:
        sys.exit("[ERROR] 找不到 samtools。")

    ensure_dir(workdir)
    ref_norm = workdir / "segments.norm.fa"
    refs = normalize_reference(Path(args.ref), segments, ref_norm)

    bam = workdir / "aligned.sorted.bam"
    log = workdir / "map.log"
    run_pipeline(
        [
            [
                minimap2,
                "-ax", "sr",
                "-t", str(args.threads),
                str(ref_norm),
                args.r1,
                args.r2,
            ],
            [samtools, "sort", "-@", str(args.threads), "-o", str(bam), "-"],
        ],
        log_path=log,
    )
    run([samtools, "index", str(bam)], log_path=log)

    records: dict[str, str] = {}
    rows: list[dict] = []
    for seg in segments:
        seq, depths = call_consensus(
            bam, seg, refs[seg], args.min_depth, args.threshold
        )
        records[f"{args.sample}|{seg}"] = seq
        mean_depth = sum(depths) / len(depths) if depths else 0.0
        rows.append(
            {
                "sample": args.sample,
                "segment": seg,
                "status": "ok",
                "length": len(seq),
                "n_count": seq.count("N"),
                "mean_depth": f"{mean_depth:.1f}",
                "source": "minimap2+pysam",
            }
        )
    return records, rows


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--sample", required=True)
    p.add_argument("--r1", required=True)
    p.add_argument("--r2", required=True)
    p.add_argument("--ref", default="",
                   help="仅 map 模式需要；IRMA 模式可留空")
    p.add_argument("--out-fasta", required=True)
    p.add_argument("--out-stats", required=True)
    p.add_argument("--assembler", choices=["irma", "map"], default="irma")
    p.add_argument("--irma-config", default="FLU")
    p.add_argument("--min-depth", type=int, default=10)
    p.add_argument("--threshold", type=float, default=0.5)
    p.add_argument("--workdir", required=True)
    p.add_argument("--segments", default="PB2,PB1,PA,HA,NP,NA,MP,NS")
    p.add_argument("--threads", type=int, default=8)
    args = p.parse_args()

    setup_logging()
    segments = [
        s.strip().upper() for s in args.segments.split(",") if s.strip()
    ]
    workdir = Path(args.workdir)
    ensure_dir(workdir.parent)

    if args.assembler == "irma":
        records, rows = run_irma(args, workdir, segments)
    else:
        if not args.ref:
            sys.exit("[ERROR] map 模式需要 --ref 提供 8 节段参考 FASTA。")
        records, rows = run_map(args, segments, workdir)

    write_fasta(records, args.out_fasta)
    write_tsv(rows, args.out_stats)

    if not records:
        Path(args.out_fasta).unlink(missing_ok=True)
        sys.exit(
            "[ERROR] 没有生成任何节段 consensus。请检查：\n"
            f"  - stats: {args.out_stats}\n"
            f"  - 组装日志: {workdir}.irma.log（IRMA）或 {workdir}/map.log（map）\n"
            "  - reads 是否与参考/预期亚型匹配、是否有足够流感 reads。"
        )
    if len(records) < len(segments):
        logging.warning(
            "只组装出 %d/%d 个节段；缺失节段请查看 stats 表。",
            len(records), len(segments),
        )


if __name__ == "__main__":
    main()
