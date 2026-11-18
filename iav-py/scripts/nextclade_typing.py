#!/usr/bin/env python3
"""对每个节段运行 Nextclade，并合并成一张样本级 TSV。"""
from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import sys
from pathlib import Path

from utils import (
    ensure_dir,
    parse_pairs,
    read_fasta,
    setup_logging,
    write_fasta,
    write_tsv,
)


def split_by_segment(records: dict[str, str]) -> dict[str, list[tuple[str, str]]]:
    out: dict[str, list[tuple[str, str]]] = {}
    for header, seq in records.items():
        if "|" in header:
            _sample, seg = header.split("|", 1)
        else:
            parts = header.replace("_", "|").split("|")
            seg = parts[-1]
        out.setdefault(seg, []).append((header, seq))
    return out


def run_nextclade(
    fasta_path: Path, dataset: str, outdir: Path, tsv_path: Path
) -> None:
    ensure_dir(outdir)
    base = [
        "nextclade", "run",
        "--input-fasta", str(fasta_path),
        "--output-dir", str(outdir),
        "--dataset-name", dataset,
    ]
    proc = subprocess.run(
        base + ["--output-tsv", str(tsv_path)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0:
        # 兼容没有 --output-tsv 的旧版本
        proc = subprocess.run(base, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr[-2000:])
    if not tsv_path.exists():
        fallback = outdir / "nextclade.tsv"
        if fallback.exists():
            shutil.copy(fallback, tsv_path)
    if not tsv_path.exists():
        raise FileNotFoundError(f"Nextclade 未生成 TSV：{tsv_path}")


def parse_first_row(path: Path) -> dict:
    with open(path, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    return rows[0] if rows else {}


def get(row: dict, *keys: str) -> str:
    for key in keys:
        value = row.get(key)
        if value not in (None, ""):
            return str(value)
    return ""


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--consensus", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--dataset-map", required=True,
                   help="逗号分隔，如 HA=flu_h1n1pdm_ha,NA=flu_h1n1pdm_na")
    p.add_argument("--outdir", required=True)
    args = p.parse_args()

    setup_logging()
    datasets = parse_pairs(args.dataset_map)
    records = read_fasta(args.consensus)
    by_seg = split_by_segment(records)

    out_rows: list[dict] = []
    for seg in sorted(by_seg, key=lambda s: list(datasets).index(s)
                      if s in datasets else 999):
        if seg not in datasets:
            logging.warning("节段 %s 没有配置 Nextclade 数据集，跳过。", seg)
            continue

        seg_dir = ensure_dir(Path(args.outdir) / seg)
        fasta = seg_dir / f"{seg}.fasta"
        write_fasta(dict(by_seg[seg]), fasta)
        tsv = seg_dir / "nextclade.tsv"

        status = "ok"
        row: dict = {}
        try:
            run_nextclade(fasta, datasets[seg], seg_dir / "nextclade_out", tsv)
            row = parse_first_row(tsv)
        except Exception as exc:  # noqa: BLE001
            logging.error("Nextclade 失败：%s / %s: %s", seg, datasets[seg], exc)
            status = "failed"

        sample = by_seg[seg][0][0].split("|", 1)[0]
        out_rows.append(
            {
                "sample": sample,
                "segment": seg,
                "dataset": datasets[seg],
                "status": status,
                "clade": get(row, "clade"),
                "subclade": get(row, "subclade"),
                "qc_status": get(row, "qc.overallStatus", "qc.overallStatus"),
                "aa_substitutions": get(row, "aaSubstitutions"),
                "aa_deletions": get(row, "aaDeletions"),
                "total_substitutions": get(row, "totalSubstitutions"),
                "total_deletions": get(row, "totalDeletions"),
                "total_frame_shifts": get(row, "totalFrameShifts"),
                "total_missing": get(row, "totalMissing"),
                "coverage": get(row, "coverage"),
            }
        )

    if not out_rows:
        sys.exit("[ERROR] 没有产生任何 Nextclade 结果，请检查 dataset-map。")
    write_tsv(out_rows, args.out)
    logging.info("写出 %d 个节段的 Nextclade 结果：%s", len(out_rows), args.out)


if __name__ == "__main__":
    main()
