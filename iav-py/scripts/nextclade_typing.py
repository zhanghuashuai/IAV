#!/usr/bin/env python3
"""对每个节段运行 Nextclade，并合并成一张样本级 TSV。"""
from __future__ import annotations

import argparse
import csv
import logging
import re
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


KNOWN_SEGMENTS = {"PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"}


def split_by_segment(records: dict[str, str]) -> dict[str, list[tuple[str, str]]]:
    """从 header 中识别节段。

    支持 sample|HA、sample_HA、sample-HA、sample.HA、sample|HA|extra 等写法。
    识别不到已知节段时，保留最后一个 token，便于后续报错诊断。
    """
    out: dict[str, list[tuple[str, str]]] = {}
    for header, seq in records.items():
        tokens = [t for t in re.split(r"[|_\-\s.]+", header) if t]
        seg = next(
            (t for t in tokens if t.upper() in KNOWN_SEGMENTS),
            tokens[-1] if tokens else "",
        )
        seg = seg.upper() if seg.upper() in KNOWN_SEGMENTS else seg
        out.setdefault(seg, []).append((header, seq))
    return out


def sample_from_header(header: str) -> str:
    if "|" in header:
        return header.split("|", 1)[0]
    tokens = [t for t in re.split(r"[|_\-\s.]+", header) if t]
    return tokens[0] if tokens else header


def is_cli_mismatch(stderr: str) -> bool:
    """判断失败是否由 Nextclade 版本间 CLI 参数差异导致。"""
    text = stderr.lower()
    return any(
        key in text
        for key in (
            "is removed in favor of",
            "unexpected argument",
            "unrecognized",
            "positional arguments",
        )
    )


def run_nextclade(
    fasta_path: Path, dataset: str, outdir: Path, tsv_path: Path
) -> None:
    ensure_dir(outdir)
    # Nextclade >= 3.29：输入 FASTA 为位置参数，--output-dir 改为 --output-all
    modern = [
        "nextclade", "run",
        "--dataset-name", dataset,
        "--output-all", str(outdir),
        "--output-tsv", str(tsv_path),
        str(fasta_path),
    ]
    # Nextclade < 3.29：--input-fasta + --output-dir
    legacy = [
        "nextclade", "run",
        "--input-fasta", str(fasta_path),
        "--output-dir", str(outdir),
        "--dataset-name", dataset,
        "--output-tsv", str(tsv_path),
    ]
    # 更旧版本可能没有 --output-tsv，运行后从输出目录复制 nextclade.tsv
    legacy_no_tsv = [
        "nextclade", "run",
        "--input-fasta", str(fasta_path),
        "--output-dir", str(outdir),
        "--dataset-name", dataset,
    ]

    errors: list[str] = []
    proc = None
    for cmd in (modern, legacy, legacy_no_tsv):
        proc = subprocess.run(cmd, capture_output=True, text=True)
        if proc.returncode == 0:
            break
        errors.append(proc.stderr[-2000:])
        # 只有疑似参数不兼容时才继续尝试旧语法
        if not is_cli_mismatch(proc.stderr):
            break
    if proc is None or proc.returncode != 0:
        raise RuntimeError(errors[-1] if errors else "Nextclade 运行失败")
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
    datasets = {
        key.upper(): value
        for key, value in parse_pairs(args.dataset_map).items()
    }
    records = read_fasta(args.consensus)
    if not records:
        sys.exit(
            f"[ERROR] consensus FASTA 为空或没有序列：{args.consensus}\n"
            "  请检查组装是否成功（results/consensus/*.stats.tsv 和 IRMA 日志），"
            "然后重新运行 consensus。"
        )
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
        errors = ""
        try:
            run_nextclade(fasta, datasets[seg], seg_dir / "nextclade_out", tsv)
            row = parse_first_row(tsv)
            if not row:
                raise RuntimeError("Nextclade TSV 没有结果行")
            errors = get(row, "errors")
            if errors:
                status = "excluded"
                logging.warning(
                    "节段 %s 未匹配数据集 %s：%s",
                    seg, datasets[seg], errors[:200],
                )
        except Exception as exc:  # noqa: BLE001
            logging.error("Nextclade 失败：%s / %s: %s", seg, datasets[seg], exc)
            status = "failed"

        sample = sample_from_header(by_seg[seg][0][0])
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
                "errors": errors,
            }
        )

    if not out_rows:
        sys.exit(
            "[ERROR] 没有产生任何 Nextclade 结果。\n"
            f"  consensus 中的节段: {sorted(by_seg)}\n"
            f"  dataset-map 中的节段: {sorted(datasets)}\n"
            "  请确认 FASTA header 形如 sample|HA，且 dataset-map 的键"
            "（HA/NA/...）与节段名一致。"
        )
    write_tsv(out_rows, args.out)
    logging.info("写出 %d 个节段的 Nextclade 结果：%s", len(out_rows), args.out)


if __name__ == "__main__":
    main()
