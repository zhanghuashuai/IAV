#!/usr/bin/env python3
"""用 TreeTime 构建时间标定树，并整理分子钟摘要。"""
from __future__ import annotations

import argparse
import csv
import shlex
import sys
from pathlib import Path

import dendropy

from utils import ensure_dir, require_tool, run, setup_logging, write_tsv


def build_dates(
    metadata: Path, name_col: str, date_col: str, out: Path
) -> None:
    with open(metadata, encoding="utf-8-sig", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"[ERROR] 元数据为空：{metadata}")
    missing = {name_col, date_col} - set(rows[0].keys())
    if missing:
        sys.exit(
            f"[ERROR] 元数据缺少列：{', '.join(sorted(missing))}；"
            f"现有列：{', '.join(rows[0].keys())}"
        )
    with open(out, "w", encoding="utf-8") as fh:
        fh.write("name\tdate\n")
        for row in rows:
            name = (row.get(name_col) or "").strip()
            date = (row.get(date_col) or "").strip()
            if name and date:
                fh.write(f"{name}\t{date}\n")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree", required=True)
    p.add_argument("--aln", required=True)
    p.add_argument("--metadata", required=True)
    p.add_argument("--outdir", default="results/phylo")
    p.add_argument("--prefix", required=True)
    p.add_argument("--out-tree", required=True)
    p.add_argument("--out-clock", required=True)
    p.add_argument("--name-col", default="name")
    p.add_argument("--date-col", default="date")
    p.add_argument("--treetime-args", default="--reroot least-squares")
    args = p.parse_args()

    setup_logging()
    treetime = require_tool("treetime")
    outdir = ensure_dir(args.outdir)

    dates = outdir / "dates.tsv"
    build_dates(Path(args.metadata), args.name_col, args.date_col, dates)

    cmd = [
        treetime,
        "--tree", args.tree,
        "--aln", args.aln,
        "--dates", str(dates),
        "--outdir", str(outdir),
    ] + shlex.split(args.treetime_args)
    run(cmd, log_path=outdir / "treetime.log")

    nexus = outdir / "timetree.nexus"
    if not nexus.exists():
        sys.exit(f"[ERROR] TreeTime 未生成 timetree.nexus：{outdir}")
    tree = dendropy.Tree.get(
        path=str(nexus), schema="nexus", preserve_underscores=True
    )
    ensure_dir(Path(args.out_tree).parent)
    tree.write(path=args.out_tree, schema="newick")

    clock = outdir / "molecular_clock.txt"
    rows: list[dict] = []
    if clock.exists():
        for line in clock.read_text(encoding="utf-8").splitlines():
            if ":" in line:
                key, value = line.split(":", 1)
                rows.append({"key": key.strip(), "value": value.strip()})
    if not rows:
        rows = [{"key": "status", "value": "molecular_clock.txt not found"}]
    write_tsv(rows, args.out_clock, fieldnames=["key", "value"])

    print(f"timetree: {args.out_tree}")
    print(f"clock:    {args.out_clock}")


if __name__ == "__main__":
    main()
