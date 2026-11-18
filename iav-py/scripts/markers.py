#!/usr/bin/env python3
"""根据 Nextclade 的 aaSubstitutions / aaDeletions 扫描功能标记。

Nextclade 的突变字符串通常形如：
  HA:Q226L,HA:N158D
  M2:S31N
  NS1:80-84
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import yaml

from utils import setup_logging, write_tsv

SUB_RE = re.compile(
    r"^(?:(?P<gene>[^:]+):)?"
    r"(?P<ref>[A-Za-z*])(?P<pos>\d+)(?P<alt>[A-Za-z*])$"
)
DEL_RE = re.compile(
    r"^(?:(?P<gene>[^:]+):)?(?P<start>\d+)-(?P<end>\d+)"
)


def parse_substitutions(text: str | None) -> list[dict]:
    out: list[dict] = []
    if not text:
        return out
    for token in text.replace(" ", "").split(","):
        if not token:
            continue
        m = SUB_RE.match(token)
        if not m:
            continue
        out.append(
            {
                "gene": (m.group("gene") or "").strip(),
                "ref": m.group("ref").upper(),
                "pos": int(m.group("pos")),
                "alt": m.group("alt").upper(),
                "raw": token,
            }
        )
    return out


def parse_deletions(text: str | None) -> list[dict]:
    out: list[dict] = []
    if not text:
        return out
    for token in text.replace(" ", "").split(","):
        if not token:
            continue
        m = DEL_RE.match(token)
        if not m:
            continue
        out.append(
            {
                "gene": (m.group("gene") or "").strip(),
                "start": int(m.group("start")),
                "end": int(m.group("end")),
                "raw": token,
            }
        )
    return out


def gene_matches(marker_gene: str, observed_gene: str) -> bool:
    if not marker_gene or not observed_gene:
        return True
    return marker_gene.lower() == observed_gene.lower()


def is_applicable(marker: dict, dataset: str) -> bool:
    """如果 marker 指定了 dataset_contains，则只在匹配的数据集上评估。"""
    needle = str(marker.get("dataset_contains") or "").lower()
    if not needle or not dataset:
        return True
    return needle in dataset.lower()


def marker_found(
    marker: dict, subs: list[dict], dels: list[dict]
) -> tuple[bool, str]:
    marker_type = str(marker.get("type", "substitution")).lower()
    gene = str(marker.get("gene") or "")

    if marker_type == "substitution":
        for sub in subs:
            if sub["pos"] != int(marker["position"]):
                continue
            if sub["ref"] != str(marker["ref"]).upper():
                continue
            if sub["alt"] != str(marker["alt"]).upper():
                continue
            if not gene_matches(gene, sub["gene"]):
                continue
            return True, sub["raw"]
    elif marker_type == "deletion":
        start = int(marker["start"])
        end = int(marker["end"])
        for dele in dels:
            if dele["start"] <= start and dele["end"] >= end:
                if gene_matches(gene, dele["gene"]):
                    return True, dele["raw"]
    return False, ""


def load_markers(path: str) -> list[dict]:
    data = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if not isinstance(data, dict) or "markers" not in data:
        sys.exit(f"[ERROR] markers 文件格式不对：{path}")
    return data["markers"]


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--nextclade", required=True)
    p.add_argument("--markers", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    setup_logging()
    markers = load_markers(args.markers)

    with open(args.nextclade, encoding="utf-8-sig") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        sys.exit(f"[ERROR] Nextclade TSV 为空：{args.nextclade}")

    sample = rows[0].get("sample") or rows[0].get("seqName", "unknown").split("|")[0]
    rows_by_segment = {row.get("segment", ""): row for row in rows}

    out_rows: list[dict] = []
    for marker in markers:
        segment = marker.get("segment", "")
        row = rows_by_segment.get(segment, {})
        dataset = str(row.get("dataset", ""))
        applicable = segment in rows_by_segment and is_applicable(marker, dataset)
        if applicable:
            subs = parse_substitutions(
                row.get("aa_substitutions") or row.get("aaSubstitutions")
            )
            dels = parse_deletions(
                row.get("aa_deletions") or row.get("aaDeletions")
            )
            found, observed = marker_found(marker, subs, dels)
        else:
            found, observed = False, ""
        position = marker.get("position")
        if position is None and "start" in marker:
            position = f"{marker['start']}-{marker['end']}"
        out_rows.append(
            {
                "sample": sample,
                "segment": segment,
                "marker_id": marker.get("id", ""),
                "category": marker.get("category", ""),
                "label": marker.get("label", ""),
                "type": marker.get("type", "substitution"),
                "position": position,
                "ref": marker.get("ref", ""),
                "alt": marker.get("alt", ""),
                "found": "yes" if found else "no",
                "applicable": "yes" if applicable else "no",
                "observed": observed,
                "note": marker.get("note", ""),
            }
        )

    write_tsv(out_rows, args.out)


if __name__ == "__main__":
    main()
