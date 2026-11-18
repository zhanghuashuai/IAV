#!/usr/bin/env python3
"""基于 8 节段树拓扑相似性的重配筛查。

方法：把每棵树转换成“共享 taxa 的 clade 集合”，两两计算
Jaccard 距离 1 - |交集| / |并集|。距离越大，拓扑越不一致。
这是筛查，不是重配的最终证据；需结合支持度和元数据。
"""
from __future__ import annotations

import argparse
import csv
import itertools
import sys
from pathlib import Path

import dendropy

from utils import ensure_dir, parse_pairs, setup_logging, write_tsv

SEGMENT_ORDER = ["PB2", "PB1", "PA", "HA", "NP", "NA", "MP", "NS"]


def clade_set(tree: dendropy.Tree, common: set[str]) -> set[frozenset]:
    out: set[frozenset] = set()
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            continue
        labels = frozenset(
            taxon.label for taxon in node.leaf_iter()
            if taxon.label in common
        )
        if 1 < len(labels) < len(common):
            out.add(labels)
    return out


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree-map", required=True,
                   help="逗号分隔，如 PB2=...,PB1=...,...")
    p.add_argument("--outdir", required=True)
    p.add_argument("--threshold", type=float, default=0.30)
    args = p.parse_args()

    setup_logging()
    tree_map = parse_pairs(args.tree_map)
    if len(tree_map) < 2:
        sys.exit("[ERROR] 至少需要两棵节段树。")

    trees: dict[str, dendropy.Tree] = {}
    for seg, path in tree_map.items():
        if not Path(path).exists():
            sys.exit(f"[ERROR] 找不到树文件：{path}")
        trees[seg] = dendropy.Tree.get(
            path=path, schema="newick", preserve_underscores=True
        )

    taxa_sets = {
        seg: {t.label for t in tree.taxon_namespace}
        for seg, tree in trees.items()
    }
    common = set.intersection(*taxa_sets.values())
    if len(common) < 4:
        sys.exit(
            f"[ERROR] 8 个节段的共享 taxa 只有 {len(common)} 个，无法比较。"
        )

    segs = [s for s in SEGMENT_ORDER if s in trees]
    segs += [s for s in trees if s not in segs]
    clades = {seg: clade_set(trees[seg], common) for seg in segs}

    dist: dict[tuple[str, str], float] = {}
    for a, b in itertools.combinations(segs, 2):
        union = clades[a] | clades[b]
        inter = clades[a] & clades[b]
        d = 1.0 - (len(inter) / len(union)) if union else 0.0
        dist[(a, b)] = dist[(b, a)] = d

    outdir = ensure_dir(args.outdir)

    matrix_path = outdir / "rf_matrix.tsv"
    with open(matrix_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["segment"] + segs)
        for a in segs:
            writer.writerow(
                [a] + [f"{dist.get((a, b), 0.0):.4f}" for b in segs]
            )

    candidates = []
    for a, b in itertools.combinations(segs, 2):
        d = dist[(a, b)]
        if d >= args.threshold:
            candidates.append(
                {
                    "segment_1": a,
                    "segment_2": b,
                    "clade_jaccard_distance": f"{d:.4f}",
                    "shared_taxa": len(common),
                }
            )
    candidates.sort(
        key=lambda r: float(r["clade_jaccard_distance"]), reverse=True
    )
    write_tsv(candidates, outdir / "candidates.tsv")

    summary = [
        {
            "segment": seg,
            "n_taxa": len(taxa_sets[seg]),
            "n_shared_taxa": len(common),
            "n_clades": len(clades[seg]),
        }
        for seg in segs
    ]
    write_tsv(summary, outdir / "taxa_summary.tsv")

    print(f"matrix:     {matrix_path}")
    print(f"candidates: {outdir / 'candidates.tsv'}")


if __name__ == "__main__":
    main()
