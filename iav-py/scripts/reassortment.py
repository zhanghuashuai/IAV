#!/usr/bin/env python3
"""基于 8 节段树拓扑相似性的重配筛查。

方法：把每棵树转换成“共享 taxa 的 clade 集合”，两两计算
Jaccard 距离 1 - |交集| / |并集|。距离越大，拓扑越不一致。
这是筛查，不是重配的最终证据；需结合支持度和元数据。

数据不足时（树缺失、共享 taxa 少于阈值）写出 not_evaluable 状态并
正常退出，不阻塞主流程。
"""
from __future__ import annotations

import argparse
import csv
import itertools
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
            leaf.taxon.label
            for leaf in node.leaf_iter()
            if leaf.taxon is not None and leaf.taxon.label in common
        )
        if 1 < len(labels) < len(common):
            out.add(labels)
    return out


def ordered_segments(tree_map: dict[str, str]) -> list[str]:
    segs = [s for s in SEGMENT_ORDER if s in tree_map]
    return segs + [s for s in tree_map if s not in segs]


def write_matrix(
    path: Path, segs: list[str], dist: dict[tuple[str, str], float] | None = None
) -> None:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["segment"] + segs)
        if dist is None:
            return
        for a in segs:
            writer.writerow(
                [a] + [f"{dist.get((a, b), 0.0):.4f}" for b in segs]
            )


def write_status(
    outdir: Path, status: str, reason: str, n_trees: int, n_shared: int
) -> None:
    write_tsv(
        [
            {
                "status": status,
                "reason": reason,
                "n_trees": n_trees,
                "n_shared_taxa": n_shared,
            }
        ],
        outdir / "status.tsv",
    )


def write_taxa_summary(
    outdir: Path,
    tree_map: dict[str, str],
    taxa_sets: dict[str, set[str]],
    common: set[str],
    clades: dict[str, set[frozenset]] | None = None,
) -> None:
    rows = []
    for seg in ordered_segments(tree_map):
        path = tree_map[seg]
        rows.append(
            {
                "segment": seg,
                "tree": path,
                "exists": "yes" if Path(path).exists() else "no",
                "n_taxa": len(taxa_sets.get(seg, set())),
                "n_shared_taxa": len(common),
                "n_clades": len(clades.get(seg, set())) if clades else 0,
            }
        )
    write_tsv(rows, outdir / "taxa_summary.tsv")


def not_evaluable(
    outdir: Path,
    reason: str,
    tree_map: dict[str, str],
    taxa_sets: dict[str, set[str]],
    common: set[str],
    clades: dict[str, set[frozenset]] | None = None,
) -> None:
    segs = ordered_segments(tree_map)
    write_matrix(outdir / "rf_matrix.tsv", segs)
    write_tsv(
        [
            {
                "status": "not_evaluable",
                "reason": reason,
                "segment_1": "",
                "segment_2": "",
                "clade_jaccard_distance": "",
                "shared_taxa": len(common),
            }
        ],
        outdir / "candidates.tsv",
    )
    write_taxa_summary(outdir, tree_map, taxa_sets, common, clades)
    write_status(outdir, "not_evaluable", reason, len(tree_map), len(common))
    print(f"[WARN] reassortment not evaluable: {reason}")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree-map", required=True,
                   help="逗号分隔，如 PB2=...,PB1=...,...")
    p.add_argument("--outdir", required=True)
    p.add_argument("--threshold", type=float, default=0.30)
    p.add_argument("--min-shared-taxa", type=int, default=4)
    args = p.parse_args()

    setup_logging()
    tree_map = parse_pairs(args.tree_map)
    outdir = ensure_dir(args.outdir)

    if len(tree_map) < 2:
        not_evaluable(
            outdir,
            "至少需要两棵节段树。",
            tree_map,
            {},
            set(),
        )
        return

    missing = [
        seg for seg, path in tree_map.items() if not Path(path).exists()
    ]
    if missing:
        # 已存在的树仍统计 taxa 数，便于查看还缺什么
        taxa_sets: dict[str, set[str]] = {}
        for seg, path in tree_map.items():
            if seg in missing:
                continue
            try:
                tree = dendropy.Tree.get(
                    path=path, schema="newick", preserve_underscores=True
                )
                taxa_sets[seg] = {
                    t.label for t in tree.taxon_namespace
                }
            except Exception:  # noqa: BLE001
                continue
        not_evaluable(
            outdir,
            f"缺少节段树：{', '.join(sorted(missing))}。",
            tree_map,
            taxa_sets,
            set(),
        )
        return

    trees: dict[str, dendropy.Tree] = {}
    for seg, path in tree_map.items():
        trees[seg] = dendropy.Tree.get(
            path=path, schema="newick", preserve_underscores=True
        )

    taxa_sets = {
        seg: {t.label for t in tree.taxon_namespace}
        for seg, tree in trees.items()
    }
    common = set.intersection(*taxa_sets.values()) if taxa_sets else set()

    if len(common) < args.min_shared_taxa:
        not_evaluable(
            outdir,
            f"8 个节段树的共享 taxa 只有 {len(common)} 个，"
            f"少于 {args.min_shared_taxa} 个；通常需要至少 "
            f"{args.min_shared_taxa} 个共同样本才能做重配筛查。",
            tree_map,
            taxa_sets,
            common,
        )
        return

    segs = ordered_segments(tree_map)
    clades = {seg: clade_set(trees[seg], common) for seg in segs}

    dist: dict[tuple[str, str], float] = {}
    for a, b in itertools.combinations(segs, 2):
        union = clades[a] | clades[b]
        inter = clades[a] & clades[b]
        d = 1.0 - (len(inter) / len(union)) if union else 0.0
        dist[(a, b)] = dist[(b, a)] = d

    matrix_path = outdir / "rf_matrix.tsv"
    write_matrix(matrix_path, segs, dist)

    candidates = []
    for a, b in itertools.combinations(segs, 2):
        d = dist[(a, b)]
        if d >= args.threshold:
            candidates.append(
                {
                    "status": "ok",
                    "reason": "",
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

    write_taxa_summary(outdir, tree_map, taxa_sets, common, clades)
    write_status(outdir, "ok", "", len(tree_map), len(common))

    print(f"matrix:     {matrix_path}")
    print(f"candidates: {outdir / 'candidates.tsv'}")


if __name__ == "__main__":
    main()
