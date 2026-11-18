#!/usr/bin/env python3
"""用 Biopython + matplotlib 快速绘制系统发育树。

出版级美化建议使用 ggtree、iTOL 或 toytree。
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from Bio import Phylo  # noqa: E402

from utils import ensure_dir, setup_logging  # noqa: E402


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--tree", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--title", default="")
    args = p.parse_args()

    setup_logging()
    tree = Phylo.read(args.tree, "newick")
    n_tips = len(tree.get_terminals())
    height = max(6.0, n_tips * 0.12)

    fig, ax = plt.subplots(figsize=(8, height))
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        label_func=lambda c: c.name if c.is_terminal() else "",
    )
    if args.title:
        ax.set_title(args.title)
    ax.set_axis_off()
    ensure_dir(Path(args.out).parent)
    fig.savefig(args.out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"figure: {args.out}")


if __name__ == "__main__":
    main()
