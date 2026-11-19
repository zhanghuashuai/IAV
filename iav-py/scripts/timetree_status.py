#!/usr/bin/env python3
"""记录 TreeTime 是否可以运行（metadata 是否存在且非空）。"""
from __future__ import annotations

import argparse
from pathlib import Path

from utils import setup_logging, write_tsv


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--metadata", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    setup_logging()
    metadata = Path(args.metadata)
    ok = metadata.is_file() and metadata.stat().st_size > 0
    write_tsv(
        [
            {
                "status": "ok" if ok else "skipped",
                "reason": (
                    ""
                    if ok
                    else f"metadata 不存在或为空：{metadata}；"
                    "TreeTime 需要 name/date 两列。"
                ),
                "metadata": str(metadata),
            }
        ],
        args.out,
    )


if __name__ == "__main__":
    main()
