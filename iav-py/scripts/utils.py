"""通用工具函数：日志、外部命令、FASTA/TSV 读写。"""
from __future__ import annotations

import csv
import gzip
import logging
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Iterable


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def require_tool(name: str) -> str:
    path = shutil.which(name)
    if not path:
        sys.exit(
            f"[ERROR] 找不到可执行文件：{name}\n"
            f"请先激活 conda 环境（envs/flubio.yml 或 envs/flurnaseq.yml），"
            f"或检查 PATH。"
        )
    return path


def run(
    cmd: Iterable[str],
    log_path: str | Path | None = None,
    cwd: str | Path | None = None,
) -> None:
    """运行一个外部命令；可将 stdout/stderr 追加到日志。"""
    cmd = [str(x) for x in cmd]
    logging.info("RUN: %s", " ".join(shlex.quote(x) for x in cmd))
    if log_path:
        log_path = Path(log_path)
        ensure_dir(log_path.parent)
        with open(log_path, "a", encoding="utf-8") as fh:
            fh.write("$ " + " ".join(shlex.quote(x) for x in cmd) + "\n")
            fh.flush()
            subprocess.run(
                cmd, check=True, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT
            )
    else:
        subprocess.run(cmd, check=True, cwd=cwd)


def run_pipeline(
    cmds: list[list[str]], log_path: str | Path | None = None
) -> None:
    """把多个命令用管道连接，例如 minimap2 | samtools sort。"""
    cmds = [[str(x) for x in c] for c in cmds]
    logging.info(
        "PIPE: %s",
        " | ".join(" ".join(shlex.quote(x) for x in c) for c in cmds),
    )
    if log_path:
        log_path = Path(log_path)
        ensure_dir(log_path.parent)
    log_fh = (
        open(log_path, "a", encoding="utf-8")
        if log_path
        else subprocess.DEVNULL
    )
    try:
        procs: list[subprocess.Popen] = []
        prev_stdout = None
        for i, cmd in enumerate(cmds):
            stdout = subprocess.PIPE if i < len(cmds) - 1 else log_fh
            p = subprocess.Popen(
                cmd, stdin=prev_stdout, stdout=stdout, stderr=log_fh
            )
            if prev_stdout is not None:
                prev_stdout.close()
            prev_stdout = p.stdout
            procs.append(p)
        for p in procs:
            if p.wait() != 0:
                raise subprocess.CalledProcessError(p.returncode, p.args)
    finally:
        if log_path:
            log_fh.close()


def parse_pairs(text: str | None) -> dict[str, str]:
    """解析 "HA=flu_h1n1pdm_ha,NA=..." 形式的参数。"""
    out: dict[str, str] = {}
    if not text:
        return out
    for item in text.split(","):
        item = item.strip()
        if not item:
            continue
        if "=" not in item:
            raise ValueError(f"无法解析键值对：{item}")
        key, value = item.split("=", 1)
        out[key.strip()] = value.strip()
    return out


def count_fastq_reads(path: str | Path) -> int:
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    n = 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for _ in fh:
            n += 1
    return n // 4


def read_fasta(path: str | Path) -> dict[str, str]:
    """极简 FASTA 读取，返回 {header: sequence}（header 不含 '>'）。"""
    seqs: dict[str, list[str]] = {}
    name: str | None = None
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                name = line[1:].strip()
                seqs[name] = []
            elif name is not None:
                seqs[name].append(line.strip())
    return {k: "".join(v).upper() for k, v in seqs.items()}


def write_fasta(
    records: dict[str, str], path: str | Path, width: int = 60
) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records.items():
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i : i + width] + "\n")


def write_tsv(
    rows: list[dict],
    path: str | Path,
    fieldnames: list[str] | None = None,
) -> None:
    path = Path(path)
    ensure_dir(path.parent)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def setup_logging(level: int = logging.INFO) -> None:
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )
