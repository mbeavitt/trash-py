"""Command-line entry point. Mirrors the upstream TRASH.R argument surface."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from .pipeline import run_pipeline


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="trash-py", description="TRASH — tandem-repeat array identifier (Python)"
    )
    p.add_argument("-f", "--fasta", required=True, type=Path, help="input fasta")
    p.add_argument("-o", "--output", required=True, type=Path, help="output directory")
    p.add_argument("-m", "--max-rep-size", type=int, default=1000)
    p.add_argument("-i", "--min-rep-size", type=int, default=7)
    p.add_argument(
        "-t",
        "--templates",
        type=Path,
        default=None,
        help="optional template fasta — assigns class names from headers",
    )
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.fasta.exists():
        print(f"fasta not found: {args.fasta}", file=sys.stderr)
        return 2
    args.output.mkdir(parents=True, exist_ok=True)
    run_pipeline(args)
    return 0
