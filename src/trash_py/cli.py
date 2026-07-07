"""Command-line entry point. Mirrors the upstream TRASH.R argument surface."""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from . import _log as log
from .pipeline import run_pipeline
from .hor_cli import add_hor_arguments, hor_requested, run_hor_stage


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="trash-py", description="TRASH — tandem-repeat array identifier (Python)"
    )
    p.add_argument(
        "-V", "--version", action="version", version=f"trash-py {__version__}"
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
    p.add_argument("-q", "--quiet", action="store_true", help="suppress progress output")
    p.add_argument(
        "-p",
        "--processes",
        type=int,
        default=1,
        help="parallel worker processes for the array-identification and "
        "repeat-mapping stages (default 1 = serial)",
    )
    add_hor_arguments(p)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.processes < 1:
        print(f"--processes must be >= 1 (got {args.processes})", file=sys.stderr)
        return 2
    log.configure(quiet=args.quiet)
    if not args.fasta.exists():
        print(f"fasta not found: {args.fasta}", file=sys.stderr)
        return 2
    args.output.mkdir(parents=True, exist_ok=True)
    run_pipeline(args)

    # Optional HOR detection stage, on the repeat table we just wrote.
    if hor_requested(args):
        repeats_with_seq = args.output / f"{Path(args.fasta).name}_repeats_with_seq.csv"
        if not repeats_with_seq.exists():
            print(f"cannot run HOR: {repeats_with_seq} not found", file=sys.stderr)
            return 1
        return 1 if run_hor_stage(args, repeats_with_seq) else 0
    return 0
