"""Command-line entry point. Mirrors the upstream TRASH.R argument surface."""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

from . import __version__
from . import _log as log
from .pipeline import run_pipeline


def _pin_to_first_cpus(processes: int, quiet: bool) -> None:
    """Confine this process (and every nhmmer/clustalo child it spawns) to the
    first ``processes`` CPU cores, so total CPU usage tracks ``-p`` instead of
    the whole machine lighting up from many short-lived nhmmer calls.

    The affinity mask is inherited across fork and exec, so pinning here — before
    any worker pool is created — covers the entire process tree. No-op on
    platforms without an affinity API (macOS, Windows).
    """
    if not hasattr(os, "sched_setaffinity"):
        return
    available = sorted(os.sched_getaffinity(0))
    cpus = set(available[:processes])
    os.sched_setaffinity(0, cpus)
    if not quiet:
        print(f"pinned to {len(cpus)} core(s): {sorted(cpus)}", file=sys.stderr)


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
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    if args.processes < 1:
        print(f"--processes must be >= 1 (got {args.processes})", file=sys.stderr)
        return 2
    log.configure(quiet=args.quiet)
    _pin_to_first_cpus(args.processes, args.quiet)
    if not args.fasta.exists():
        print(f"fasta not found: {args.fasta}", file=sys.stderr)
        return 2
    args.output.mkdir(parents=True, exist_ok=True)
    run_pipeline(args)
    return 0
