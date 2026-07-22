"""Command-line entry point. Mirrors the upstream TRASH.R argument surface.

Two entry paths:
* ``trash-py -f in.fasta -o out ...`` — the tandem-repeat pipeline (default).
* ``trash-py hor <repeats_with_seq.csv> ...`` — HOR detection on an existing
  repeat table (the port of the upstream ``HORT.R`` module).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import __version__
from . import _log as log
from .pipeline import output_stem, run_pipeline
from .hor_cli import (
    add_hor_arguments,
    build_hor_parser,
    hor_requested,
    run_hor_after_pipeline,
    run_hor_cli,
    warn_deprecated_aliases,
)


def fasta_path(value: str) -> Path:
    """Resolve the `-f` argument, mapping the conventional `-` to stdin.

    /dev/stdin keeps the rest of the pipeline working on a plain path, and
    its basename doubles as a sensible default output prefix (`stdin_*`).
    """
    return Path("/dev/stdin") if value == "-" else Path(value)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="trash-py",
        description="TRASH — tandem-repeat array identifier (Python)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="The --hor-* options above configure an optional HOR-detection stage "
               "that runs only AFTER array identification finishes.\n"
               "To run HOR detection on its own on an existing repeat table, use the "
               "`hor` subcommand:  trash-py hor --help",
    )
    p.add_argument(
        "-V", "--version", action="version", version=f"trash-py {__version__}"
    )
    p.add_argument("-f", "--fasta", required=True, type=fasta_path,
                   help="input fasta; `-` reads from stdin")
    p.add_argument("-o", "--output", required=True, type=Path, help="output directory")
    p.add_argument("-n", "--name", default=None,
                   help="prefix for output filenames (default: the input filename, "
                        "or `stdin` when reading from stdin)")
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

    # Register `hor` so it shows up under `trash-py --help`. Actual parsing of
    # `trash-py hor ...` is handled by the manual dispatch in main() (which owns
    # the full HORT.R-compatible parser); this stub exists only for discoverability.
    sub = p.add_subparsers(title="subcommands")
    sub.add_parser(
        "hor", add_help=False,
        help="detect higher-order repeats on an existing repeat table "
             "(standalone; see `trash-py hor --help`)")
    return p


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv

    # `trash-py hor ...` is a separate subcommand for HOR detection; everything
    # else is the original tandem-repeat pipeline, unchanged.
    if argv and argv[0] == "hor":
        hor_parser = build_hor_parser()
        ns = hor_parser.parse_args(argv[1:])
        log.configure(quiet=ns.quiet)
        warn_deprecated_aliases(argv[1:])
        return run_hor_cli(ns, hor_parser)

    args = build_parser().parse_args(argv)
    if args.processes < 1:
        print(f"--processes must be >= 1 (got {args.processes})", file=sys.stderr)
        return 2
    log.configure(quiet=args.quiet)
    if not args.fasta.exists():
        print(f"fasta not found: {args.fasta}", file=sys.stderr)
        return 2
    if args.name is not None and Path(args.name).name != args.name:
        print(f"--name must be a bare filename prefix, not a path: {args.name}",
              file=sys.stderr)
        return 2
    args.output.mkdir(parents=True, exist_ok=True)
    run_pipeline(args)

    # "Run and done": optionally detect HORs on the table we just produced.
    if hor_requested(args):
        repeats_with_seq = args.output / f"{output_stem(args)}_repeats_with_seq.csv"
        if not repeats_with_seq.exists():
            print(f"cannot run HOR: {repeats_with_seq} not found", file=sys.stderr)
            return 1
        return run_hor_after_pipeline(args, repeats_with_seq)
    return 0
