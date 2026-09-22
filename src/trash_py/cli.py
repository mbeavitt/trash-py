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
from . import checkpoint as ckpt_mod
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

    ck = p.add_argument_group(
        "checkpointing",
        "Survive an HPC wall-clock limit: save progress as the run goes, stop "
        "cleanly when the scheduler warns us, and resume where it stopped when "
        "the same command runs again.",
    )
    ck.add_argument(
        "--checkpoint",
        nargs="?",
        const="",
        default=None,
        metavar="DIR",
        help="record resumable progress in DIR and resume from it if it "
             "already holds state for this run (default DIR: "
             "<output>/<name>.checkpoint)",
    )
    ck.add_argument(
        "--restart", action="store_true",
        help="discard any existing checkpoint and start from scratch",
    )
    ck.add_argument(
        "--keep-checkpoint", action="store_true",
        help="keep the checkpoint directory after a successful run "
             "(it is deleted by default)",
    )
    ck.add_argument(
        "--checkpoint-interval", type=float,
        default=ckpt_mod.DEFAULT_INTERVAL, metavar="SECONDS",
        help=f"how often to commit progress to disk "
             f"(default {ckpt_mod.DEFAULT_INTERVAL:.0f}s)",
    )
    ck.add_argument(
        "--time-limit", default=None, metavar="DURATION",
        help="wall-clock budget for this run, e.g. 11h, 690m, 11:30:00, or a "
             "bare number of seconds; when unset, SLURM_JOB_END_TIME is used "
             "if the scheduler exported it",
    )
    ck.add_argument(
        "--time-margin", type=str, default=None, metavar="DURATION",
        help=f"stop this long before the deadline, to leave room for saving "
             f"state and for in-flight tasks to finish (default "
             f"{ckpt_mod.DEFAULT_MARGIN / 60:.0f}m)",
    )
    ck.add_argument(
        "--checkpoint-signal", default=ckpt_mod.DEFAULT_SIGNALS, metavar="LIST",
        help=f"comma-separated signals that request a checkpoint-and-stop "
             f"(default {ckpt_mod.DEFAULT_SIGNALS})",
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

    try:
        checkpoint = _build_checkpoint(args)
    except ValueError as exc:
        print(f"{exc}", file=sys.stderr)
        return 2

    try:
        run_pipeline(args, checkpoint)

        # "Run and done": optionally detect HORs on the table we just produced.
        if hor_requested(args):
            repeats_with_seq = args.output / f"{output_stem(args)}_repeats_with_seq.csv"
            if not repeats_with_seq.exists():
                print(f"cannot run HOR: {repeats_with_seq} not found", file=sys.stderr)
                return 1
            status = run_hor_after_pipeline(args, repeats_with_seq)
            if status != 0:
                return status
    except ckpt_mod.Checkpointed as stop:
        # State is already on disk; re-running the same command resumes.
        print(
            f"trash-py stopped early ({stop.reason}) and saved its state to "
            f"{stop.directory}.\nRe-run the same command to resume.",
            file=sys.stderr,
        )
        return ckpt_mod.EXIT_CHECKPOINTED

    checkpoint.complete()
    return 0


def _build_checkpoint(args: argparse.Namespace):
    """Turn the --checkpoint/--time-* flags into a Checkpoint (or the no-op)."""
    if args.checkpoint is None:
        return ckpt_mod.NO_CHECKPOINT

    if not args.fasta.is_file():
        raise ValueError(
            f"--checkpoint needs a regular fasta file to resume from, "
            f"not {args.fasta}"
        )

    if args.checkpoint == "":
        args.checkpoint = args.output / f"{output_stem(args)}.checkpoint"

    margin = ckpt_mod.DEFAULT_MARGIN
    if args.time_margin is not None:
        margin = ckpt_mod.parse_duration(args.time_margin)
    args.time_margin = margin
    if args.time_limit is not None:
        ckpt_mod.parse_duration(args.time_limit)  # fail fast on a typo
    ckpt_mod.parse_signals(args.checkpoint_signal)

    return ckpt_mod.from_args(args, output_stem(args))
