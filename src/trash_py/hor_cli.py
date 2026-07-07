"""HOR options for the main CLI, plus the post-pipeline HOR stage.

The HOR module runs as the tail end of the normal pipeline: once the repeat
table (``…_repeats_with_seq.csv``) has been written, if the user asked for HOR
detection (``--chr-list`` for the common self-comparison case, or ``--ChrA``,
optionally with ``--ChrB`` for a cross-region comparison) we detect HORs on that
table. Mirrors the upstream ``HORT.R`` argument surface.
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from . import _log as log
from .hor import HorArgs, run_hor_single, run_hor_pair


def add_hor_arguments(p: argparse.ArgumentParser) -> None:
    g = p.add_argument_group("HOR detection (optional; runs after the main pipeline)")
    g.add_argument("--chr-list",
                   help="comma-separated sequence IDs — detect HORs (self-comparison) "
                        "for each. The common way to turn HOR detection on.")
    g.add_argument("-A", "--ChrA", dest="chrA",
                   help="single sequence ID for HOR detection (region A if --ChrB given)")
    g.add_argument("-B", "--ChrB", dest="chrB",
                   help="region B sequence ID — enables cross-region comparison")
    g.add_argument("-c", "--class", dest="hor_class",
                   help="repeat class to analyse for HORs (e.g. 178_1); required for HOR")
    g.add_argument("--hor-threshold", type=int, default=25,
                   help="HOR divergence threshold %% (default 25)")
    g.add_argument("--hor-min-len", type=int, default=3,
                   help="minimum HOR length in repeat units (default 3)")
    g.add_argument("--repeatsB", type=Path,
                   help="region-B repeats table (default: the pipeline's own output)")
    g.add_argument("--classB", help="region-B repeat class (default: same as --class)")
    g.add_argument("--genomeA", default="A", help="label for genome A (--ChrB mode)")
    g.add_argument("--genomeB", default="B", help="label for genome B (--ChrB mode)")
    g.add_argument("--no-hor-plot", action="store_true", help="skip the HOR line plot")
    g.add_argument("--no-saveR", action="store_true",
                   help="--ChrB mode only: skip writing repeats_with_hors")


def hor_requested(args: argparse.Namespace) -> bool:
    return bool(getattr(args, "chr_list", None) or getattr(args, "chrA", None))


def run_hor_stage(args: argparse.Namespace, repeats_path: Path) -> int:
    """Run the HOR module on `repeats_path`. Returns the number of failures."""
    if not args.hor_class:
        log.warn("HOR requested but no --class given; skipping HOR detection")
        return 1
    if shutil.which("mafft") is None:
        log.warn("mafft not found on PATH — skipping HOR detection "
                 "(install e.g. `conda install -c bioconda mafft`)")
        return 1

    hor_args = HorArgs(
        repeats=repeats_path,
        output_folder=Path(args.output),
        hor_class=args.hor_class,
        hor_threshold=args.hor_threshold,
        hor_min_len=args.hor_min_len,
        make_plot=not args.no_hor_plot,
    )

    failures = 0

    # cross-region comparison (method 2)
    if args.chrB:
        if not args.chrA:
            log.warn("--ChrB requires --ChrA; skipping HOR detection")
            return 1
        try:
            run_hor_pair(
                hor_args, args.chrA, args.chrB,
                classB=args.classB or args.hor_class,
                repeatsB=args.repeatsB or repeats_path,
                genomeA=args.genomeA, genomeB=args.genomeB,
                saveR=not args.no_saveR,
            )
        except Exception as e:
            log.warn(str(e))
            failures += 1
        return failures

    # one or more self-comparisons (method 1)
    chrs: list[str] = []
    if args.chr_list:
        chrs.extend(c.strip() for c in args.chr_list.split(",") if c.strip())
    if args.chrA:
        chrs.append(args.chrA)
    seen: set[str] = set()
    for chrom in chrs:
        if chrom in seen:
            continue
        seen.add(chrom)
        try:
            run_hor_single(hor_args, chrom)
        except Exception as e:
            log.warn(f"{chrom}: {e}")
            failures += 1
    return failures
