"""HOR detection CLI — usable two ways, sharing one core runner:

* As a stage of the main pipeline: ``trash-py -f g.fasta -o out --hor-chr-list …``
  runs the pipeline and then detects HORs on the repeat table it just wrote.
  The "run and done" path.
* As a standalone subcommand on an existing table:
  ``trash-py hor g.fasta_repeats_with_seq.csv --chr-list … -o out``. Handy to
  re-run HOR finding for a different class or set of sequences without redoing
  the whole pipeline.

Mirrors the upstream ``HORT.R`` surface (``--ChrA``/``--ChrB``,
``--class``/``--classB``, ``-t``/``-l``) and adds two conveniences: a chromosome
list for the common self-comparison case, and an optional ``--class`` that
auto-selects the most abundant class on the target sequence(s) — useful on a new
organism where the major satellite isn't known yet.
"""
from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

from . import _log as log
from .hor import (
    HorArgs,
    available_seqids,
    class_abundance,
    run_hor_single,
    run_hor_pair,
)


def _truncate(ids: list[str], limit: int = 5) -> str:
    if len(ids) <= limit:
        return ", ".join(ids)
    return ", ".join(ids[:limit]) + f", … (+{len(ids) - limit} more)"


# ---------------------------------------------------------------------------
# Main-pipeline integration: HOR options on the top-level parser
# ---------------------------------------------------------------------------

def add_hor_arguments(p: argparse.ArgumentParser) -> None:
    """Add the (all ``--hor-``-prefixed) HOR options to the main pipeline parser.
    Passing ``--hor-chr-list`` or ``--hor-ChrA`` turns HOR detection on."""
    g = p.add_argument_group("HOR detection (optional; runs after the main pipeline)")
    g.add_argument("--hor-chr-list", dest="chr_list",
                   help="comma-separated sequence IDs — detect HORs (self-comparison) "
                        "for each. The common way to turn HOR detection on.")
    g.add_argument("--hor-ChrA", dest="chrA",
                   help="single sequence ID for HOR detection (region A if --hor-ChrB given)")
    g.add_argument("--hor-ChrB", dest="chrB",
                   help="region B sequence ID — enables cross-region comparison")
    g.add_argument("--hor-class", dest="hor_class",
                   help="repeat class to analyse (default: most abundant on the targets)")
    g.add_argument("--hor-classB", dest="classB",
                   help="region-B repeat class (default: same as --hor-class)")
    g.add_argument("--hor-threshold", dest="hor_threshold", type=int, default=25,
                   help="HOR divergence threshold %% (default 25)")
    g.add_argument("--hor-min-len", dest="hor_min_len", type=int, default=3,
                   help="minimum HOR length in repeat units (default 3)")
    g.add_argument("--hor-genomeA", dest="genomeA", default="A", help="label for genome A")
    g.add_argument("--hor-genomeB", dest="genomeB", default="B", help="label for genome B")
    g.add_argument("--no-hor-plot", dest="no_plot", action="store_true",
                   help="skip the HOR line plot")
    g.add_argument("--no-hor-saveR", dest="no_saveR", action="store_true",
                   help="--hor-ChrB mode only: skip writing repeats_with_hors")


def hor_requested(args: argparse.Namespace) -> bool:
    return bool(getattr(args, "chr_list", None) or getattr(args, "chrA", None))


def run_hor_after_pipeline(args: argparse.Namespace, repeats_path: Path) -> int:
    """Entry point for the main-pipeline path. `repeats_path` is the table the
    pipeline just wrote. Returns a process exit code."""
    if shutil.which("mafft") is None:
        log.warn("mafft not found on PATH — skipping HOR detection "
                 "(install e.g. `conda install -c bioconda mafft`)")
        return 1
    return _run(
        repeats=repeats_path,
        repeatsB=repeats_path,
        output=Path(args.output),
        chr_list=args.chr_list,
        chrA=args.chrA,
        chrB=args.chrB,
        hor_class=args.hor_class,
        classB=args.classB,
        hor_threshold=args.hor_threshold,
        hor_min_len=args.hor_min_len,
        genomeA=args.genomeA,
        genomeB=args.genomeB,
        make_plot=not args.no_plot,
        saveR=not args.no_saveR,
        threads=getattr(args, "processes", 1) or 1,
    )


# ---------------------------------------------------------------------------
# Standalone subcommand: `trash-py hor <repeats_with_seq.csv> ...`
# ---------------------------------------------------------------------------

def build_hor_parser() -> argparse.ArgumentParser:
    """The standalone HOR parser. Accepts the full upstream ``HORT.R`` getopt
    surface (short flags and the underscore long-names) so that replacing
    ``Rscript HORT.R`` with ``trash-py hor`` in an existing script Just Works,
    while also accepting a positional repeats file and a ``--chr-list``."""
    p = argparse.ArgumentParser(
        prog="trash-py hor",
        description="TRASH HOR module — detect higher-order repeats in a repeat table.",
    )
    # The repeats table: positional *or* -r/--repeats (HORT.R used -r).
    p.add_argument("repeats_pos", nargs="?", type=Path, metavar="repeats",
                   help="repeat table with a sequence column (…_repeats_with_seq.csv); "
                        "may also be given with -r/--repeats")
    p.add_argument("-r", "--repeats", dest="repeats", type=Path, default=None,
                   help="repeat table (HORT.R-compatible alternative to the positional arg)")
    p.add_argument("-o", "--output", "--output_folder", dest="output",
                   type=Path, default=Path.cwd(),
                   help="output folder (default: current directory)")
    p.add_argument("-c", "--class", dest="hor_class", default=None,
                   help="repeat class to analyse (e.g. 178_1); "
                        "default: the most abundant class on the target sequences")
    p.add_argument("-t", "--hor-threshold", "--hor_threshold", dest="hor_threshold",
                   type=int, default=25, help="divergence threshold %% (default 25)")
    p.add_argument("-l", "--hor-min-len", "--hor_min_len", dest="hor_min_len",
                   type=int, default=3, help="minimum HOR length in repeat units (default 3)")
    p.add_argument("-m", "--method", type=int, choices=(1, 2), default=None,
                   help="1 = self-comparison, 2 = cross-region (usually inferred from --ChrB)")
    p.add_argument("--chr-list", "--hor-chr-list", dest="chr_list",
                   help="comma-separated sequence IDs — self-comparison (method 1) for each")
    p.add_argument("-A", "--ChrA", "--chrA", dest="chrA",
                   help="single sequence ID for method 1 (or region A of method 2)")
    p.add_argument("-B", "--ChrB", "--chrB", dest="chrB",
                   help="region B sequence ID — enables cross-region comparison (method 2)")
    p.add_argument("-b", "--repeatsB", dest="repeatsB", type=Path,
                   help="region-B repeats table (default: same as the main repeats file)")
    p.add_argument("-C", "--classB", dest="classB",
                   help="region-B repeat class (default: same as --class)")
    p.add_argument("-g", "--genomeA", dest="genomeA", default="A",
                   help="label for genome A (--ChrB mode)")
    p.add_argument("-G", "--genomeB", dest="genomeB", default="B",
                   help="label for genome B (--ChrB mode)")
    p.add_argument("-s", "--saveR", dest="saveR", choices=("y", "n"), default=None,
                   help="--ChrB mode: write repeats_with_hors (y/n)")
    p.add_argument("-p", "--plot-simple", "--plot_simple", dest="plot_simple",
                   choices=("y", "n"), default=None,
                   help="accepted for HORT.R compatibility (a plot is always produced)")
    p.add_argument("--no-plot", dest="no_plot", action="store_true",
                   help="skip the HOR plot")
    p.add_argument("--no-saveR", dest="no_saveR", action="store_true",
                   help="--ChrB mode only: skip writing repeats_with_hors")
    p.add_argument("-T", "--threads", type=int, default=None,
                   help="MAFFT threads (default: all available cores)")
    p.add_argument("-q", "--quiet", action="store_true", help="suppress progress output")
    return p


def run_hor_cli(ns: argparse.Namespace) -> int:
    """Entry point for the standalone `trash-py hor` subcommand."""
    repeats = ns.repeats or ns.repeats_pos
    if repeats is None:
        print("no repeats table given (pass it positionally or with -r/--repeats)",
              file=sys.stderr)
        return 2
    if not repeats.exists():
        print(f"repeats file not found: {repeats}", file=sys.stderr)
        return 2
    if shutil.which("mafft") is None:
        print("mafft not found on PATH — install it (e.g. `conda install -c bioconda mafft`)",
              file=sys.stderr)
        return 2

    # Reconcile the explicit --method with the presence of --ChrB.
    if ns.method == 2 and not ns.chrB:
        print("method 2 (-m 2) requires a region B (--ChrB/-B)", file=sys.stderr)
        return 2
    if ns.method == 1 and ns.chrB:
        print("method 1 (-m 1) is a self-comparison; drop --ChrB or use -m 2", file=sys.stderr)
        return 2

    # saveR: -s y/n wins, else --no-saveR, else default on.
    if ns.saveR is not None:
        save = ns.saveR == "y"
    else:
        save = not ns.no_saveR

    ns.output.mkdir(parents=True, exist_ok=True)
    return _run(
        repeats=repeats,
        repeatsB=ns.repeatsB or repeats,
        output=ns.output,
        chr_list=ns.chr_list,
        chrA=ns.chrA,
        chrB=ns.chrB,
        hor_class=ns.hor_class,
        classB=ns.classB,
        hor_threshold=ns.hor_threshold,
        hor_min_len=ns.hor_min_len,
        genomeA=ns.genomeA,
        genomeB=ns.genomeB,
        make_plot=not ns.no_plot,
        saveR=save,
        threads=ns.threads if ns.threads is not None else (os.cpu_count() or 1),
    )


# ---------------------------------------------------------------------------
# Shared core
# ---------------------------------------------------------------------------

def _resolve_class(repeats: Path, seq_ids: set[str], explicit: str | None) -> str | None:
    """Return the class to analyse: the user's choice, or the most abundant class
    on the target sequences (logged so the auto-selection is never silent)."""
    if explicit:
        return explicit
    ranked = class_abundance(repeats, seq_ids)
    if not ranked:
        return None
    top_class, top_n = ranked[0]
    runners = ", ".join(f"{c} ({n})" for c, n in ranked[1:4])
    log.info(f"auto-selected class '{top_class}' ({top_n} repeats on "
             f"{len(seq_ids)} sequence{'s' if len(seq_ids) != 1 else ''})"
             + (f"; next: {runners}" if runners else ""))
    return top_class


def _run(*, repeats: Path, repeatsB: Path, output: Path, chr_list: str | None,
         chrA: str | None, chrB: str | None, hor_class: str | None, classB: str | None,
         hor_threshold: int, hor_min_len: int, genomeA: str, genomeB: str,
         make_plot: bool, saveR: bool, threads: int = 1) -> int:
    if not chr_list and not chrA:
        print("specify --hor-chr-list, or --hor-ChrA (optionally with --hor-ChrB)",
              file=sys.stderr)
        return 2

    # cross-region comparison (method 2)
    if chrB:
        if not chrA:
            print("--ChrB requires --ChrA", file=sys.stderr)
            return 2
        missing = ([chrA] if chrA not in available_seqids(repeats) else []) \
            + ([chrB] if chrB not in available_seqids(repeatsB) else [])
        if missing:
            print(f"sequence ID(s) not in the repeats table: {_truncate(missing)}",
                  file=sys.stderr)
            return 2
        resolved = _resolve_class(repeats, {chrA}, hor_class)
        if not resolved:
            print("could not determine a repeat class (empty table?)", file=sys.stderr)
            return 2
        args = HorArgs(repeats=repeats, output_folder=output, hor_class=resolved,
                       hor_threshold=hor_threshold, hor_min_len=hor_min_len,
                       make_plot=make_plot, threads=threads)
        try:
            run_hor_pair(args, chrA, chrB, classB=classB or resolved, repeatsB=repeatsB,
                         genomeA=genomeA, genomeB=genomeB, saveR=saveR)
        except Exception as e:
            log.warn(str(e))
            return 1
        return 0

    # one or more self-comparisons (method 1)
    chrs: list[str] = []
    if chr_list:
        chrs.extend(c.strip() for c in chr_list.split(",") if c.strip())
    if chrA:
        chrs.append(chrA)
    seen: set[str] = set()
    chrs = [c for c in chrs if not (c in seen or seen.add(c))]

    # Check requested IDs against what's actually in the table, up front, and
    # hard-fail (running nothing) if ANY are missing — so the user can't
    # silently miss sequences. The list is truncated to stay readable.
    missing = [c for c in chrs if c not in available_seqids(repeats)]
    if missing:
        print(f"not running HOR identification for {len(missing)} ID(s) "
              f"not in the repeats table: {_truncate(missing)}", file=sys.stderr)
        return 2

    resolved = _resolve_class(repeats, set(chrs), hor_class)
    if not resolved:
        print("could not determine a repeat class (empty table?)", file=sys.stderr)
        return 2
    args = HorArgs(repeats=repeats, output_folder=output, hor_class=resolved,
                   hor_threshold=hor_threshold, hor_min_len=hor_min_len, make_plot=make_plot, threads=threads)

    failures = 0
    for chrom in chrs:
        try:
            run_hor_single(args, chrom)
        except Exception as e:
            log.warn(f"{chrom}: {e}")
            failures += 1
    return 1 if failures else 0
