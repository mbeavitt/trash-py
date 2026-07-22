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


# Redundant spellings that trash-py itself introduced (a canonical form is
# preferred). We keep parsing them but nudge users toward the canonical name.
# The genuine upstream HORT.R getopt surface (-r, --output_folder, --ChrA,
# -m/--method, -s/--saveR, --plot_simple) is deliberately NOT listed here: it is
# the drag-and-drop contract and stays silent.
_DEPRECATED_ALIASES = {
    "--chrA": "--ChrA",
    "--chrB": "--ChrB",
    "--hor-chr-list": "--chr-list",
    "--hor-threshold": "--threshold",
    "--hor_threshold": "--threshold",
    "--hor-min-len": "--min-len",
    "--hor_min_len": "--min-len",
}


def warn_deprecated_aliases(argv: list[str]) -> None:
    """Emit a one-line nudge for each redundant alias present in `argv`. Purely
    advisory — the aliases still work; this just points at the canonical name."""
    seen: set[str] = set()
    for tok in argv:
        # handle both `--alias value` and `--alias=value`
        name = tok.split("=", 1)[0]
        canonical = _DEPRECATED_ALIASES.get(name)
        if canonical and name not in seen:
            seen.add(name)
            log.warn(f"{name} is deprecated; use {canonical} instead "
                     "(the old name still works for now)")


# ---------------------------------------------------------------------------
# Main-pipeline integration: HOR options on the top-level parser
# ---------------------------------------------------------------------------

def add_hor_arguments(p: argparse.ArgumentParser) -> None:
    """Add the (all ``--hor-``-prefixed) HOR options to the main pipeline parser.
    Passing ``--hor-chr-list`` or ``--hor-ChrA`` turns HOR detection on.

    The ``--hor-`` prefix is deliberate: it makes clear that these options
    configure a *separate second stage* that runs only after array identification
    finishes — they do not change the tandem-repeat steps above."""
    g = p.add_argument_group(
        "HOR detection — a SEPARATE second stage",
        "These options configure higher-order-repeat detection, which runs only "
        "AFTER array identification finishes. They do not affect the tandem-repeat "
        "steps above. Pass --hor-chr-list (or --hor-ChrA) to turn the stage on. "
        "To run HOR detection on its own on an existing table, see `trash-py hor --help`.")
    g.add_argument("--hor-chr-list", dest="chr_list",
                   help="comma-separated sequence IDs — detect HORs (self-comparison) "
                        "for each. The common way to turn the HOR stage on.")
    g.add_argument("--hor-ChrA", dest="chrA",
                   help="HOR stage: single sequence ID (region A if --hor-ChrB given)")
    g.add_argument("--hor-ChrB", dest="chrB",
                   help="HOR stage: region B sequence ID — enables cross-region comparison")
    g.add_argument("--hor-class", dest="hor_class",
                   help="HOR stage: repeat class to analyse (default: most abundant on the targets)")
    g.add_argument("--hor-classB", dest="classB",
                   help="HOR stage: region-B repeat class (default: same as --hor-class)")
    g.add_argument("--hor-threshold", dest="hor_threshold", type=int, default=4,
                   help="HOR stage: divergence threshold %% (default 4)")
    g.add_argument("--hor-min-len", dest="hor_min_len", type=int, default=3,
                   help="HOR stage: minimum HOR length in repeat units (default 3)")
    g.add_argument("--hor-genomeA", dest="genomeA", default="A",
                   help="HOR stage: label for genome A")
    g.add_argument("--hor-genomeB", dest="genomeB", default="B",
                   help="HOR stage: label for genome B")
    g.add_argument("--hor-sweep", dest="sweep", action="store_true",
                   help="HOR stage: also write an interactive HTML dot-plot with a "
                        "threshold slider, sweeping thresholds 1..--hor-sweep-max "
                        "(self-comparison only)")
    g.add_argument("--hor-sweep-max", dest="sweep_max", type=int, default=30,
                   help="HOR stage: top threshold %% for --hor-sweep (default 30)")
    g.add_argument("--no-hor-plot", dest="no_plot", action="store_true",
                   help="HOR stage: skip the HOR line plot")
    g.add_argument("--no-hor-saveR", dest="no_saveR", action="store_true",
                   help="HOR stage (--hor-ChrB mode only): skip writing repeats_with_hors")


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
        sweep=getattr(args, "sweep", False),
        sweep_max=getattr(args, "sweep_max", 30),
        chrA_flag="--hor-ChrA", chrB_flag="--hor-ChrB", chr_list_flag="--hor-chr-list",
    )


# ---------------------------------------------------------------------------
# Standalone subcommand: `trash-py hor <repeats_with_seq.csv> ...`
# ---------------------------------------------------------------------------

_HOR_EPILOG = """\
examples:
  # detect HORs on one or more sequences (self-comparison)
  trash-py hor out/genome.fasta_repeats_with_seq.csv --chr-list Chr1,Chr2 -o out

  # a single sequence, explicit repeat class
  trash-py hor out/genome.fasta_repeats_with_seq.csv --ChrA Chr1 -c 178_1 -o out

  # cross-region comparison (region A vs region B)
  trash-py hor out/genome.fasta_repeats_with_seq.csv --ChrA Chr1 --ChrB Chr2 -o out

Legacy HORT.R flags (-r, --output_folder, --chrA, -m/--method, -s/--saveR,
--plot_simple, …) are still accepted for drag-and-drop compatibility.
"""


def build_hor_parser() -> argparse.ArgumentParser:
    """The standalone HOR parser. Accepts the full upstream ``HORT.R`` getopt
    surface (short flags and the underscore long-names) so that replacing
    ``Rscript HORT.R`` with ``trash-py hor`` in an existing script Just Works,
    while also accepting a positional repeats file and a ``--chr-list``."""
    p = argparse.ArgumentParser(
        prog="trash-py hor",
        description="TRASH HOR module — detect higher-order repeats in a repeat table.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=_HOR_EPILOG,
    )
    # The repeats table: positional (preferred) or the hidden HORT.R -r/--repeats.
    p.add_argument("repeats_pos", nargs="?", type=Path, metavar="repeats",
                   help="repeat table with a sequence column (…_repeats_with_seq.csv)")
    p.add_argument("-o", "--output", dest="output", type=Path, default=Path.cwd(),
                   help="output folder (default: current directory)")
    p.add_argument("-c", "--class", dest="hor_class", default=None,
                   help="repeat class to analyse (e.g. 178_1); "
                        "default: the most abundant class on the target sequences")
    p.add_argument("--chr-list", dest="chr_list",
                   help="comma-separated sequence IDs — one self-comparison per ID")
    p.add_argument("-A", "--ChrA", dest="chrA",
                   help="single sequence ID (region A if --ChrB is also given)")
    p.add_argument("-B", "--ChrB", dest="chrB",
                   help="region B sequence ID — enables cross-region comparison")
    p.add_argument("-C", "--classB", dest="classB",
                   help="region-B repeat class (default: same as --class)")
    p.add_argument("-b", "--repeatsB", dest="repeatsB", type=Path,
                   help="region-B repeats table (default: same as the main repeats file)")
    p.add_argument("-t", "--threshold", dest="hor_threshold", type=int, default=4,
                   help="divergence threshold %% (default 4)")
    p.add_argument("-l", "--min-len", dest="hor_min_len", type=int, default=3,
                   help="minimum HOR length in repeat units (default 3)")
    p.add_argument("-g", "--genomeA", dest="genomeA", default="A",
                   help="label for genome A (--ChrB mode)")
    p.add_argument("-G", "--genomeB", dest="genomeB", default="B",
                   help="label for genome B (--ChrB mode)")
    p.add_argument("--no-plot", dest="no_plot", action="store_true",
                   help="skip the HOR plot")
    p.add_argument("--sweep", dest="sweep", action="store_true",
                   help="also write an interactive HTML dot-plot with a threshold "
                        "slider, sweeping thresholds 1..--sweep-max (self-comparison only)")
    p.add_argument("--sweep-max", dest="sweep_max", type=int, default=30,
                   help="top threshold %% for --sweep (default 30)")
    p.add_argument("--no-saveR", dest="no_saveR", action="store_true",
                   help="--ChrB mode only: skip writing repeats_with_hors")
    p.add_argument("-T", "--threads", type=int, default=None,
                   help="MAFFT threads (default: all available cores)")
    p.add_argument("-q", "--quiet", action="store_true", help="suppress progress output")

    # Hidden aliases: the upstream HORT.R getopt surface plus earlier trash-py
    # spellings. All still parse (drag-and-drop compatibility) but are kept out
    # of --help to reduce clutter. default=SUPPRESS so they don't clobber the
    # canonical option's default when absent.
    S = argparse.SUPPRESS
    p.add_argument("-r", "--repeats", dest="repeats", type=Path, default=None, help=S)
    p.add_argument("--output_folder", dest="output", type=Path, default=S, help=S)
    p.add_argument("--hor-threshold", "--hor_threshold", dest="hor_threshold",
                   type=int, default=S, help=S)
    p.add_argument("--hor-min-len", "--hor_min_len", dest="hor_min_len",
                   type=int, default=S, help=S)
    p.add_argument("--hor-chr-list", dest="chr_list", default=S, help=S)
    p.add_argument("--chrA", dest="chrA", default=S, help=S)
    p.add_argument("--chrB", dest="chrB", default=S, help=S)
    p.add_argument("-m", "--method", type=int, choices=(1, 2), default=None, help=S)
    p.add_argument("-s", "--saveR", dest="saveR", choices=("y", "n"), default=None, help=S)
    p.add_argument("-p", "--plot-simple", "--plot_simple", dest="plot_simple",
                   choices=("y", "n"), default=None, help=S)
    return p


def run_hor_cli(ns: argparse.Namespace,
                parser: argparse.ArgumentParser | None = None) -> int:
    """Entry point for the standalone `trash-py hor` subcommand."""
    repeats = ns.repeats or ns.repeats_pos
    if repeats is None:
        # No table at all: the user most likely just typed `trash-py hor` — show
        # help rather than a terse one-liner.
        if parser is not None:
            parser.print_help(sys.stderr)
        else:
            print("no repeats table given (pass it positionally or with -r/--repeats)",
                  file=sys.stderr)
        return 2
    if not repeats.exists():
        print(f"repeats file not found: {repeats}", file=sys.stderr)
        return 2
    if not (ns.chr_list or ns.chrA):
        print("specify --chr-list, or --ChrA (optionally with --ChrB)", file=sys.stderr)
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
        sweep=ns.sweep,
        sweep_max=ns.sweep_max,
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
    top_class, top_n, top_bp = ranked[0]
    runners = ", ".join(f"{c} ({bp / 1e6:.2f} Mb)" for c, _, bp in ranked[1:4])
    log.info(f"auto-selected class '{top_class}' ({top_bp / 1e6:.2f} Mb of repeats, "
             f"{top_n} monomers, on {len(seq_ids)} "
             f"sequence{'s' if len(seq_ids) != 1 else ''})"
             + (f"; next by coverage: {runners}" if runners else ""))
    return top_class


def _run(*, repeats: Path, repeatsB: Path, output: Path, chr_list: str | None,
         chrA: str | None, chrB: str | None, hor_class: str | None, classB: str | None,
         hor_threshold: int, hor_min_len: int, genomeA: str, genomeB: str,
         make_plot: bool, saveR: bool, threads: int = 1,
         sweep: bool = False, sweep_max: int = 30,
         chrA_flag: str = "--ChrA", chrB_flag: str = "--ChrB",
         chr_list_flag: str = "--chr-list") -> int:
    if not chr_list and not chrA:
        print(f"specify {chr_list_flag}, or {chrA_flag} (optionally with {chrB_flag})",
              file=sys.stderr)
        return 2

    # cross-region comparison (method 2)
    if chrB:
        if sweep:
            log.warn("--sweep applies to self-comparison only; ignoring it for the "
                     f"{chrB_flag} cross-region comparison")
        if not chrA:
            print(f"{chrB_flag} requires {chrA_flag}", file=sys.stderr)
            return 2
        missing = ([chrA] if chrA not in available_seqids(repeats) else []) \
            + ([chrB] if chrB not in available_seqids(repeatsB) else [])
        resolved = _resolve_class(repeats, {chrA}, hor_class)
        if missing:
            cls_str = f" of class '{resolved}'" if resolved else ""
            log.error(f"skipping all HOR identification{cls_str}: sequence ID(s) "
                      f"not in the repeats table: {_truncate(missing)}")
            return 2
        if not resolved:
            log.error("skipping all HOR identification: could not determine a repeat class (empty table?)")
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

    # Resolve the repeat class (species of HOR) to analyse
    resolved = _resolve_class(repeats, set(chrs), hor_class)

    # Check requested IDs against what's actually in the table, up front, and
    # hard-fail (running nothing) if ANY are missing — so the user can't
    # silently miss sequences. The list is truncated to stay readable.
    missing = [c for c in chrs if c not in available_seqids(repeats)]
    if missing:
        cls_str = f" of class '{resolved}'" if resolved else ""
        log.error(f"skipping all HOR identification{cls_str}: the following {len(missing)} "
                  f"requested sequence ID(s) are not in the repeats table: {_truncate(missing)}")
        return 2

    if not resolved:
        log.error("skipping all HOR identification: could not determine a repeat class (empty table?)")
        return 2
    args = HorArgs(repeats=repeats, output_folder=output, hor_class=resolved,
                   hor_threshold=hor_threshold, hor_min_len=hor_min_len, make_plot=make_plot,
                   threads=threads, sweep=sweep, sweep_max=sweep_max)

    failures = 0
    for chrom in chrs:
        try:
            run_hor_single(args, chrom)
        except Exception as e:
            log.warn(f"{chrom}: {e}")
            failures += 1
    return 1 if failures else 0
