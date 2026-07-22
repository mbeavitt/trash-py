"""HOR (higher-order repeat) detection pipeline — Python port of `HORT.R`.

This mirrors the upstream TRASH HOR module: it filters a repeat table to one
class/sequence, aligns the repeats with MAFFT, runs the HOR scan (the native
`_ext.find_hors`, a byte-for-byte reimplementation of the `HOR.V3.3` binary),
then writes the annotated HOR table and per-repeat "repeats_with_hors" table.

Two modes, exactly as upstream:
* **method 1** (``hort``): all repeats of one class on one sequence, compared
  against themselves — the common case, driven by ``--chr-list`` / ``--ChrA``.
* **method 2** (``horb``): repeats of one region vs another (``--ChrB``),
  split-and-compare.

The HOR tables are reproduced byte-for-byte against the reference tool. See
``docs/HOR_source_bugs.md`` for the upstream quirks we deliberately preserve.
"""
from __future__ import annotations

import math
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

from . import _log as log
from . import _ext
from .hor_bins import calc_plot_regs
from .io_csv import _format


# Column layouts (match R's read.csv/write.csv output exactly).
REPEATS_IN_COLUMNS = [
    "seqID", "arrayID", "start", "end", "strand", "score", "eval",
    "width", "class", "score_template", "sequence",
]
HORS_COLUMNS = [
    "start_A", "end_A", "start_B", "end_B", "direction.1.para_2.perp.",
    "total_variant", "start.A.bp", "start.B.bp", "end.A.bp", "end.B.bp",
    "chrA", "chrB", "block.size.in.units", "block.A.size.bp",
    "block.B.size.bp", "SNV_per_kbp",
]
REPEATS_WITH_HORS_COLUMNS = REPEATS_IN_COLUMNS + [
    "start.adjusted", "hors_formed_count", "hors_formed_tot_rep_normalised",
]

WINDOW = 100000  # 100 Kbp binning window for plot regions


@dataclass
class Repeat:
    raw: dict[str, str]      # original CSV cell strings, for byte-faithful passthrough
    seqID: str
    start: int
    end: int
    width: int
    strand: str              # normalised to "1" / "2"
    sequence: str


# --------------------------------------------------------------------------
# CSV I/O
# --------------------------------------------------------------------------

def read_repeats(path: Path) -> list[Repeat]:
    import csv

    out: list[Repeat] = []
    with Path(path).open(newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            strand = row["strand"]
            strand = "1" if strand == "+" else "2" if strand == "-" else strand
            out.append(
                Repeat(
                    raw=row,
                    seqID=row["seqID"],
                    start=int(row["start"]),
                    end=int(row["end"]),
                    width=int(row["width"]),
                    strand=strand,
                    sequence=row["sequence"],
                )
            )
    return out


def class_abundance(path: Path, seq_ids: set[str] | None = None) -> "list[tuple[str, int, int]]":
    """Rank repeat classes (optionally restricted to `seq_ids`) by total repeat
    bp — i.e. Σ monomer width, a proxy for how much sequence the class covers.

    Returns `(class, monomer_count, total_bp)` sorted by `total_bp` (then count)
    descending. Used to auto-pick the *major* class for HOR detection: ranking by
    coverage rather than raw count avoids picking a short microsatellite that
    happens to have the most monomers over the genome's real satellite."""
    import csv
    from collections import defaultdict

    count: dict[str, int] = defaultdict(int)
    total_bp: dict[str, int] = defaultdict(int)
    with Path(path).open(newline="") as f:
        for row in csv.DictReader(f):
            if seq_ids is None or row["seqID"] in seq_ids:
                cls = row["class"]
                count[cls] += 1
                total_bp[cls] += int(row["width"])
    return sorted(((c, count[c], total_bp[c]) for c in count),
                  key=lambda t: (-t[2], -t[1], t[0]))


def available_seqids(path: Path) -> set[str]:
    """The set of sequence IDs present in a repeats table."""
    import csv

    with Path(path).open(newline="") as f:
        return {row["seqID"] for row in csv.DictReader(f)}


def _median(values: Sequence[int]) -> float:
    s = sorted(values)
    n = len(s)
    if n == 0:
        return 0.0
    mid = n // 2
    if n % 2:
        return float(s[mid])
    return (s[mid - 1] + s[mid]) / 2.0


def _fmt_cell(raw: str) -> str:
    """Format an original numeric cell the way R's read.csv -> write.csv would:
    integers unchanged, decimals re-rendered at 15 significant digits."""
    if any(c in raw for c in ".eE") and _looks_numeric(raw):
        return _format(float(raw))
    if _looks_int(raw):
        return raw
    return f'"{raw}"'


def _looks_int(s: str) -> bool:
    try:
        int(s)
        return True
    except ValueError:
        return False


def _looks_numeric(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def _write_repeats_with_hors_csv(path: Path, repeats: list[Repeat],
                                 start_adjusted: list[float],
                                 hors_count: list[int], n_rep: int) -> None:
    """repeats_with_hors — R's `write.csv(..., row.names = FALSE)`."""
    with Path(path).open("w", newline="") as f:
        f.write(",".join(f'"{c}"' for c in REPEATS_WITH_HORS_COLUMNS))
        f.write("\n")
        for k, rep in enumerate(repeats):
            cells = [
                f'"{rep.seqID}"',
                _fmt_cell(rep.raw["arrayID"]),
                _fmt_cell(rep.raw["start"]),
                _fmt_cell(rep.raw["end"]),
                f'"{rep.strand}"',
                _fmt_cell(rep.raw["score"]),
                _fmt_cell(rep.raw["eval"]),
                _fmt_cell(rep.raw["width"]),
                f'"{rep.raw["class"]}"',
                _fmt_cell(rep.raw["score_template"]),
                f'"{rep.sequence}"',
                _format(float(start_adjusted[k])),
                str(hors_count[k]),
                _format(hors_count[k] / n_rep),
            ]
            f.write(",".join(cells))
            f.write("\n")


# --------------------------------------------------------------------------
# Alignment
# --------------------------------------------------------------------------

def _wrap60(seq: str) -> str:
    return "\n".join(seq[i:i + 60] for i in range(0, len(seq), 60))


def align_repeats(repeats: list[Repeat], out_dir: Path, name: str, threads: int = 1) -> Path:
    """Reproduce `write_align_read`: write a FASTA (names ``<i>_D<strand>``),
    run ``mafft --retree 2 --inputorder``, and lower-case the result. The HOR
    scanner strips newlines, so wrapping width of the aligned file is
    irrelevant — only order, names and (lower-cased) columns matter.

    `threads` is passed to ``mafft --thread``; it parallelises FFT-NS-2 and
    (verified) leaves the alignment — and therefore the HOR table — byte-for-byte
    identical, so it's a free speedup (MAFFT is ~65-70%% of HOR runtime)."""
    if shutil.which("mafft") is None:
        raise RuntimeError(
            "mafft not found on PATH — install it (e.g. `conda install -c bioconda mafft`)"
        )
    in_fasta = out_dir / f"{name}temp.fasta"
    aligned = out_dir / f"{name}temp.aligned.fasta"
    with in_fasta.open("w") as f:
        for i, rep in enumerate(repeats, start=1):
            f.write(f">{i}_D{rep.strand}\n{_wrap60(rep.sequence)}\n")

    with aligned.open("w") as out:
        proc = log.run_external(
            "mafft",
            ["mafft", "--thread", str(threads), "--retree", "2", "--inputorder", str(in_fasta)],
            stdout=out, stderr=__import__("subprocess").DEVNULL,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"mafft failed (exit {proc.returncode})")
    in_fasta.unlink(missing_ok=True)

    # lower-case sequence lines in place (forceToLower = TRUE upstream)
    tmp = aligned.with_suffix(".lc")
    with aligned.open() as src, tmp.open("w") as dst:
        for line in src:
            dst.write(line if line.startswith(">") else line.lower())
    tmp.replace(aligned)
    if aligned.stat().st_size == 0:
        raise RuntimeError("mafft produced an empty alignment")
    return aligned


# --------------------------------------------------------------------------
# Derived columns
# --------------------------------------------------------------------------

def _stream_hors_to_csv(raw_path: Path, repeats: list[Repeat], out_path: Path,
                        both_blocks: bool, total: int, plot_cap: int = 0
                        ) -> "tuple[list[int], list[list]]":
    """Stream the raw HOR file (6 native int32 per HOR, as written by
    ``_ext.find_hors_stream``) into the final HORs CSV, computing the derived
    bp/size/divergence columns per chunk. Peak memory is O(N + chunk),
    independent of the number of HORs — the point of the exercise when a single
    array yields tens of millions of them.

    Also returns the per-repeat ``hors_formed_count`` (via an O(N) difference
    array) and a strided sample of annotated rows for plotting (≤ ~plot_cap).

    NOTE: `SNV_per_kbp` divides by block A's size twice — reproducing an upstream
    typo (`block.A.size.bp + block.A.size.bp`) instead of A + B. Kept byte-exact
    against the batch path via `io_csv._format`."""
    import numpy as np

    n_rep = len(repeats)
    starts = np.fromiter((r.start for r in repeats), np.int64, n_rep)
    ends = np.fromiter((r.end for r in repeats), np.int64, n_rep)
    qseq = [f'"{r.seqID}"' for r in repeats]
    single_seq = len({r.seqID for r in repeats}) == 1
    const_q = qseq[0] if single_seq else None

    diff = np.zeros(n_rep + 1, dtype=np.int64)          # A-block range increments
    diffB = np.zeros(n_rep + 1, dtype=np.int64) if both_blocks else None
    stride = max(1, total // plot_cap) if plot_cap else 0
    plot_rows: list[list] = []

    CHUNK = 1 << 18   # HOR rows per chunk (bounds the per-chunk memory spike)
    header = ",".join(['""'] + [f'"{c}"' for c in HORS_COLUMNS]) + "\n"
    idx = 0
    with Path(out_path).open("w", newline="") as out, Path(raw_path).open("rb") as raw:
        out.write(header)
        while True:
            buf = raw.read(CHUNK * 24)          # 6 int32 = 24 bytes / HOR
            if not buf:
                break
            a = np.frombuffer(buf, dtype=np.int32).reshape(-1, 6).astype(np.int64)
            sa, ea, sb, eb, dv, tv = (a[:, k] for k in range(6))
            sabp = starts[sa - 1]; sbbp = starts[sb - 1]
            eabp = ends[ea - 1]; ebbp = ends[eb - 1]
            bsu = ea - sa + 1
            ba = eabp - sabp + 1
            bb = ebbp - sbbp + 1
            snvk = 1000.0 * tv / ba              # (ba+ba)/2 == ba; matches the batch path

            # hors_formed_count via difference array (each HOR marks a unit range)
            diff += np.bincount(sa - 1, minlength=n_rep + 1)
            diff -= np.bincount(sa - 1 + bsu, minlength=n_rep + 1)
            if both_blocks:
                diffB += np.bincount(sb - 1, minlength=n_rep + 1)
                diffB -= np.bincount(sb - 1 + bsu, minlength=n_rep + 1)

            saL, eaL, sbL, ebL = sa.tolist(), ea.tolist(), sb.tolist(), eb.tolist()
            dvL, tvL = dv.tolist(), tv.tolist()
            sabpL, sbbpL = sabp.tolist(), sbbp.tolist()
            eabpL, ebbpL = eabp.tolist(), ebbp.tolist()
            bsuL, baL, bbL, snvL = bsu.tolist(), ba.tolist(), bb.tolist(), snvk.tolist()

            lines = []
            for t in range(len(saL)):
                idx += 1
                if single_seq:
                    ca = cb = const_q
                else:
                    ca = qseq[saL[t] - 1]; cb = qseq[sbL[t] - 1]
                lines.append(
                    f'"{idx}",{saL[t]},{eaL[t]},{sbL[t]},{ebL[t]},{dvL[t]},{tvL[t]},'
                    f'{sabpL[t]},{sbbpL[t]},{eabpL[t]},{ebbpL[t]},{ca},{cb},'
                    f'{bsuL[t]},{baL[t]},{bbL[t]},{_format(snvL[t])}'
                )
                if stride and (idx - 1) % stride == 0:
                    plot_rows.append([saL[t], eaL[t], sbL[t], ebL[t], dvL[t], tvL[t],
                                      sabpL[t], sbbpL[t], eabpL[t], ebbpL[t], ca, cb,
                                      bsuL[t], baL[t], bbL[t], snvL[t]])
            out.write("\n".join(lines))
            out.write("\n")

    counts = np.cumsum(diff[:n_rep])
    if both_blocks:
        counts = counts + np.cumsum(diffB[:n_rep])
    return counts.tolist(), plot_rows


def _start_adjusted(repeats: list[Repeat]) -> tuple[list[float], "object"]:
    starts = [r.start for r in repeats]
    bins = calc_plot_regs(starts, WINDOW)
    adjusted = [float(r.start) for r in repeats]
    for bs, be, corr in zip(bins.starts, bins.ends, bins.corrections):
        for k, r in enumerate(repeats):
            if bs <= r.start < be:
                adjusted[k] = r.start - corr
    return adjusted, bins


# --------------------------------------------------------------------------
# Orchestration
# --------------------------------------------------------------------------

@dataclass
class HorArgs:
    repeats: Path
    output_folder: Path
    hor_class: str
    hor_threshold: int = 4
    hor_min_len: int = 3
    make_plot: bool = True
    threads: int = 1


def run_hor_single(args: HorArgs, chrA: str, rng_tag: float = 0.0) -> None:
    """method 1: one class on one sequence, self-comparison (the ``hort`` path)."""
    log.section(f"HOR: class {args.hor_class} on {chrA}")
    repeats = [r for r in read_repeats(args.repeats)
               if r.raw["class"] == args.hor_class]
    if len(repeats) < args.hor_min_len * 2:
        raise ValueError(
            f"Not enough repeats of class {args.hor_class}: {len(repeats)}")
    repeats = [r for r in repeats if r.seqID == chrA]
    if len(repeats) < args.hor_min_len * 2:
        raise ValueError(f"Not enough repeats on sequence {chrA}: {len(repeats)}")
    repeats.sort(key=lambda r: r.start)

    threshold_snv = math.floor(args.hor_threshold * _median([r.width for r in repeats]) / 100)
    log.detail(f"{len(repeats)} repeats; threshold_SNV={threshold_snv}, min_len={args.hor_min_len}")

    name = f"{args.hor_class}__{chrA}_{rng_tag}_"
    aligned = align_repeats(repeats, args.output_folder, name, threads=args.threads)
    log.tool_summary("mafft")

    raw = args.output_folder / f"{name}raw.hors"
    total = _ext.find_hors_stream(str(aligned), str(raw), 1, threshold_snv, args.hor_min_len, 1)
    log.detail(f"HORs identified: {total}")
    aligned.unlink(missing_ok=True)

    if total <= 1:
        log.detail("No HORs identified")
        raw.unlink(missing_ok=True)
        return

    out_hors = args.output_folder / f"HORs_{args.hor_class}_{chrA}.csv"
    plot_cap = 400_000 if args.make_plot else 0
    counts, plot_rows = _stream_hors_to_csv(
        raw, repeats, out_hors, both_blocks=False, total=total, plot_cap=plot_cap)
    raw.unlink(missing_ok=True)

    adjusted, bins = _start_adjusted(repeats)
    out_rep = args.output_folder / f"repeats_with_hors_{args.hor_class}_{chrA}.csv"
    _write_repeats_with_hors_csv(out_rep, repeats, adjusted, counts, len(repeats))

    if args.make_plot:
        try:
            from .hor_plot import plot_hors
            plot_hors(plot_rows,
                      args.output_folder / f"HORs_dotplot_{args.hor_class}_{chrA}.png",
                      chrA, args.hor_class, args.hor_threshold)
        except Exception as e:  # plots are best-effort
            log.warn(f"HOR plot failed: {e}")
    log.detail(f"wrote {out_hors.name}, {out_rep.name}")


def run_hor_pair(args: HorArgs, chrA: str, chrB: str, classB: str,
                 repeatsB: Path, genomeA: str = "A", genomeB: str = "B",
                 saveR: bool = True) -> None:
    """method 2: repeats of (classA, chrA) vs (classB, chrB) — the ``horb`` path."""
    log.section(f"HOR: {args.hor_class}/{chrA} vs {classB}/{chrB}")
    repeats = [r for r in read_repeats(args.repeats)
               if r.raw["class"] == args.hor_class and r.seqID == chrA]
    if not repeats:
        raise ValueError(f"No repeats of class {args.hor_class} on {chrA}")
    repeatsB_all = [r for r in read_repeats(repeatsB)
                    if r.raw["class"] == classB and r.seqID == chrB]
    if not repeatsB_all:
        raise ValueError(f"No repeatsB of class {classB} on {chrB}")

    split_after = len(repeats)
    combined = repeats + repeatsB_all
    threshold_snv = math.floor(args.hor_threshold * _median([r.width for r in combined]) / 100)

    name = f"{args.hor_class}_{classB}_{chrA}_{chrB}_{genomeA}_{genomeB}_"
    aligned = align_repeats(combined, args.output_folder, name, threads=args.threads)
    log.tool_summary("mafft")

    raw = args.output_folder / f"{name}raw.hors"
    total = _ext.find_hors_stream(str(aligned), str(raw), split_after, threshold_snv,
                                  args.hor_min_len, 2)
    log.detail(f"HORs identified: {total}")
    aligned.unlink(missing_ok=True)

    suffix = f"{args.hor_class}_{classB}_{chrA}_{chrB}_{genomeA}_{genomeB}"
    if total <= 1:
        log.detail("No HORs identified")
        raw.unlink(missing_ok=True)
        _write_summary(args.output_folder / f"summary_of_hors_{suffix}.csv",
                       genomeA, genomeB, chrA, chrB, args.hor_class, classB,
                       split_after, len(combined) - split_after, 0, 0)
        return

    counts, _ = _stream_hors_to_csv(
        raw, combined, args.output_folder / f"HORs_{suffix}.csv",
        both_blocks=True, total=total, plot_cap=0)
    raw.unlink(missing_ok=True)

    if saveR:
        adjusted, _ = _start_adjusted(combined)
        _write_repeats_with_hors_csv(
            args.output_folder / f"repeats_with_hors_{suffix}.csv",
            combined, adjusted, counts, len(combined))

    a_in_b = sum(1 for c in counts[:split_after] if c > 0)
    b_in_a = sum(1 for c in counts[split_after - 1:] if c > 0)
    _write_summary(args.output_folder / f"summary_of_hors_{suffix}.csv",
                   genomeA, genomeB, chrA, chrB, args.hor_class, classB,
                   split_after, len(combined) - split_after, a_in_b, b_in_a)


def _write_summary(path: Path, genomeA, genomeB, chrA, chrB, repA, repB,
                   repA_total, repB_total, a_in_b, b_in_a) -> None:
    cols = ["genomeA", "genomeB", "chrA", "chrB", "repA", "repB",
            "repA_total", "repB_total", "A_repeats_found_in_B", "B_repeats_found_in_A"]
    with Path(path).open("w", newline="") as f:
        f.write(",".join(f'"{c}"' for c in cols) + "\n")
        vals = [f'"{genomeA}"', f'"{genomeB}"', f'"{chrA}"', f'"{chrB}"',
                f'"{repA}"', f'"{repB}"', str(repA_total), str(repB_total),
                str(a_in_b), str(b_in_a)]
        f.write(",".join(vals) + "\n")
