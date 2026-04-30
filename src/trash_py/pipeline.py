"""End-to-end pipeline orchestrator."""
from __future__ import annotations

import shutil
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

from . import _log as log
from . import __version__
from .arrays import ArrayRow, split_and_check_arrays
from .classify import classify_arrays
from .io_csv import write_csv_r_style
from .io_fasta import read_fasta_and_list
from .io_gff import export_gff
from .mapping import (
    clean_overlaps,
    handle_edge_repeat,
    map_repeats,
    rescore_repeats,
)
from .merge_windows import merge_windows
from .shift import shift_and_compare
from .summarise import (
    ARRAYS_COLUMNS,
    REPEATS_COLUMNS,
    REPEATS_WITH_SEQ_COLUMNS,
    append_sequence_column,
    strip_repeats_columns,
    summarise_arrays,
)
from .window_score import sequence_window_score


STAGE8_COLUMNS = [
    "start", "end", "seqID", "numID", "score", "top_N", "top_5_N",
    "representative", "class", "array_num_ID",
]

ARRAYS_PER_CHUNK = 100
KMER = 10


# External programs the pipeline shells out to. Both must be on PATH.
# Install via bioconda (works on macOS, Linux, Windows-WSL):
#     conda install -c bioconda clustalo hmmer
_EXTERNAL_TOOLS: dict[str, str] = {
    "clustalo": "Clustal Omega — `conda install -c bioconda clustalo`",
    "nhmmer": "HMMER (provides nhmmer) — `conda install -c bioconda hmmer`",
}


def check_external_tools(tools: dict[str, str] = _EXTERNAL_TOOLS) -> list[str]:
    """Return the names of tools missing from PATH (empty list if all present)."""
    return [name for name in tools if shutil.which(name) is None]


def _require_external_tools() -> None:
    missing = check_external_tools()
    if not missing:
        return
    lines = [
        "trash-py needs the following external programs on your PATH:",
        "",
    ]
    for name in missing:
        lines.append(f"  - {name}: not found")
        lines.append(f"      install: {_EXTERNAL_TOOLS[name]}")
    print("\n".join(lines), file=sys.stderr)
    raise SystemExit(2)


def _log_genome_summary(fasta: list[tuple[str, str]]) -> None:
    name_w = min(max((len(n) for n, _ in fasta), default=0), 32)
    if len(fasta) <= 12:
        for name, seq in fasta:
            log.detail(f"{name:<{name_w}}  {len(seq):>15,} bp")
    else:
        # Show the longest 10 + a "(N more)" line.
        ranked = sorted(fasta, key=lambda x: -len(x[1]))
        for name, seq in ranked[:10]:
            log.detail(f"{name:<{name_w}}  {len(seq):>15,} bp")
        log.detail(f"({len(fasta) - 10} more sequence{'s' if len(fasta) - 10 != 1 else ''})")
    total = sum(len(s) for _, s in fasta)
    plural = "s" if len(fasta) != 1 else ""
    log.info()
    if total >= 1_000_000:
        log.detail(f"{len(fasta)} sequence{plural}, {total:,} bp ({total / 1e6:.2f} Mbp)")
    else:
        log.detail(f"{len(fasta)} sequence{plural}, {total:,} bp")


def _format_pct(part: int, whole: int) -> str:
    return f"{100.0 * part / whole:.1f}%" if whole else "0.0%"


def _format_span(bp: int) -> str:
    if bp >= 1_000_000:
        return f"{bp / 1e6:.2f}Mbp"
    if bp >= 1_000:
        return f"{bp / 1e3:.1f}kbp"
    return f"{bp}bp"


def _log_class_breakdown(classarrays: list[dict[str, Any]], top_n: int = 8) -> None:
    """Tabulate the dominant classes: name, member count, median rep width,
    and total genomic span (sum of array end-start). Ranked by span desc."""
    from collections import defaultdict
    from statistics import median

    by_class: dict[str, list[dict[str, Any]]] = defaultdict(list)
    for r in classarrays:
        by_class[r["class"]].append(r)

    def class_span(members: list[dict[str, Any]]) -> int:
        return sum(int(m["end"]) - int(m["start"]) + 1 for m in members)

    ranked = sorted(by_class.items(), key=lambda kv: -class_span(kv[1]))
    if not ranked:
        return

    rows: list[tuple[str, int, str, str]] = []
    for cls, members in ranked[:top_n]:
        widths = [len(m["representative"]) for m in members]
        rows.append((
            cls,
            len(members),
            f"{int(median(widths))}bp",
            _format_span(class_span(members)),
        ))

    name_w = min(max(max(len(c) for c, *_ in rows), len("class")), 24)
    count_w = max(max(len(str(n)) for _, n, *_ in rows), len("count"))
    width_w = max(max(len(w) for _, _, w, _ in rows), len("width"))
    span_w = max(max(len(s) for *_, s in rows), len("span"))

    log.detail("top classes by span:")
    log.detail(
        f"    {'class':<{name_w}}  {'count':>{count_w}}  "
        f"{'width':>{width_w}}  {'span':>{span_w}}"
    )
    for cls, n, w, s in rows:
        log.detail(
            f"    {cls:<{name_w}}  {n:>{count_w}}  {w:>{width_w}}  {s:>{span_w}}"
        )
    if len(ranked) > top_n:
        log.detail(
            f"    (+{len(ranked) - top_n} more "
            f"class{'es' if len(ranked) - top_n != 1 else ''})"
        )


def _worker_init() -> None:
    """ProcessPool worker bootstrap: silence per-worker logging and switch
    `_log.run_external` into worker mode so the one-shot `running {tool}...`
    line is emitted only by the parent."""
    log.configure(quiet=True)
    log.set_worker_mode(True)


def _split_arrays_worker(
    task: tuple[int, str, str, int, int, int, int],
) -> tuple[list[ArrayRow], dict]:
    region_start, region_seq, seqID, numID, max_repeat, min_repeat, kmer = task
    rows = split_and_check_arrays(
        region_start=region_start,
        sequence=region_seq,
        seqID=seqID,
        numID=numID,
        max_repeat=max_repeat,
        min_repeat=min_repeat,
        kmer=kmer,
    )
    return rows, log.pop_stats()


def _map_array(task: tuple[dict, str, str, str, int, dict, str]) -> list[dict]:
    (arr, seqID, array_sequence, sequence_substring,
     adjust_start, templates_by_name, output_folder) = task

    rows = map_repeats(
        representative=arr["representative"],
        top_N=int(arr["top_N"]),
        array_sequence=array_sequence,
        seqID=seqID,
        arrayID=int(arr["array_num_ID"]),
        start_offset=int(arr["start"]),
        output_folder=Path(output_folder),
    )
    if len(rows) < 2:
        return []

    rows = clean_overlaps(
        rows=rows,
        representative_len=int(arr["top_N"]),
        arr_class=arr["class"],
    )
    if len(rows) < 3:
        return []

    rows = rescore_repeats(
        rows=rows,
        arr_representative=arr["representative"],
        arr_class=arr["class"],
        top_N=int(arr["top_N"]),
        array_sequence=array_sequence,
        array_start=int(arr["start"]),
        templates=templates_by_name,
    )

    rows = handle_edge_repeat(
        rows=rows,
        sequence_vector=sequence_substring,
        sequence_vector_start=adjust_start,
        template_sequence=templates_by_name.get(arr["class"], ""),
    )
    if len(rows) < 3:
        return []

    return [{k: r[k] for k in (*REPEATS_COLUMNS, "representative")} for r in rows]


def _map_array_worker(
    task: tuple[dict, str, str, str, int, dict, str],
) -> tuple[list[dict], dict]:
    out = _map_array(task)
    return out, log.pop_stats()


def run_pipeline(args: Any) -> None:
    """Drive the pipeline. `args` must expose `.fasta`, `.output`,
    `.max_rep_size`, `.min_rep_size`, and (optionally) `.templates`,
    `.processes`.
    """
    _require_external_tools()
    processes = max(1, int(getattr(args, "processes", 1) or 1))
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_name = Path(args.fasta).name
    window_size = round((args.max_rep_size + KMER) * 1.1)

    log.header(f"trash-py v{__version__}")
    log.detail(f"input:  {args.fasta}")
    log.detail(f"output: {output_dir}/")

    # Load fasta + dedupe sequence headers (suffix any duplicates with 1, 2, ...)
    fasta = read_fasta_and_list(Path(args.fasta))
    if not fasta:
        log.warn(f"empty fasta: {args.fasta}")
        return
    log.info()
    _log_genome_summary(fasta)

    names = [name for name, _ in fasta]
    if len(names) != len(set(names)):
        counts: dict[str, int] = {}
        for n in names:
            counts[n] = counts.get(n, 0) + 1
        next_idx: dict[str, int] = {n: 1 for n, c in counts.items() if c > 1}
        renamed: list[tuple[str, str]] = []
        for name, seq in fasta:
            if counts[name] > 1:
                new = f"{name}{next_idx[name]}"
                next_idx[name] += 1
                renamed.append((new, seq))
            else:
                renamed.append((name, seq))
        fasta = renamed

    # Per-sequence kmer window scoring.
    log.section("scanning windows for repeat content")
    all_scores: list[list[float]] = []
    for _, seq in fasta:
        all_scores.append(sequence_window_score(seq, window_size, KMER))

    # Merge low-singleton windows into repetitive regions.
    rows = []
    per_seq_regions: list[list] = []
    for idx, ((name, seq), scores) in enumerate(zip(fasta, all_scores), start=1):
        seq_regions = list(merge_windows(scores, window_size, len(seq)))
        per_seq_regions.append(seq_regions)
        for region in seq_regions:
            rows.append([region.start, region.end, region.score, name, idx])

    write_csv_r_style(
        output_dir / f"{fasta_name}_regarrays.csv",
        columns=["starts", "ends", "scores", "seqID", "numID"],
        rows=rows,
    )

    region_total_bp = sum(r.end - r.start + 1 for regs in per_seq_regions for r in regs)
    region_count = sum(len(regs) for regs in per_seq_regions)
    genome_total = sum(len(s) for _, s in fasta)
    region_plural = "s" if region_count != 1 else ""
    log.detail(f"{region_count} repetitive region{region_plural} found")
    log.detail(f"{region_total_bp:,} bp covered ({_format_pct(region_total_bp, genome_total)} of input)")
    log.elapsed_marker()

    # Per-region: split into individual arrays + extract a representative.
    log.section("identifying arrays in repetitive regions")
    region_tasks: list[tuple[int, str, str, int, int, int, int]] = []
    for (name, seq), seq_regions, idx in zip(
        fasta, per_seq_regions, range(1, len(fasta) + 1)
    ):
        for region in seq_regions:
            region_seq = seq[region.start - 1:region.end]
            region_tasks.append((
                region.start, region_seq, name, idx,
                args.max_rep_size, args.min_rep_size, KMER,
            ))

    arr_rows: list[list] = []

    def _emit_array_rows(arrs: list[ArrayRow]) -> None:
        for arr in arrs:
            arr_rows.append([
                arr.start, arr.end, arr.seqID, arr.numID,
                arr.score, arr.top_N, arr.top_5_N, arr.representative,
            ])

    if processes > 1 and region_tasks:
        log.announce_tool("clustalo")
        with ProcessPoolExecutor(max_workers=processes, initializer=_worker_init) as ex:
            for arrs, stats in ex.map(_split_arrays_worker, region_tasks):
                log.merge_stats(stats)
                _emit_array_rows(arrs)
    else:
        for region_start, region_seq, name, idx, max_rep, min_rep, kmer in region_tasks:
            _emit_array_rows(split_and_check_arrays(
                region_start=region_start,
                sequence=region_seq,
                seqID=name,
                numID=idx,
                max_repeat=max_rep,
                min_repeat=min_rep,
                kmer=kmer,
            ))

    write_csv_r_style(
        output_dir / f"{fasta_name}_aregarrays.csv",
        columns=["start", "end", "seqID", "numID", "score", "top_N", "top_5_N", "representative"],
        rows=arr_rows,
    )

    arr_plural = "s" if len(arr_rows) != 1 else ""
    log.detail(f"{len(arr_rows)} candidate array{arr_plural}")
    log.tool_summary("clustalo")
    log.elapsed_marker()

    # Canonicalise each representative; if a templates fasta was supplied,
    # also classify each array against the templates.
    templates: list[tuple[str, str]] = []
    if getattr(args, "templates", None) is not None:
        templates = read_fasta_and_list(Path(args.templates))
        if not templates:
            raise SystemExit(f"No templates found in {args.templates}")
    template_names = {n for n, _ in templates}
    templates_by_name = dict(templates)

    log.section("canonicalising array representatives")
    if templates:
        log.detail(
            f"templates: {len(templates)} ({', '.join(n for n, _ in templates)})"
        )
    stage7_rows = []
    for start, end, seqID, numID, score, top_N, top_5_N, representative in arr_rows:
        cls, shifted = shift_and_compare(representative, templates or None)
        stage7_rows.append({
            "start": start, "end": end, "seqID": seqID, "numID": numID,
            "score": score, "top_N": top_N, "top_5_N": top_5_N,
            "representative": shifted, "class": cls,
        })
    if templates:
        matched = sum(1 for r in stage7_rows if r["class"])
        log.detail(f"{matched} of {len(stage7_rows)} arrays matched a template")
    log.elapsed_marker()

    # Class assignment + per-class shift + sort + split off the unclassified.
    log.section("classifying arrays")
    classarrays, no_repeats = classify_arrays(stage7_rows, template_names=template_names)
    classes = sorted({r["class"] for r in classarrays})
    log.detail(
        f"{len(classes)} class{'es' if len(classes) != 1 else ''}; "
        f"{len(classarrays)} arrays classified, {len(no_repeats)} unclassified"
    )
    if classarrays:
        log.info()
        _log_class_breakdown(classarrays)
    log.elapsed_marker()

    def _rows_for_write(rows):
        return [[r[c] for c in STAGE8_COLUMNS] for r in rows]

    write_csv_r_style(
        output_dir / f"{fasta_name}_classarrays.csv",
        columns=STAGE8_COLUMNS,
        rows=_rows_for_write(classarrays),
    )
    write_csv_r_style(
        output_dir / f"{fasta_name}_no_repeats_arrays.csv",
        columns=STAGE8_COLUMNS,
        rows=_rows_for_write(no_repeats),
    )

    if not classarrays:
        log.section("done")
        log.detail("no classifiable arrays — nothing to map")
        log.detail(f"total elapsed: {log.format_elapsed(log.elapsed_since_start())}")
        return

    # Per-array repeat mapping with overlap/gap cleanup, rescore, and
    # edge-repeat refinement.
    log.section("mapping repeats")
    fasta_by_seqID = {name: seq for name, seq in fasta}
    repeats_rows: list[dict] = []

    by_seq: dict[str, list[dict]] = {}
    for arr in classarrays:
        by_seq.setdefault(arr["seqID"], []).append(arr)

    map_tasks: list[tuple[dict, str, str, str, int, dict, str]] = []
    for seqID, arrs_in_seq in by_seq.items():
        fasta_seq = fasta_by_seqID[seqID]
        # Match upstream chunking: chunk_edges = [0, 100, 200, ..., n-1]
        # so the boundary array on each chunk is processed (and appended)
        # twice. The downstream sort+dedup hides this from the output.
        n = len(arrs_in_seq)
        chunk_edges = list(range(0, n, ARRAYS_PER_CHUNK))
        chunk_edges.append(n - 1)
        for j in range(len(chunk_edges) - 1):
            lo = chunk_edges[j]
            hi = chunk_edges[j + 1]
            chunk = arrs_in_seq[lo:hi + 1]
            first_start = int(chunk[0]["start"])
            last_end = int(chunk[-1]["end"])
            sequence_substring = fasta_seq[first_start - 1:last_end - 1]
            adjust_start = first_start - 1

            for arr in chunk:
                if arr["representative"] == "":
                    continue
                a_start = int(arr["start"])
                a_end = int(arr["end"])
                array_sequence = fasta_seq[a_start - 1:a_end]
                map_tasks.append((
                    arr, seqID, array_sequence, sequence_substring,
                    adjust_start, templates_by_name, str(output_dir),
                ))

    if processes > 1 and map_tasks:
        log.announce_tool("nhmmer")
        log.announce_tool("clustalo")
        with ProcessPoolExecutor(max_workers=processes, initializer=_worker_init) as ex:
            for rows, stats in ex.map(_map_array_worker, map_tasks):
                log.merge_stats(stats)
                repeats_rows.extend(rows)
    else:
        for task in map_tasks:
            repeats_rows.extend(_map_array(task))

    log.detail(
        f"{len(repeats_rows):,} repeats found across {len(classarrays)} arrays"
    )
    log.tool_summary("nhmmer")
    log.tool_summary("clustalo")
    log.elapsed_marker()

    # Aggregate into the final arrays table; trim the repeats table.
    arrays_final = summarise_arrays(classarrays, repeats_rows)
    repeats_final = strip_repeats_columns(repeats_rows)

    write_csv_r_style(
        output_dir / f"{fasta_name}_arrays.csv",
        columns=ARRAYS_COLUMNS,
        rows=[[a[c] for c in ARRAYS_COLUMNS] for a in arrays_final],
    )
    export_gff(
        arrays_final,
        output_dir / f"{fasta_name}_arrays.gff",
        seqid=3, source="TRASH", type_="Satellite_array",
        start=1, end=2, score=5,
        attributes=[9, 10, 11],
        attribute_names=["Name=", "Repeat_no=", "Repeat_median_width="],
    )

    write_csv_r_style(
        output_dir / f"{fasta_name}_repeats.csv",
        columns=REPEATS_COLUMNS,
        rows=[[r[c] for c in REPEATS_COLUMNS] for r in repeats_final],
    )
    export_gff(
        repeats_final,
        output_dir / f"{fasta_name}_repeats.gff",
        seqid=1, source="TRASH", type_="Satellite_DNA",
        start=3, end=4, strand=5,
        attributes=[9, 6, 10],
        attribute_names=["Name=", "Arry_EDS=", "Family_EDS="],
    )

    repeats_with_seq = append_sequence_column(repeats_final, fasta_by_seqID)
    write_csv_r_style(
        output_dir / f"{fasta_name}_repeats_with_seq.csv",
        columns=REPEATS_WITH_SEQ_COLUMNS,
        rows=[[r[c] for c in REPEATS_WITH_SEQ_COLUMNS] for r in repeats_with_seq],
    )

    log.section("output")
    expected = [
        f"{fasta_name}_regarrays.csv",
        f"{fasta_name}_aregarrays.csv",
        f"{fasta_name}_classarrays.csv",
        f"{fasta_name}_no_repeats_arrays.csv",
        f"{fasta_name}_arrays.csv",
        f"{fasta_name}_arrays.gff",
        f"{fasta_name}_repeats.csv",
        f"{fasta_name}_repeats.gff",
        f"{fasta_name}_repeats_with_seq.csv",
    ]
    written = [f for f in expected if (output_dir / f).exists()]
    for f in written:
        log.detail(f"{output_dir}/{f}")
    log.info()
    log.detail(f"{len(written)} files in {output_dir}/")

    log.info()
    log.info(f"done in {log.format_elapsed(log.elapsed_since_start())}")
