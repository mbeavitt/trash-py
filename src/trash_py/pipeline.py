"""End-to-end pipeline orchestrator."""
from __future__ import annotations

from pathlib import Path
from typing import Any

from .arrays import split_and_check_arrays
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


def run_pipeline(args: Any) -> None:
    """Drive the pipeline. `args` must expose `.fasta`, `.output`,
    `.max_rep_size`, `.min_rep_size`, and (optionally) `.templates`.
    """
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    fasta_name = Path(args.fasta).name
    window_size = round((args.max_rep_size + KMER) * 1.1)

    # Load fasta + dedupe sequence headers (suffix any duplicates with 1, 2, ...)
    fasta = read_fasta_and_list(Path(args.fasta))
    if not fasta:
        return

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

    # Per-region: split into individual arrays + extract a representative.
    arr_rows = []
    for (name, seq), seq_regions, idx in zip(
        fasta, per_seq_regions, range(1, len(fasta) + 1)
    ):
        for region in seq_regions:
            region_seq = seq[region.start - 1:region.end]
            for arr in split_and_check_arrays(
                region_start=region.start,
                sequence=region_seq,
                seqID=name,
                numID=idx,
                max_repeat=args.max_rep_size,
                min_repeat=args.min_rep_size,
                kmer=KMER,
            ):
                arr_rows.append([
                    arr.start, arr.end, arr.seqID, arr.numID,
                    arr.score, arr.top_N, arr.top_5_N, arr.representative,
                ])

    write_csv_r_style(
        output_dir / f"{fasta_name}_aregarrays.csv",
        columns=["start", "end", "seqID", "numID", "score", "top_N", "top_5_N", "representative"],
        rows=arr_rows,
    )

    # Canonicalise each representative; if a templates fasta was supplied,
    # also classify each array against the templates.
    templates: list[tuple[str, str]] = []
    if getattr(args, "templates", None) is not None:
        templates = read_fasta_and_list(Path(args.templates))
        if not templates:
            raise SystemExit(f"No templates found in {args.templates}")
    template_names = {n for n, _ in templates}
    templates_by_name = dict(templates)

    stage7_rows = []
    for start, end, seqID, numID, score, top_N, top_5_N, representative in arr_rows:
        cls, shifted = shift_and_compare(representative, templates or None)
        stage7_rows.append({
            "start": start, "end": end, "seqID": seqID, "numID": numID,
            "score": score, "top_N": top_N, "top_5_N": top_5_N,
            "representative": shifted, "class": cls,
        })

    # Class assignment + per-class shift + sort + split off the unclassified.
    classarrays, no_repeats = classify_arrays(stage7_rows, template_names=template_names)

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
        return

    # Per-array repeat mapping with overlap/gap cleanup, rescore, and
    # edge-repeat refinement.
    fasta_by_seqID = {name: seq for name, seq in fasta}
    repeats_rows: list[dict] = []

    by_seq: dict[str, list[dict]] = {}
    for arr in classarrays:
        by_seq.setdefault(arr["seqID"], []).append(arr)

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

                rows = map_repeats(
                    representative=arr["representative"],
                    top_N=int(arr["top_N"]),
                    array_sequence=array_sequence,
                    seqID=seqID,
                    arrayID=int(arr["array_num_ID"]),
                    start_offset=a_start,
                    output_folder=output_dir,
                )
                if len(rows) < 2:
                    continue

                rows = clean_overlaps(
                    rows=rows,
                    representative_len=int(arr["top_N"]),
                    arr_class=arr["class"],
                )
                if len(rows) < 3:
                    continue

                rows = rescore_repeats(
                    rows=rows,
                    arr_representative=arr["representative"],
                    arr_class=arr["class"],
                    top_N=int(arr["top_N"]),
                    array_sequence=array_sequence,
                    array_start=a_start,
                    templates=templates_by_name,
                )

                rows = handle_edge_repeat(
                    rows=rows,
                    sequence_vector=sequence_substring,
                    sequence_vector_start=adjust_start,
                    template_sequence=templates_by_name.get(arr["class"], ""),
                )
                if len(rows) < 3:
                    continue

                for r in rows:
                    repeats_rows.append({k: r[k] for k in (*REPEATS_COLUMNS, "representative")})

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
