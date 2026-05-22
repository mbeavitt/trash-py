"""Programmatic API for annotating individual sequences (reads) with TRASH.

:func:`annotate` and :func:`annotate_sequences` are higher-level wrappers
around the TRASH pipeline that operate directly on in-memory strings rather
than files, making it straightforward to embed repeat annotation in larger
Python workflows — for example, annotating individual insertion sequences
extracted from long-read alignments.

Example::

    from trash_py.reads import annotate

    result = annotate("ATCG..." * 12, name="read1")
    for unit in result.repeats:
        print(unit.start, unit.end, unit.width, unit.class_)

    for arr in result.arrays:
        print(arr.n_repeats, "×", arr.class_, "at", arr.start, "-", arr.end)
"""
from __future__ import annotations

import tempfile
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from . import _log as log
from .arrays import ArrayRow, split_and_check_arrays
from .classify import classify_arrays
from .mapping import (
    clean_overlaps,
    handle_edge_repeat,
    map_repeats,
    rescore_repeats,
)
from .merge_windows import merge_windows
from .pipeline import (
    ARRAYS_PER_CHUNK,
    KMER,
    _map_array_worker,
    _worker_init,
    _require_external_tools,
)
from .shift import shift_and_compare
from .summarise import (
    REPEATS_COLUMNS,
    strip_repeats_columns,
    summarise_arrays,
)
from .window_score import sequence_window_score


@dataclass
class RepeatUnit:
    """A single tandem repeat unit detected by nhmmer.

    Coordinates are 1-based and inclusive on both ends, matching the
    convention used in the TRASH CSV output.
    """

    start: int    #: 1-based start position within the input sequence.
    end: int      #: 1-based end position within the input sequence (inclusive).
    width: int    #: Span in bp (``end - start + 1``).
    strand: str   #: ``'+'`` or ``'-'``.
    array_id: int #: Identifier of the parent :class:`RepeatArray`.
    class_: str   #: TRASH class string, e.g. ``"178_1"``.
    score: float  #: nhmmer bit-score.


@dataclass
class RepeatArray:
    """A contiguous run of tandem repeat units.

    Coordinates are 1-based and inclusive on both ends.
    """

    start: int          #: 1-based start of the array.
    end: int            #: 1-based end of the array (inclusive).
    class_: str         #: TRASH class string, e.g. ``"178_1"``.
    n_repeats: int      #: Number of individual repeat units in this array.
    median_width: int   #: Median repeat unit width (bp).
    representative: str #: Representative monomer sequence.
    array_id: int       #: Unique numeric ID for this array.


@dataclass
class SequenceAnnotation:
    """TRASH annotation result for one input sequence.

    When no tandem structure is detected both *arrays* and *repeats* are
    empty lists.
    """

    name: str                              #: Sequence identifier.
    length: int                            #: Input sequence length (bp).
    arrays: list[RepeatArray] = field(default_factory=list)
    repeats: list[RepeatUnit] = field(default_factory=list)


def annotate(
    sequence: str,
    name: str = "seq",
    *,
    max_rep_size: int = 250,
    min_rep_size: int = 100,
    templates: Optional[list[tuple[str, str]]] = None,
    quiet: bool = True,
) -> SequenceAnnotation:
    """Annotate tandem repeats in a single DNA sequence.

    Parameters
    ----------
    sequence:
        DNA sequence string.  Any case is accepted; ambiguous bases are
        tolerated by the underlying tools.
    name:
        Identifier used for this sequence in the returned annotation.
    max_rep_size:
        Upper bound (bp) for the repeat unit size to search for.
    min_rep_size:
        Lower bound (bp) for the repeat unit size to search for.
    templates:
        Optional list of ``(name, sequence)`` pairs used to classify
        arrays.  When *None* classes are inferred from the unit length
        (e.g. ``"178_1"`` for a 178 bp unit).
    quiet:
        Suppress all trash-py log output (default *True*).

    Returns
    -------
    SequenceAnnotation
        Annotation result; *arrays* and *repeats* are empty when no
        tandem structure is found.

    Raises
    ------
    SystemExit
        If ``clustalo`` or ``nhmmer`` are not found on the PATH.
    """
    return annotate_sequences(
        [(name, sequence)],
        max_rep_size=max_rep_size,
        min_rep_size=min_rep_size,
        templates=templates,
        quiet=quiet,
    )[0]


def annotate_sequences(
    sequences: list[tuple[str, str]],
    *,
    max_rep_size: int = 250,
    min_rep_size: int = 100,
    templates: Optional[list[tuple[str, str]]] = None,
    processes: int = 1,
    quiet: bool = True,
) -> list[SequenceAnnotation]:
    """Annotate tandem repeats in multiple DNA sequences.

    Runs the full TRASH pipeline on all provided sequences and returns
    one :class:`SequenceAnnotation` per input, in input order.

    Parameters
    ----------
    sequences:
        List of ``(name, sequence)`` pairs.  Duplicate names are
        disambiguated automatically by appending ``1``, ``2``, … suffixes.
    max_rep_size:
        Upper bound (bp) for the repeat unit size to search for.
    min_rep_size:
        Lower bound (bp) for the repeat unit size to search for.
    templates:
        Optional list of ``(name, sequence)`` pairs for array
        classification.
    processes:
        Number of worker processes for parallel execution.  Values above 1
        are only beneficial when many sequences or arrays are processed at
        once.
    quiet:
        Suppress all trash-py log output (default *True*).

    Returns
    -------
    list[SequenceAnnotation]
        One result per input sequence, in the same order as *sequences*.

    Raises
    ------
    SystemExit
        If ``clustalo`` or ``nhmmer`` are not found on the PATH.
    """
    if quiet:
        log.configure(quiet=True)

    _require_external_tools()

    # Disambiguate duplicate sequence names (mirrors pipeline.run_pipeline)
    names = [n for n, _ in sequences]
    if len(names) != len(set(names)):
        counts: dict[str, int] = {}
        for n in names:
            counts[n] = counts.get(n, 0) + 1
        next_idx: dict[str, int] = {n: 1 for n, c in counts.items() if c > 1}
        deduped: list[tuple[str, str]] = []
        for name, seq in sequences:
            if counts[name] > 1:
                deduped.append((f"{name}{next_idx[name]}", seq))
                next_idx[name] += 1
            else:
                deduped.append((name, seq))
        sequences = deduped

    template_names = {n for n, _ in (templates or [])}
    templates_by_name: dict[str, str] = dict(templates or [])
    window_size = round((max_rep_size + KMER) * 1.1)

    # ── 1. Window scoring ──────────────────────────────────────────────────
    all_scores = [
        sequence_window_score(seq, window_size, KMER) for _, seq in sequences
    ]

    # ── 2. Merge windows → repetitive regions ─────────────────────────────
    per_seq_regions: list[list] = [
        list(merge_windows(scores, window_size, len(seq)))
        for (_, seq), scores in zip(sequences, all_scores)
    ]

    # ── 3. Split regions into candidate arrays ────────────────────────────
    arr_rows: list[list] = []
    for idx, ((name, seq), seq_regions) in enumerate(
        zip(sequences, per_seq_regions), start=1
    ):
        for region in seq_regions:
            region_seq = seq[region.start - 1:region.end]
            for arr in split_and_check_arrays(
                region_start=region.start,
                sequence=region_seq,
                seqID=name,
                numID=idx,
                max_repeat=max_rep_size,
                min_repeat=min_rep_size,
                kmer=KMER,
            ):
                arr_rows.append([
                    arr.start, arr.end, arr.seqID, arr.numID,
                    arr.score, arr.top_N, arr.top_5_N, arr.representative,
                ])

    if not arr_rows:
        return [
            SequenceAnnotation(name=n, length=len(s)) for n, s in sequences
        ]

    # ── 4. Canonicalise representatives + classify ────────────────────────
    stage7_rows = []
    for start, end, seqID, numID, score, top_N, top_5_N, representative in arr_rows:
        cls, shifted = shift_and_compare(representative, templates or None)
        stage7_rows.append({
            "start": start, "end": end, "seqID": seqID, "numID": numID,
            "score": score, "top_N": top_N, "top_5_N": top_5_N,
            "representative": shifted, "class": cls,
        })

    classarrays, _ = classify_arrays(
        stage7_rows, template_names=template_names
    )
    if not classarrays:
        return [
            SequenceAnnotation(name=n, length=len(s)) for n, s in sequences
        ]

    # ── 5. Map individual repeat units (nhmmer needs a temp dir) ──────────
    fasta_by_seqID = {name: seq for name, seq in sequences}
    repeats_rows: list[dict] = []

    with tempfile.TemporaryDirectory() as tmpdir:
        by_seq: dict[str, list[dict]] = {}
        for arr in classarrays:
            by_seq.setdefault(arr["seqID"], []).append(arr)

        map_tasks: list[tuple] = []
        for seqID, arrs_in_seq in by_seq.items():
            fasta_seq = fasta_by_seqID[seqID]
            n = len(arrs_in_seq)
            chunk_edges = list(range(0, n, ARRAYS_PER_CHUNK))
            chunk_edges.append(n - 1)
            for j in range(len(chunk_edges) - 1):
                lo, hi = chunk_edges[j], chunk_edges[j + 1]
                chunk = arrs_in_seq[lo:hi + 1]
                first_start = int(chunk[0]["start"])
                last_end = int(chunk[-1]["end"])
                sequence_substring = fasta_seq[first_start - 1:last_end - 1]
                adjust_start = first_start - 1
                for arr in chunk:
                    if arr["representative"] == "":
                        continue
                    a_start, a_end = int(arr["start"]), int(arr["end"])
                    array_sequence = fasta_seq[a_start - 1:a_end]
                    map_tasks.append((
                        arr, seqID, array_sequence, sequence_substring,
                        adjust_start, templates_by_name, tmpdir,
                    ))

        if processes > 1 and map_tasks:
            with ProcessPoolExecutor(
                max_workers=processes, initializer=_worker_init
            ) as ex:
                for rows, stats in ex.map(_map_array_worker, map_tasks):
                    log.merge_stats(stats)
                    repeats_rows.extend(rows)
        else:
            for task in map_tasks:
                repeats_rows.extend(_map_array_single(task))

    # ── 6. Aggregate ──────────────────────────────────────────────────────
    arrays_final = summarise_arrays(classarrays, repeats_rows)
    repeats_final = strip_repeats_columns(repeats_rows)

    # ── 7. Build per-sequence result objects ──────────────────────────────
    arrays_by_seq: dict[str, list[RepeatArray]] = {}
    for a in arrays_final:
        arrays_by_seq.setdefault(a["seqID"], []).append(RepeatArray(
            start=int(a["start"]),
            end=int(a["end"]),
            class_=a["class"],
            n_repeats=int(a["repeats_number"]),
            median_width=int(a["median_repeat_width"]),
            representative=a["representative"],
            array_id=int(a["array_num_ID"]),
        ))

    repeats_by_seq: dict[str, list[RepeatUnit]] = {}
    for r in repeats_final:
        repeats_by_seq.setdefault(r["seqID"], []).append(RepeatUnit(
            start=int(r["start"]),
            end=int(r["end"]),
            width=int(r["width"]),
            strand=r["strand"],
            array_id=int(r["arrayID"]),
            class_=r["class"],
            score=float(r["score"]),
        ))

    return [
        SequenceAnnotation(
            name=name,
            length=len(seq),
            arrays=arrays_by_seq.get(name, []),
            repeats=repeats_by_seq.get(name, []),
        )
        for name, seq in sequences
    ]


def _map_array_single(
    task: tuple[dict, str, str, str, int, dict, str],
) -> list[dict]:
    """Single-process variant of the array→repeat mapping step.

    Mirrors the logic of the private ``_map_array`` helper in
    ``pipeline.py`` but is local to this module so the public API does
    not depend on private pipeline internals.
    """
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
