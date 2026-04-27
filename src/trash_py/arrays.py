"""Identify individual repeat arrays inside a repetitive region.

`split_and_check_arrays` is the public entry point. Internally it
runs four chunks that mirror the upstream R source:

* `chunk_a_split_arrays` — detect array-break positions.
* `chunk_b_collapse_kmers` — per-array kmer counting, filtering, and
  Levenshtein-based clustering.
* `chunk_c_top_n` — distance analysis yielding `top_N` and the
  comma-joined `top_5_N` string.
* `chunk_d_consensus` — MSA + majority-base consensus to produce the
  array representative and finalise `top_N`.
"""
from __future__ import annotations

import subprocess
import tempfile
from bisect import bisect_left
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

from .window_score import seq_win_score_int

try:
    from ._ext import window_compare_scores as _window_compare_scores_fast
except ImportError:
    _window_compare_scores_fast = None


@dataclass
class ArrayBreaks:
    """Initial per-array record — boundaries only, before any kmer analysis.
    Coordinates are 1-based, inclusive, **window-relative** (start=1 at
    the first base of the region we were handed)."""
    start: int
    end: int
    seqID: str
    numID: int


@dataclass
class ArrayRow:
    """One row of the final `_aregarrays.csv` table."""
    start: int
    end: int
    seqID: str
    numID: int
    score: float
    top_N: int
    top_5_N: str
    representative: str


@dataclass
class CollapsedKmer:
    """One kmer cluster: its member strings, total occurrence count, and
    1-based absolute-window-relative anchor positions of every
    occurrence."""
    kmers: list[str]
    count: int
    locations: list[int] = field(default_factory=list)
    distances: list[int] = field(default_factory=list)


def chunk_a_split_arrays(
    sequence: str,
    seqID: str,
    numID: int,
    max_repeat: int = 1000,
    kmer: int = 10,
) -> list[ArrayBreaks]:
    window_step = -(-max_repeat // 20)  # ceil(max_repeat / 20)
    min_detach = 0.08
    min_split = 0.04
    array_overlaps = 0

    start = 1
    end = len(sequence)

    window_size = max(max_repeat, 200)
    if window_size > end:
        window_size = end

    window_starts: list[int] = []
    window_ends: list[int] = []
    window_starts_compare: list[int] = []
    windows_comparison_score: list[float] = []

    if (end - window_size) > start:
        window_starts = list(
            range(start, end - window_size - window_step + 1, window_step)
        ) or [start]
        if len(window_starts) > 1 and (end - window_size) - window_starts[-1] < window_step / 2:
            window_starts.pop()
        if len(window_starts) < 2:
            window_ends = [end - window_size]
        else:
            window_ends = [s + window_size - 1 for s in window_starts]

        cap = end - window_size / 2
        keep = [i for i, e in enumerate(window_ends) if e <= cap]
        window_starts = [window_starts[i] for i in keep]
        window_ends = [window_ends[i] for i in keep]

        window_starts_compare = [e + 1 for e in window_ends]
        window_ends_compare = [s + window_size - 1 for s in window_starts_compare]
        window_ends_compare = [end if ec <= end else ec for ec in window_ends_compare]

        windows_comparison_score = _window_comparison_scores(
            sequence,
            kmer,
            window_starts,
            window_ends,
            window_starts_compare,
            window_ends_compare,
        )

    single = [ArrayBreaks(start=start, end=end, seqID=seqID, numID=numID)]

    if not windows_comparison_score:
        return single
    if not any(s > min_detach for s in windows_comparison_score):
        return single

    window_event = ["new_top"]
    window_event_position = [1]
    first_above = next(
        j for j, s in enumerate(windows_comparison_score, start=1) if s > min_detach
    )
    i = first_above + 1
    n = len(window_starts_compare)
    while i <= n:
        while i <= n and window_event[-1] == "new_top":
            if windows_comparison_score[i - 1] <= min_split:
                window_event.append("new_bottom")
                window_event_position.append(i)
            i += 1
        while i <= n and window_event[-1] == "new_bottom":
            if windows_comparison_score[i - 1] >= min_detach:
                window_event.append("new_top")
                window_event_position.append(i)
            i += 1

    if len(window_event) <= 2:
        return single

    array_breaks: list[int] = []
    for i in range(len(window_event) - 1):
        if window_event[i] == "new_bottom" and window_event[i + 1] == "new_top":
            s_idx = window_event_position[i]
            e_idx = window_event_position[i + 1]
            subset = windows_comparison_score[s_idx - 1:e_idx]
            argmin = min(range(len(subset)), key=lambda k: subset[k])
            array_breaks.append(s_idx + argmin)

    break_coords = [window_ends[b - 1] + 1 for b in array_breaks]

    arrays = [ArrayBreaks(
        start=1,
        end=break_coords[0] - 1 + array_overlaps,
        seqID=seqID,
        numID=numID,
    )]
    for i in range(len(array_breaks) - 1):
        arrays.append(ArrayBreaks(
            start=break_coords[i] - array_overlaps,
            end=break_coords[i + 1] - 1 + array_overlaps,
            seqID=seqID,
            numID=numID,
        ))
    arrays.append(ArrayBreaks(
        start=break_coords[-1],
        end=end,
        seqID=seqID,
        numID=numID,
    ))
    return arrays


def _slice_1based(lst: list, start: int, end: int) -> list:
    """Mimic R's `lst[start:end]` with 1-based inclusive indices; out-of-range
    positions become None (matching R's NA behaviour that filters to FALSE
    in `%in%` membership)."""
    return [lst[i - 1] if 1 <= i <= len(lst) else None for i in range(start, end + 1)]


def _window_comparison_scores(
    sequence: str,
    kmer: int,
    window_starts: list[int],
    window_ends: list[int],
    window_starts_compare: list[int],
    window_ends_compare: list[int],
) -> list[float]:
    if not window_starts_compare:
        return []
    if _window_compare_scores_fast is not None:
        return _window_compare_scores_fast(
            sequence,
            kmer,
            window_starts,
            window_ends,
            window_starts_compare,
            window_ends_compare,
        )
    return _window_comparison_scores_python(
        sequence,
        kmer,
        window_starts,
        window_ends,
        window_starts_compare,
        window_ends_compare,
    )


def _window_comparison_scores_python(
    sequence: str,
    kmer: int,
    window_starts: list[int],
    window_ends: list[int],
    window_starts_compare: list[int],
    window_ends_compare: list[int],
) -> list[float]:
    kmers_list = [sequence[i:i + kmer] for i in range(len(sequence) - kmer)]
    scores = [0.0] * len(window_starts_compare)
    for i in range(len(window_starts_compare)):
        a = _slice_1based(kmers_list, window_starts[i], window_ends[i])
        b = _slice_1based(kmers_list, window_starts_compare[i], window_ends_compare[i])
        a_set = {x for x in a if x is not None}
        scores[i] = sum(1 for x in b if x in a_set) / len(b)
    return scores


def chunk_b_collapse_kmers(
    sequence: str,
    array: ArrayBreaks,
    max_repeat: int,
    min_repeat: int,
    kmer: int = 10,
) -> list[CollapsedKmer]:
    global_min = 2
    max_edit = 2

    kmers_list_local = [
        sequence[i:i + kmer]
        for i in range(array.start - 1, array.end - kmer)
    ]
    if not kmers_list_local:
        return []

    positions_by_kmer: dict[str, list[int]] = {}
    for pos, kmer_name in enumerate(kmers_list_local, start=array.start):
        positions_by_kmer.setdefault(kmer_name, []).append(pos)

    items = sorted(
        ((name, len(locations)) for name, locations in positions_by_kmer.items()),
        key=lambda x: (-x[1], x[0]),
    )
    items = [(n, c) for n, c in items if "N" not in n and "n" not in n]
    if not items:
        return []
    items = [(n, c) for n, c in items if c >= global_min]
    if not items:
        return []

    min_kmers_count = (len(items) // 1000) + 1
    if min_kmers_count > global_min:
        filtered = [(n, c) for n, c in items if c >= min_kmers_count]
        if not filtered:
            top_n = -(-len(items) // 10)
            items = items[:top_n]
        elif len(filtered) > 5000:
            items = items[:5000]
        else:
            items = filtered

    clusters = _collapse_kmers([n for n, _ in items], [c for _, c in items], max_edit)

    result: list[CollapsedKmer] = []
    for cluster_names, cluster_count in clusters:
        locations: list[int] = []
        for name in cluster_names:
            locations.extend(positions_by_kmer.get(name, ()))
        locations.sort()
        result.append(CollapsedKmer(
            kmers=cluster_names,
            count=cluster_count,
            locations=locations,
        ))
    return result


def _collapse_kmers(names: list[str], counts: list[int], max_edit: int) -> list[tuple[list[str], int]]:
    """Greedy cluster by Hamming distance ≤ max_edit on same-length kmers.

    Anchors are chosen in ascending-count order; ties broken alphabetically.
    """
    pairs = sorted(zip(names, counts), key=lambda p: (p[1], p[0]))
    names = [p[0] for p in pairs]
    counts = [p[1] for p in pairs]

    clusters: list[tuple[list[str], int]] = []
    i = 0
    while i < len(names):
        cluster_names = [names[i]]
        cluster_count = counts[i]
        if i + 1 < len(names):
            anchor = names[i]
            to_merge = [j for j in range(i + 1, len(names))
                        if _hamming_le(anchor, names[j], max_edit)]
            for j in to_merge:
                cluster_count += counts[j]
                cluster_names.append(names[j])
            for j in reversed(to_merge):
                del names[j]
                del counts[j]
        clusters.append((cluster_names, cluster_count))
        i += 1
    return clusters


def _hamming_le(a: str, b: str, k: int) -> bool:
    if len(a) != len(b):
        return False
    d = 0
    for x, y in zip(a, b):
        if x != y:
            d += 1
            if d > k:
                return False
    return True


def chunk_c_top_n(
    collapsed: list[CollapsedKmer],
    array: ArrayBreaks,
    max_repeat: int,
    min_repeat: int,
) -> tuple[int, str, list[int]]:
    """Return `(top_N, top_5_N_string, top_N_distances)`."""
    small_window_size = 1000
    small_window_step = 100
    small_window_min = small_window_size / small_window_step  # 10

    window_starts = list(
        range(array.start, array.end - small_window_step + 1, small_window_step)
    ) or [array.start]
    if len(window_starts) > 1 and array.end - window_starts[-1] < small_window_step / 2:
        window_starts.pop()
    if len(window_starts) < 2:
        window_ends = [array.end]
    else:
        window_ends = [window_starts[k + 1] - 1 for k in range(len(window_starts) - 1)] + [array.end]
    if len(window_ends) != len(window_starts):
        window_ends = [array.end]
    if len(window_ends) > 1 and (window_ends[-1] - window_starts[-1]) < (small_window_step / 2):
        window_ends[-2] = window_ends[-1]
        window_ends.pop()
        window_starts.pop()

    extend = small_window_size - small_window_step
    window_ends = [min(e + extend, array.end) for e in window_ends]

    window_distance_counts: list[dict[int, int] | None] = [None] * len(window_starts)
    window_distance_totals = [0] * len(window_starts)

    for cluster in collapsed:
        locs = cluster.locations
        if len(locs) < 2:
            cluster.locations = []
            cluster.distances = []
            continue

        filtered_starts: list[int] = []
        filtered_distances: list[int] = []
        prev = locs[0]
        for curr in locs[1:]:
            distance = curr - prev
            if min_repeat <= distance <= max_repeat:
                filtered_starts.append(prev)
                filtered_distances.append(distance)

                first_window = _first_matching_window(prev, window_starts, window_ends)
                second_window = _first_matching_window(prev + distance, window_starts, window_ends)
                if first_window == -1:
                    window_idx = second_window
                elif second_window == -1 or first_window < second_window:
                    window_idx = first_window
                else:
                    window_idx = second_window

                if window_idx != -1:
                    counts = window_distance_counts[window_idx]
                    if counts is None:
                        counts = {}
                        window_distance_counts[window_idx] = counts
                    counts[distance] = counts.get(distance, 0) + 1
                    window_distance_totals[window_idx] += 1
            prev = curr

        cluster.locations = filtered_starts
        cluster.distances = filtered_distances

    moving_top_distance = [0] * len(window_starts)
    counts_by_value: dict[int, int] = {}
    for j, total in enumerate(window_distance_totals):
        if total <= small_window_min:
            continue

        counts = window_distance_counts[j]
        if counts is None:
            continue

        best_distance = 0
        best_count = -1
        for distance, count in counts.items():
            if count > best_count or (count == best_count and distance < best_distance):
                best_distance = distance
                best_count = count

        moving_top_distance[j] = best_distance
        counts_by_value[best_distance] = counts_by_value.get(best_distance, 0) + 1

    if not counts_by_value:
        return 0, "", []

    sorted_by_name = sorted(counts_by_value.items(), key=lambda x: x[0])

    merged: list[tuple[str, int]] = [(str(sorted_by_name[0][0]), sorted_by_name[0][1])]
    last_name = sorted_by_name[0][0]
    for name_val, cnt in sorted_by_name[1:]:
        if name_val == last_name + 1:
            prev_name, prev_count = merged[-1]
            merged[-1] = (f"{prev_name},{name_val}", prev_count + cnt)
            last_name = name_val
        else:
            merged.append((str(name_val), cnt))
            last_name = name_val

    merged.sort(key=lambda x: -x[1])

    count_ns = min(5, len(merged))
    top_5_n = "".join(f"N_{merged[j][0]}_Count_{merged[j][1]}." for j in range(count_ns))

    first_name = merged[0][0]
    top_n_distances = [int(x) for x in first_name.split(",")]
    top_n = int(sum(top_n_distances) / len(top_n_distances))  # floor mean (positive values)

    return top_n, top_5_n, top_n_distances


def _first_matching_window(point: int, window_starts: list[int], window_ends: list[int]) -> int:
    idx = bisect_left(window_ends, point)
    if idx < len(window_starts) and window_starts[idx] <= point:
        return idx
    return -1


def chunk_d_consensus(
    sequence: str,
    collapsed: list[CollapsedKmer],
    top_N_distances: list[int],
    array: ArrayBreaks,
    clustalo_exe: str = "clustalo",
) -> tuple[float, int, str]:
    """Return `(score, top_N_final, representative)`."""
    max_repeats_to_align = 10
    sub = sequence[array.start - 1:array.end]
    win = seq_win_score_int(1, array.end - array.start + 1, 10, sub)
    score = 100.0 - win

    if not collapsed or not top_N_distances:
        return score, 0, ""

    top_N_set = set(top_N_distances)
    topN_counts = [sum(1 for d in c.distances if d in top_N_set) for c in collapsed]

    best_idx = 0
    best_count = topN_counts[0]
    for k in range(1, len(topN_counts)):
        if topN_counts[k] > best_count:
            best_idx = k
            best_count = topN_counts[k]

    top_kmer = collapsed[best_idx]
    paired = [(loc, d) for loc, d in zip(top_kmer.locations, top_kmer.distances) if d in top_N_set]
    if not paired:
        return score, 0, ""

    locations = [p[0] for p in paired]
    distances = [p[1] for p in paired]

    sorted_d = sorted(distances)
    n = len(sorted_d)
    median = sorted_d[n // 2] if n % 2 else (sorted_d[n // 2 - 1] + sorted_d[n // 2]) / 2
    top_N = int(median)

    top_kmer_list = [sequence[loc - 1:loc - 1 + dist] for loc, dist in zip(locations, distances)]
    if not top_kmer_list:
        return score, top_N, ""

    if len(top_kmer_list) > max_repeats_to_align:
        indices = _evenly_spaced_indices(len(top_kmer_list), max_repeats_to_align)
        top_kmer_list = [top_kmer_list[i] for i in indices]

    if len(top_kmer_list) == 1:
        consensus = top_kmer_list[0]
    else:
        alignment = _clustalo_align(top_kmer_list, clustalo_exe)
        consensus = _consensus_N(alignment, top_N)

    return score, top_N, consensus


def _evenly_spaced_indices(n: int, k: int) -> list[int]:
    """Mirror `round(seq(1, n, length.out = k))` then convert to 0-based."""
    if k == 1:
        return [0]
    step = (n - 1) / (k - 1)
    return [round(1 + i * step) - 1 for i in range(k)]


def _clustalo_align(sequences: list[str], clustalo_exe: str) -> list[str]:
    with tempfile.TemporaryDirectory() as tmp:
        in_path = Path(tmp) / "in.fa"
        out_path = Path(tmp) / "out.fa"
        with in_path.open("w") as f:
            for i, seq in enumerate(sequences, start=1):
                f.write(f">{i}\n{seq}\n")
        subprocess.run(
            [
                clustalo_exe,
                "-i", str(in_path),
                "-o", str(out_path),
                "--force",
                "--threads=1",
                "--output-order=input-order",
            ],
            check=True,
            capture_output=True,
        )
        aligned: list[str] = []
        chunks: list[str] = []
        with out_path.open() as f:
            for line in f:
                if line.startswith(">"):
                    if chunks:
                        aligned.append("".join(chunks))
                    chunks = []
                else:
                    chunks.append(line.strip())
            if chunks:
                aligned.append("".join(chunks))
    return [s.lower() for s in aligned]


def split_and_check_arrays(
    region_start: int,
    sequence: str,
    seqID: str,
    numID: int,
    max_repeat: int = 1000,
    min_repeat: int = 7,
    kmer: int = 10,
    clustalo_exe: str = "clustalo",
) -> list[ArrayRow]:
    """One call = one repetitive region → one or more arrays, each with
    all columns of the final `_aregarrays.csv` row populated.

    `region_start` is the 1-based fasta-absolute start of the region;
    window-relative coordinates from chunks A-D are shifted back to
    fasta-absolute before returning.
    """
    arrays = chunk_a_split_arrays(sequence, seqID, numID, max_repeat, kmer)
    rows: list[ArrayRow] = []
    for arr in arrays:
        collapsed = chunk_b_collapse_kmers(sequence, arr, max_repeat, min_repeat, kmer)
        score: float | None = None
        top_N, top_5_N, representative = 0, "", ""
        if collapsed:
            _, top_5_N, top_N_distances = chunk_c_top_n(collapsed, arr, max_repeat, min_repeat)
            if top_N_distances:
                score, top_N, representative = chunk_d_consensus(
                    sequence, collapsed, top_N_distances, arr, clustalo_exe=clustalo_exe,
                )
        if score is None:
            sub = sequence[arr.start - 1:arr.end]
            win = seq_win_score_int(1, arr.end - arr.start + 1, kmer, sub)
            score = 100.0 - win

        rows.append(ArrayRow(
            start=arr.start + region_start - 1,
            end=arr.end + region_start - 1,
            seqID=seqID,
            numID=numID,
            score=score,
            top_N=top_N,
            top_5_N=top_5_N,
            representative=representative,
        ))
    return rows


def _consensus_N(alignment: list[str], N: int) -> str:
    if not alignment or N <= 0:
        return ""
    ncol = len(alignment[0])
    frequencies = [sum(1 for s in alignment if s[c] != "-") for c in range(ncol)]
    order = sorted(range(ncol), key=lambda c: (-frequencies[c], c))
    top_cols = set(order[:N])

    bases = ("g", "c", "t", "a")
    result: list[str] = []
    for c in range(ncol):
        if c not in top_cols:
            continue
        col = [s[c] for s in alignment]
        counts = [col.count(b) for b in bases]
        best = 0
        best_c = counts[0]
        for k in range(1, 4):
            if counts[k] > best_c:
                best = k
                best_c = counts[k]
        result.append(bases[best])
    return "".join(result)
