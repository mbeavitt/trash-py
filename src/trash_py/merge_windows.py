"""Port of stage 05: `src/array/merge_windows.R`."""
from __future__ import annotations

from dataclasses import dataclass


@dataclass
class Region:
    start: int
    end: int
    score: float


def merge_windows(scores: list[float], window_size: int, sequence_full_length: int) -> list[Region]:
    threshold = 90

    if not any(s < threshold for s in scores):
        return []

    n = sequence_full_length
    starts = list(range(1, max(n - window_size + 1, 2), window_size))
    if len(starts) > 1 and n - starts[-1] < window_size / 2:
        starts.pop()
    # Each window's region extends 2 * window_size from its start, capped at n;
    # the final region always ends at n.
    ends = [min(s + 2 * window_size - 1, n) for s in starts]
    ends[-1] = n

    if len(starts) == 1:
        return [Region(starts[0], n, scores[0])]

    if not (len(starts) == len(ends) == len(scores)):
        raise RuntimeError(
            f"merge_windows length mismatch: starts={len(starts)} ends={len(ends)} scores={len(scores)}"
        )

    regions = [
        Region(s, e, sc)
        for s, e, sc in zip(starts, ends, scores)
        if sc < threshold
    ]
    if len(regions) < 2:
        return regions

    i = 0
    while i < len(regions) - 1:
        if (regions[i].end + 1) >= regions[i + 1].start:
            len_i = regions[i].end - regions[i].start
            len_j = regions[i + 1].end - regions[i + 1].start
            merged_score = (regions[i].score * len_i + regions[i + 1].score * len_j) / (len_i + len_j)
            regions[i] = Region(regions[i].start, regions[i + 1].end, merged_score)
            del regions[i + 1]
            i = max(i - 1, 0)
        else:
            i += 1
    return regions
