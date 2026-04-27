"""Port of stage 04: `src/array/sequence_window_score.R` and
`src/array/seq_win_score_int.R`."""
from __future__ import annotations

from collections import Counter

from .genomic_bins import genomic_bins_starts


def seq_win_score_int(start: int, end: int, kmer: int, seq: str) -> float:
    """Score a window by the proportion of singleton kmers.

    Lower score = more repetition. Returns 100 when the window is too
    short or too thinly sampled to be informative.

    `seq` is a 0-based substring spanning the window. `start` and `end`
    are 1-based bounds on which kmer start positions inside that
    substring are considered.
    """
    if (end - start) <= kmer:
        return 100.0

    # R: X in start:(end-kmer), kmer = seq[X:X+kmer-1] (1-based, inclusive).
    # Python 0-based equivalent: i in range(start-1, end-kmer).
    counts: Counter[str] = Counter()
    for i in range(start - 1, end - kmer):
        counts[seq[i:i + kmer]] += 1

    for key in [k for k in counts if "n" in k or "N" in k]:
        del counts[key]

    total = sum(counts.values())
    if total < kmer * 2:
        return 100.0

    singletons = sum(v for v in counts.values() if v == 1)
    return 100.0 * singletons / total


def sequence_window_score(seq: str, window_size: int, kmer: int = 10) -> list[float]:
    """Produce a per-window singleton-kmer score across the whole sequence."""
    if (window_size // 2) <= kmer:
        raise ValueError("sequence_window_score: window size is too small")

    n = len(seq)
    sliding = window_size  # R: ceiling(window_size / 1)

    if sliding >= n:
        starts = [1]
        ends = [n]
    else:
        starts = genomic_bins_starts(start=1, end=n, bin_size=sliding)
        starts = [s for s in starts if s < n]
        if len(starts) == 1:
            ends = [n]
        else:
            ends = [starts[i + 1] - 1 + window_size for i in range(len(starts) - 1)] + [n + window_size]
            ends = [min(e, n) for e in ends]

    scores: list[float] = []
    for s, e in zip(starts, ends):
        sub = seq[s - 1:e]  # 1-based inclusive → 0-based slice
        scores.append(seq_win_score_int(1, window_size, kmer, sub))

    if not scores or any(sc != sc for sc in scores):  # NaN check
        raise RuntimeError("sequence_window_score did not produce valid result")
    return scores
