"""Port of stage 04: `src/array/sequence_window_score.R` and
`src/array/seq_win_score_int.R`."""
from __future__ import annotations

from collections import Counter


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
    """Produce a per-window singleton-kmer score across the whole sequence.

    Windows of `window_size` tile the sequence non-overlapping; the trailing
    window is dropped if less than half a window of sequence remains after it.
    """
    if (window_size // 2) <= kmer:
        raise ValueError("sequence_window_score: window size is too small")

    n = len(seq)
    starts = list(range(0, max(n - window_size, 1), window_size))
    if len(starts) > 1 and n - starts[-1] - 1 < window_size / 2:
        starts.pop()

    scores = [
        seq_win_score_int(1, window_size, kmer, seq[s : s + window_size])
        for s in starts
    ]

    if not scores or any(sc != sc for sc in scores):  # NaN check
        raise RuntimeError("sequence_window_score did not produce valid result")
    return scores
