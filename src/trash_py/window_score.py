"""Port of stage 04: `src/array/sequence_window_score.R` and
`src/array/seq_win_score_int.R`."""
from __future__ import annotations

from ._ext import seq_win_score_int


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
