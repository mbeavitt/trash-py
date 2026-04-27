from __future__ import annotations

import random

from trash_py import arrays


def test_window_comparison_scores_matches_python_reference_generic_case() -> None:
    rng = random.Random(0)
    alphabet = "ACGT"

    for _ in range(20):
        sequence = "".join(rng.choice(alphabet) for _ in range(160))
        kmer = 7
        n_kmers = len(sequence) - kmer

        window_starts = [rng.randint(-10, n_kmers + 10) for _ in range(12)]
        window_ends = [s + rng.randint(1, 25) for s in window_starts]
        window_starts_compare = [rng.randint(-10, n_kmers + 10) for _ in range(12)]
        window_ends_compare = [s + rng.randint(1, 25) for s in window_starts_compare]

        got = arrays._window_comparison_scores(
            sequence,
            kmer,
            window_starts,
            window_ends,
            window_starts_compare,
            window_ends_compare,
        )
        expected = arrays._window_comparison_scores_python(
            sequence,
            kmer,
            window_starts,
            window_ends,
            window_starts_compare,
            window_ends_compare,
        )
        assert got == expected


def test_window_comparison_scores_matches_python_reference_sliding_suffix_case() -> None:
    sequence = ("ACGTACGTAA" * 160) + ("TTTTCCCCGG" * 80)
    kmer = 10
    window_size = 100
    step = 5
    n_kmers = len(sequence) - kmer

    window_starts = list(range(1, 1 + 16 * step, step))
    window_ends = [s + window_size - 1 for s in window_starts]
    window_starts_compare = [e + 1 for e in window_ends]
    window_ends_compare = [n_kmers + 100] * len(window_starts)

    got = arrays._window_comparison_scores(
        sequence,
        kmer,
        window_starts,
        window_ends,
        window_starts_compare,
        window_ends_compare,
    )
    expected = arrays._window_comparison_scores_python(
        sequence,
        kmer,
        window_starts,
        window_ends,
        window_starts_compare,
        window_ends_compare,
    )
    assert got == expected


def test_chunk_a_split_arrays_fast_path_matches_fallback(monkeypatch) -> None:
    sequence = ("A" * 600) + ("CGTACGTAGC" * 220) + ("T" * 600)

    expected = arrays.chunk_a_split_arrays(sequence, seqID="seq1", numID=7)

    monkeypatch.setattr(arrays, "_window_compare_scores_fast", None)
    got = arrays.chunk_a_split_arrays(sequence, seqID="seq1", numID=7)

    assert got == expected
