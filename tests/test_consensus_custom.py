"""Unit tests for the edlib-seeded consensus that replaces clustalo when
``--fast`` is set.

The ``_consensus_custom`` function in ``trash_py.arrays`` is meant to
match the behaviour of ``_clustalo_align`` + ``_consensus_N`` on the
near-identical inputs the pipeline actually feeds it.
"""
from __future__ import annotations

import shutil

import pytest

from trash_py.arrays import (
    _clustalo_align,
    _consensus_custom,
    _consensus_N,
    _nw_project_to_seed,
)


def _majority_consensus(seqs: list[str]) -> str:
    """Reference (test-only) consensus for equal-length inputs: per-column
    majority vote with the same (g,c,t,a) tie-break as the production
    code."""
    bases = ("g", "c", "t", "a")
    L = len(seqs[0])
    out = []
    for c in range(L):
        col = [s[c] for s in seqs]
        counts = [col.count(b) for b in bases]
        best = 0
        for k in range(1, 4):
            if counts[k] > counts[best]:
                best = k
        out.append(bases[best])
    return "".join(out)


# ---- equal-length fast path -------------------------------------------------

def test_equal_length_majority_vote_full_length():
    seqs = [
        "acgtacgtacgt",
        "acgtacgtacgt",
        "acgtacgtacgt",
        "acgtacgtacgg",  # one SNP at last col
    ]
    got = _consensus_custom(seqs, len(seqs[0]))
    assert got == _majority_consensus(seqs)


def test_equal_length_top_N_picks_first_columns_on_tie():
    # Equal occupancy across all 8 cols → tie-break is column index.
    seqs = [
        "acgtacgt",
        "acgtacgt",
        "acgtacgt",
    ]
    # top_N=4 should pick the first 4 cols (ascending index on tie).
    assert _consensus_custom(seqs, 4) == "acgt"


def test_uppercase_input_is_normalised():
    # The production pipeline passes lowercase, but be defensive.
    seqs = ["ACGT", "ACGT", "ACGT"]
    assert _consensus_custom(seqs, 4) == "acgt"


def test_empty_input_returns_empty():
    assert _consensus_custom([], 10) == ""
    assert _consensus_custom(["acgt"], 0) == ""


def test_single_input_uses_itself_as_consensus():
    # Caller in arrays.py:431 short-circuits len==1, but defensively the
    # function still returns sensible output.
    assert _consensus_custom(["acgt"], 4) == "acgt"


# ---- variable-length (median-seed projection) path -------------------------

def test_variable_length_drops_minority_long_input():
    # Four 5-bp + one 6-bp. Median length is 5 → seed is the first 5-bp
    # input. The 6-bp query gets its extra base treated as an insertion
    # vs seed and dropped, so the consensus is just the majority of the
    # 5-bp inputs.
    seqs = ["acgta", "acgta", "acgta", "acgta", "acgtag"]
    assert _consensus_custom(seqs, 5) == "acgta"


def test_variable_length_short_outlier_falls_below_occupancy():
    # 4 length-5 + 1 length-3. Median = 5. The length-3 input projects
    # as "acg--" — cols 3 and 4 get only 4 votes vs the majority 5.
    # top_N=5 keeps all 5 cols since they're the only ones available.
    seqs = ["acgta", "acgta", "acgta", "acgta", "acg"]
    assert _consensus_custom(seqs, 5) == "acgta"


# ---- seed-projection mode --------------------------------------------------

def test_seed_mode_handles_query_insertion():
    # Query has an extra base in the middle relative to the seed;
    # NW projection should drop the inserted base.
    seed = "acgtac"
    assert _nw_project_to_seed("acgxxtac".replace("xx", "tg"), seed) == "acgtac"


def test_seed_mode_handles_query_deletion():
    # Query is missing the middle base relative to the seed; NW projection
    # places a gap at that seed column.
    seed = "acgtac"
    out = _nw_project_to_seed("acgac", seed)
    assert len(out) == len(seed)
    assert out.count("-") == 1


def test_seed_mode_consensus_ignores_short_outlier():
    # 9 near-identical 6-bp samples + 1 length-3 partial outlier. With
    # left-anchor padding the outlier would skew positions 3-5; with
    # seed-projection the outlier becomes "acg---" which only weakens
    # occupancy for cols 3-5 (still 9 votes from the others).
    seed = "acgtac"
    sample = "acgtac"
    samples = [sample] * 9 + ["acg"]
    assert _consensus_custom(samples, 6, seed=seed) == "acgtac"


def test_seed_mode_consensus_corrects_seed_snp():
    # Seed has the "wrong" base at col 2 (g vs c); 8 of 8 samples vote
    # for "c". The new consensus should reflect the majority, not the seed.
    seed = "acgtac"
    samples = ["acctac"] * 8
    out = _consensus_custom(samples, 6, seed=seed)
    assert out == "acctac"


# ---- parity vs clustalo (only runs if clustalo is installed) ---------------

@pytest.mark.skipif(
    shutil.which("clustalo") is None, reason="clustalo not on PATH"
)
def test_matches_clustalo_within_5pct_on_realistic_repeats():
    # Fabricate a small set of ~178bp near-identical CEN178-like
    # monomers, then check that the custom consensus is within the
    # 5% Hamming budget the parity suite already tolerates.
    base = (
        "agtcagtgacgaagaggcttgtttagagcacgaaagcgactagcaattgaaagtagccacc"
        "tttccgaaaagtgaaagactaacgtactcggttcagtactgctaaaactttctgttgaaca"
        "atgtcaagtatctgaacgtactagcatactcggtttgctcgtacgg"
    )
    seqs = [base] * 8
    # Sprinkle a couple of SNPs into half of them.
    seqs[1] = base[:30] + "t" + base[31:]
    seqs[3] = base[:90] + "g" + base[91:]
    seqs[5] = base[:150] + "a" + base[151:]

    top_N = len(base)
    custom = _consensus_custom(seqs, top_N)
    clu = _consensus_N(_clustalo_align(seqs, "clustalo"), top_N)

    assert len(custom) == len(clu) == top_N
    hamming = sum(a != b for a, b in zip(custom, clu))
    assert hamming / top_N <= 0.05, (
        f"custom vs clustalo Hamming {hamming}/{top_N} exceeds 5% budget"
    )
