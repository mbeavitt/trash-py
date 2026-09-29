"""The C classification loop must reproduce the Python reference exactly.

`_classify_loop` is the original implementation, kept as the fallback for
non-ASCII representatives; `_classify_native` is the C port. Cases cover
non-ACGT characters, representatives shorter than the kmer, zero scores and
near-identical classes where fw/rv ties decide the strand.
"""
from __future__ import annotations

import copy
import random

import pytest

from trash_py import classify as C


def _random_rows(rng: random.Random, n: int) -> list[dict]:
    alpha = rng.choice(["ac", "acgt", "acgtn", "aacgtN-", "at"])
    bases = [
        "".join(rng.choice(alpha) for _ in range(rng.randint(1, 40)))
        for _ in range(rng.randint(1, 6))
    ]
    rows = []
    for _ in range(n):
        rep = list(rng.choice(bases) * rng.randint(1, 3))
        for _ in range(rng.randint(0, 3)):
            rep[rng.randrange(len(rep))] = rng.choice(alpha)
        cut = rng.randrange(len(rep))
        if rng.random() < 0.3:
            rep = rep[cut:] + rep[:cut]
        if rng.random() < 0.3:
            rep = list(C.rev_comp_string("".join(rep)))
        start = rng.randint(0, 1000)
        # width >= 1 and score >= 0 keep every importance positive, so the
        # reference loop always terminates.
        width = rng.choice([1, 5, 50, 500, rng.randint(1, 5000)])
        score = rng.choice([0.0, 1.0, 12.5, rng.random() * 100])
        rows.append({
            "start": start, "end": start + width, "score": score,
            "representative": "".join(rep)[: rng.randint(1, 60)], "class": "",
        })
    return rows


@pytest.mark.parametrize("seed", range(20))
def test_native_matches_reference(seed: int) -> None:
    rng = random.Random(seed)
    for _ in range(50):
        rows = _random_rows(rng, rng.randint(1, 120))
        ref, nat = copy.deepcopy(rows), copy.deepcopy(rows)
        C._classify_loop(ref)
        C._classify_native(nat)
        assert [(r["class"], r["representative"]) for r in nat] == \
               [(r["class"], r["representative"]) for r in ref]


def test_no_progress_raises() -> None:
    # Row 1 has importance 0 and no width-compatible top, so the reference
    # loop re-picks the already-classified row 0 forever. Native fails loudly.
    rows = [
        {"start": 0, "end": 100, "score": 1.0, "representative": "acgtacgtac", "class": ""},
        {"start": 0, "end": 0, "score": 1.0, "representative": "a" * 40, "class": ""},
    ]
    with pytest.raises(RuntimeError, match="no progress"):
        C._classify_native(rows)


def _first_match_scan(pattern: str, kmers: list[str], start_idx: int) -> int:
    """The original linear `first_match` that compare_kmer_grep replaced."""
    for j in range(start_idx, len(kmers)):
        if kmers[j] == pattern:
            return j - start_idx + 1
    return 0


def _compare_kmer_grep_scan(sequence_kmers, seq, max_size_dif, string_length, kmer):
    n = len(seq)
    lo = C.math.floor(string_length * (1 - max_size_dif))
    hi = C.math.ceil(string_length * (1 + max_size_dif))
    if not (lo <= n <= hi) or n == 0:
        return seq
    copies_base = -(-(n + kmer) // n)
    copies = 1 + copies_base
    n_kmers = n * copies_base
    ext_fw = seq * copies
    ext_rv = C.rev_comp_string(seq) * copies
    kfw = [ext_fw[j:j + kmer] for j in range(n_kmers)]
    krv = [ext_rv[j:j + kmer] for j in range(n_kmers)]
    dfw = [_first_match_scan(sequence_kmers[i], kfw, i) for i in range(len(sequence_kmers))]
    drv = [_first_match_scan(sequence_kmers[i], krv, i) for i in range(len(sequence_kmers))]
    if sum(dfw) + sum(drv) == 0:
        return seq
    if sum(1 for d in dfw if d) > sum(1 for d in drv if d):
        shift, ext = C._mode_smallest([d for d in dfw if d]), ext_fw
    else:
        shift, ext = C._mode_smallest([d for d in drv if d]), ext_rv
    return ext[shift - 1:shift - 1 + n]


@pytest.mark.parametrize("seed", range(10))
def test_compare_kmer_grep_matches_scan(seed: int) -> None:
    rng = random.Random(seed)
    for _ in range(200):
        rows = _random_rows(rng, 2)
        cls, target = rows[0]["representative"], rows[1]["representative"]
        kmers = C._circular_kmers(cls, C.KMER_SHIFT)
        assert C.compare_kmer_grep(kmers, target, 1, len(cls)) == \
               _compare_kmer_grep_scan(kmers, target, 1, len(cls), C.KMER_SHIFT)
