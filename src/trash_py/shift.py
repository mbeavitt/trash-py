"""Per-array representative canonicalisation and template matching.

`shift_sequence` rotates (and optionally rev-complements) a representative
to a deterministic orientation. `shift_and_compare` then compares against
an optional template fasta — the best-scoring template (under the 0.4
adist threshold) becomes the array's class.
"""
from __future__ import annotations

import math
from typing import Sequence

from rapidfuzz.distance import Levenshtein

from .sequence import rev_comp_string


def kmer_hash_score(kmer: str) -> int:
    """Position-weighted base-4 hash of a kmer (`a=0, c=1, t=2, g=3`)."""
    values = {"a": 0, "c": 1, "t": 2, "g": 3}
    total = 0
    for i, ch in enumerate(kmer, start=1):
        total += (4 ** i) * values.get(ch, 0)
    return total


def _pearson_corr(x: list[int], y: list[int]) -> float:
    n = len(x)
    mx = sum(x) / n
    my = sum(y) / n
    sxy = 0.0
    sxx = 0.0
    syy = 0.0
    for xi, yi in zip(x, y):
        dx = xi - mx
        dy = yi - my
        sxy += dx * dy
        sxx += dx * dx
        syy += dy * dy
    if sxx == 0 or syy == 0:
        return math.nan
    return sxy / math.sqrt(sxx * syy)


def _argmin_first(values: list[float]) -> int:
    best_idx = -1
    best_val = math.inf
    for i, v in enumerate(values):
        if math.isnan(v):
            continue
        if v < best_val:
            best_val = v
            best_idx = i
    return best_idx


def shift_sequence(sequence: str, k: int = 6) -> str:
    if len(set(sequence)) <= 1:
        return sequence

    n = len(sequence)
    repeats_needed = -(-(n + k) // n)

    def shift_scores(seq: str) -> list[float]:
        extended = seq * repeats_needed
        kmers = [extended[i:i + k] for i in range(n)]
        scores = [kmer_hash_score(km) for km in kmers]
        scores_doubled = scores + scores
        x = list(range(1, n + 1))
        return [_pearson_corr(x, scores_doubled[j:j + n]) for j in range(n)]

    fw = shift_scores(sequence)
    rev = shift_scores(rev_comp_string(sequence))

    all_shifts = fw + rev
    best_idx = _argmin_first(all_shifts)
    shift = best_idx + 1  # 1-based

    if shift > n:
        sequence = rev_comp_string(sequence)
        shift -= n

    if shift != 1:
        sequence = sequence[shift - 1:] + sequence[:shift - 1]
    return sequence


def compare_circular(
    sequence: str, template: str, max_size_dif: float = 0.15
) -> tuple[float, str]:
    """Best-scoring rotation (or rev-comp rotation) of `sequence` against
    `template`. Returns `(adist / max(len), best_rotation)`. Skips work
    and returns `(1.0, sequence)` if lengths differ by more than `max_size_dif`."""
    n_seq, n_tpl = len(sequence), len(template)
    lo = math.floor(n_tpl * (1 - max_size_dif))
    hi = math.ceil(n_tpl * (1 + max_size_dif))
    if not (lo <= n_seq <= hi):
        return 1.0, sequence

    rotations: list[str] = [sequence[x:] + sequence[:x] for x in range(1, n_seq + 1)]
    rotations.extend([rev_comp_string(r) for r in rotations])

    denom = max(n_seq, n_tpl)
    best_score = math.inf
    best_idx = 0
    for i, r in enumerate(rotations):
        d = Levenshtein.distance(template, r) / denom
        if d < best_score:
            best_score = d
            best_idx = i
    return best_score, rotations[best_idx]


def shift_and_compare(
    sequence: str, templates: Sequence[tuple[str, str]] | None = None
) -> tuple[str, str]:
    """Canonicalise `sequence` and (optionally) classify against templates.

    `templates` is `[(name, seq), ...]`. Without templates, returns
    `("", shifted)`. With templates, returns `(best_name, best_rotation)`
    when the best score is ≤ 0.4, else `("", shifted)`.
    """
    if sequence == "":
        return "", ""

    shifted = shift_sequence(sequence)
    if not templates:
        return "", shifted

    score_threshold = 0.4
    max_size_dif = 0.15
    best_score = math.inf
    best_name = ""
    best_seq = shifted
    for name, tseq in templates:
        score, shifted_for_template = compare_circular(shifted, tseq, max_size_dif)
        if score < best_score:
            best_score = score
            best_name = name
            best_seq = shifted_for_template

    if best_score <= score_threshold:
        return best_name, best_seq
    return "", shifted
