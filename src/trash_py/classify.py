"""Group arrays into classes and align each class to a single orientation.

`classify_arrays` is the public entry point: greedy kmer-Jaccard
classification of each array's representative, followed by a per-class
shift so every member is rotated/oriented to match the most-important
representative in its class.
"""
from __future__ import annotations

import math
from collections import Counter
from typing import Any

from .sequence import rev_comp_string


KMER_CLASSIFY = 7
SIZE_DIF_TO_CHECK = 0.15
MAX_DISTANCE_TO_CLASSIFY = 0.85
KMER_SHIFT = 6


def _circular_kmers(sequence: str, k: int) -> list[str]:
    n = len(sequence)
    repeats_needed = -(-(n + k) // n)
    extended = sequence * repeats_needed
    return [extended[i:i + k] for i in range(n)]


def _coverage_count(source: list[str], target_set: set[str]) -> int:
    """Count kmers in `source` that appear anywhere in `target_set`.

    Mirrors R `sum(source %in% target)` — elements at each source position
    are checked against the target membership, so duplicates in source count
    each occurrence.
    """
    return sum(1 for km in source if km in target_set)


def _argmax_first(values: list[float]) -> int:
    """Mirror R `which.max`: 0-based index of first maximum (NaN skipped)."""
    best_idx = -1
    best_val = -math.inf
    for i, v in enumerate(values):
        if v > best_val:
            best_val = v
            best_idx = i
    return best_idx


def classify_repeats(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Greedy classify-by-kmer-Jaccard in descending importance order.

    `rows` keys must include start, end, score, representative, class (class
    may be empty). Returns a shallow copy of each row with class populated
    and representative rev-comped when the best match is on the rv strand.
    """
    out: list[dict[str, Any]] = [dict(r) for r in rows]

    # Stage 7 leaves representative="" for arrays that didn't resolve; those
    # are tagged here rather than going through classify.
    for row in out:
        if row["representative"] == "":
            row["class"] = "none_identified"

    which_to_classify = [i for i, r in enumerate(out) if r["class"] == ""]
    if not which_to_classify:
        return out

    temp = [out[i] for i in which_to_classify]
    n = len(temp)

    width = [int(r["end"]) - int(r["start"]) for r in temp]
    score = [float(r["score"]) for r in temp]
    importance: list[float] = [w * s if s != 0 else float(w) for w, s in zip(width, score)]
    rep_width = [len(r["representative"]) for r in temp]
    class_vec: list[str] = [r["class"] for r in temp]

    kmers_fw: list[list[str]] = []
    kmers_rv: list[list[str]] = []
    for r in temp:
        rep = r["representative"]
        kmers_fw.append(_circular_kmers(rep, KMER_CLASSIFY))
        kmers_rv.append(_circular_kmers(rev_comp_string(rep), KMER_CLASSIFY))

    names_iterator = 1
    while any(c == "" for c in class_vec):
        # `which.max` over importance — ties break on first occurrence.
        # Already-classified rows have importance 0 (set below), so they
        # can only be picked when every remaining row is also 0; they're
        # then filtered out by the class=="" mask on `which_to_compare`.
        which_top = _argmax_first(importance)
        top_rep_width = rep_width[which_top]
        new_class_name = f"{top_rep_width}_{names_iterator}"

        lo = math.floor(top_rep_width * (1 - SIZE_DIF_TO_CHECK))
        hi = math.ceil(top_rep_width * (1 + SIZE_DIF_TO_CHECK))
        which_to_compare = [
            j for j in range(n)
            if lo <= rep_width[j] <= hi and class_vec[j] == "" and j != which_top
        ]

        if not which_to_compare:
            class_vec[which_top] = new_class_name
            importance[which_top] = 0.0
            names_iterator += 1
            continue

        top_fw_list = kmers_fw[which_top]
        top_len = len(top_fw_list)

        # scores_fw/rv ask "how many of top's fw kmers appear in candidate's
        # fw (or rv) kmers?" — the comparison direction is top→candidate.
        scores_fw = [_coverage_count(top_fw_list, set(kmers_fw[j])) for j in which_to_compare]
        scores_rv = [_coverage_count(top_fw_list, set(kmers_rv[j])) for j in which_to_compare]

        distances = [1 - s / top_len for s in scores_fw]
        distances_rv = [1 - s / top_len for s in scores_rv]

        # Per-candidate: keep the stronger direction; disable the other.
        for i in range(len(distances)):
            if distances[i] >= distances_rv[i]:
                distances[i] = 1.0
            else:
                distances_rv[i] = 1.0

        similar = [i for i, d in enumerate(distances) if d <= MAX_DISTANCE_TO_CLASSIFY]
        similar_rv = [i for i, d in enumerate(distances_rv) if d <= MAX_DISTANCE_TO_CLASSIFY]

        for i in similar:
            j = which_to_compare[i]
            class_vec[j] = new_class_name
            importance[j] = 0.0
        for i in similar_rv:
            j = which_to_compare[i]
            class_vec[j] = new_class_name
            importance[j] = 0.0
            temp[j]["representative"] = rev_comp_string(temp[j]["representative"])

        class_vec[which_top] = new_class_name
        importance[which_top] = 0.0
        names_iterator += 1

    for idx, orig_idx in enumerate(which_to_classify):
        out[orig_idx]["class"] = class_vec[idx]
        out[orig_idx]["representative"] = temp[idx]["representative"]

    return out


def _mode_smallest(values: list[int]) -> int:
    """R `as.numeric(names(table(x))[which.max(table(x))])`.

    `table()` sorts by value ascending, so ties on count break on the
    smallest value.
    """
    counts = Counter(values)
    return min(counts, key=lambda k: (-counts[k], k))


def compare_kmer_grep(
    sequence_kmers: list[str],
    sequence_to_realign: str,
    max_size_dif: float,
    string_length: int,
    kmer: int = KMER_SHIFT,
) -> str:
    """Shift `sequence_to_realign` so its kmers align with `sequence_kmers`.

    Picks the modal offset at which class kmers first appear in the target's
    circular kmer stream (fw vs rv), then returns the target rotated to
    that offset. Returns the input unchanged if the length gate rejects
    it or no kmers match.
    """
    realign_length = len(sequence_to_realign)
    lo = math.floor(string_length * (1 - max_size_dif))
    hi = math.ceil(string_length * (1 + max_size_dif))
    if not (lo <= realign_length <= hi):
        return sequence_to_realign
    if realign_length == 0:
        return sequence_to_realign

    copies_base = -(-(realign_length + kmer) // realign_length)  # ceil((L+k)/L)
    copies = 1 + copies_base
    n_kmers = realign_length * copies_base

    extended_fw = sequence_to_realign * copies
    kmers_fw = [extended_fw[j:j + kmer] for j in range(n_kmers)]

    rev = rev_comp_string(sequence_to_realign)
    extended_rv = rev * copies
    kmers_rv = [extended_rv[j:j + kmer] for j in range(n_kmers)]

    def first_match(pattern: str, kmers: list[str], start_idx: int) -> int:
        """1-based index of first equal kmer in kmers[start_idx:], or 0."""
        for j in range(start_idx, len(kmers)):
            if kmers[j] == pattern:
                return j - start_idx + 1
        return 0

    distances_fw = [first_match(sequence_kmers[i], kmers_fw, i) for i in range(len(sequence_kmers))]
    distances_rv = [first_match(sequence_kmers[i], kmers_rv, i) for i in range(len(sequence_kmers))]

    if sum(distances_fw) + sum(distances_rv) == 0:
        return sequence_to_realign

    fw_nonzero = sum(1 for d in distances_fw if d != 0)
    rv_nonzero = sum(1 for d in distances_rv if d != 0)

    if fw_nonzero > rv_nonzero:
        nonzero = [d for d in distances_fw if d != 0]
        shift = _mode_smallest(nonzero)
        extended = extended_fw
    else:
        nonzero = [d for d in distances_rv if d != 0]
        shift = _mode_smallest(nonzero)
        extended = extended_rv

    return extended[shift - 1:shift - 1 + realign_length]


def shift_classes(rows: list[dict[str, Any]], kmer: int = KMER_SHIFT) -> list[str]:
    """Return shifted representatives for all rows in one class.

    Anchors on the row with highest importance (width * score, falling back
    to width when score == 0), builds its circular kmer stream, and aligns
    every row's representative to it via compare_kmer_grep.
    """
    if not rows:
        return []

    widths = [int(r["end"]) - int(r["start"]) for r in rows]
    scores = [float(r["score"]) for r in rows]
    importance: list[float] = [w * s if s != 0 else float(w) for w, s in zip(widths, scores)]

    top_idx = _argmax_first(importance)
    class_sequence = rows[top_idx]["representative"]
    string_length = len(class_sequence)
    if string_length == 0:
        return [r["representative"] for r in rows]

    class_kmers = _circular_kmers(class_sequence, kmer)

    return [
        compare_kmer_grep(class_kmers, r["representative"], 1, string_length, kmer=kmer)
        for r in rows
    ]


def classify_arrays(
    rows: list[dict[str, Any]],
    template_names: set[str] | None = None,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Classify + shift + sort + split.

    Returns `(classarrays, no_repeats_arrays)`. Each row is a fresh dict
    with `array_num_ID` populated; sort order matches R: by seqID, then
    by start (both stable).

    `template_names` are class names that came from a templates fasta —
    those skip `shift_classes` (their representatives were already aligned
    to the template in stage 07) and rbind to the tail like none_identified.
    """
    classified = classify_repeats(rows)
    template_names = template_names or set()

    # R uses `unique()` which preserves first-occurrence order.
    seen: dict[str, None] = {}
    for r in classified:
        c = r["class"]
        if c != "none_identified" and c not in template_names and c not in seen:
            seen[c] = None
    class_order = list(seen.keys())

    for cls in class_order:
        idx = [i for i, r in enumerate(classified) if r["class"] == cls]
        subset = [classified[i] for i in idx]
        shifted = shift_classes(subset)
        for k, new_rep in zip(idx, shifted):
            classified[k]["representative"] = new_rep

    # R: rbind(processed classes in order, then template/none_identified rows)
    processed: list[dict[str, Any]] = []
    for cls in class_order:
        processed.extend(r for r in classified if r["class"] == cls)
    tail_rows = [
        r for r in classified
        if r["class"] in template_names or r["class"] == "none_identified"
    ]
    ordered = processed + tail_rows

    ordered.sort(key=lambda r: int(r["start"]))
    ordered.sort(key=lambda r: r["seqID"])
    for i, r in enumerate(ordered, start=1):
        r["array_num_ID"] = i

    classarrays = [r for r in ordered if r["class"] != "none_identified"]
    no_repeats = [r for r in ordered if r["class"] == "none_identified"]
    return classarrays, no_repeats
