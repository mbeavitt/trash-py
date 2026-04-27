from __future__ import annotations

import random
from collections import Counter

from trash_py.arrays import ArrayBreaks, CollapsedKmer, _collapse_kmers, collapse_array_kmers


def _collapse_kmers_reference(
    sequence: str,
    array: ArrayBreaks,
    max_repeat: int,
    min_repeat: int,
    kmer: int = 10,
) -> list[CollapsedKmer]:
    global_min = 2
    max_edit = 2

    region_kmers = [sequence[i:i + kmer] for i in range(len(sequence) - kmer)]
    kmers_list_local = region_kmers[array.start - 1:array.end - kmer]
    if not kmers_list_local:
        return []

    counts = Counter(kmers_list_local)
    items = sorted(counts.items(), key=lambda x: (-x[1], x[0]))
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
        cluster_set = set(cluster_names)
        locations = [i + array.start for i, k in enumerate(kmers_list_local) if k in cluster_set]
        result.append(CollapsedKmer(
            kmers=cluster_names,
            count=cluster_count,
            locations=locations,
        ))
    return result


def _snapshot(collapsed: list[CollapsedKmer]) -> list[tuple[list[str], int, list[int]]]:
    return [
        (list(cluster.kmers), cluster.count, list(cluster.locations))
        for cluster in collapsed
    ]


def test_collapse_array_kmers_matches_reference_randomized() -> None:
    rng = random.Random(0)
    alphabet = "ACGTNn"

    for _ in range(30):
        sequence = "".join(rng.choice(alphabet) for _ in range(rng.randint(120, 260)))
        start = rng.randint(1, 20)
        end = rng.randint(max(start + 20, 40), len(sequence))
        array = ArrayBreaks(start=start, end=end, seqID="seq", numID=1)

        got = collapse_array_kmers(sequence, array, max_repeat=1000, min_repeat=7, kmer=10)
        expected = _collapse_kmers_reference(sequence, array, max_repeat=1000, min_repeat=7, kmer=10)

        assert _snapshot(got) == _snapshot(expected)
