from __future__ import annotations

import random
from collections import Counter

from trash_py.arrays import ArrayBreaks, CollapsedKmer, find_top_repeat_distances


def _clone_collapsed(collapsed: list[CollapsedKmer]) -> list[CollapsedKmer]:
    return [
        CollapsedKmer(
            kmers=list(cluster.kmers),
            count=cluster.count,
            locations=list(cluster.locations),
            distances=list(cluster.distances),
        )
        for cluster in collapsed
    ]


def _snapshot(collapsed: list[CollapsedKmer]) -> list[tuple[list[int], list[int]]]:
    return [(list(cluster.locations), list(cluster.distances)) for cluster in collapsed]


def _find_top_repeat_distances_reference(
    collapsed: list[CollapsedKmer],
    array: ArrayBreaks,
    max_repeat: int,
    min_repeat: int,
) -> tuple[int, str, list[int]]:
    small_window_size = 1000
    small_window_step = 100
    small_window_min = small_window_size / small_window_step

    kmer_starts: list[int] = []
    distances: list[int] = []
    for cluster in collapsed:
        locs = cluster.locations
        if len(locs) < 2:
            cluster.locations = []
            cluster.distances = []
            continue
        diffs = [locs[k + 1] - locs[k] for k in range(len(locs) - 1)]
        starts = locs[:-1]
        paired = [(s, d) for s, d in zip(starts, diffs) if min_repeat <= d <= max_repeat]
        cluster.locations = [p[0] for p in paired]
        cluster.distances = [p[1] for p in paired]
        kmer_starts.extend(cluster.locations)
        distances.extend(cluster.distances)

    kmer_starts_2 = [s + d for s, d in zip(kmer_starts, distances)]

    window_starts = list(
        range(array.start, array.end - small_window_step + 1, small_window_step)
    ) or [array.start]
    if len(window_starts) > 1 and array.end - window_starts[-1] < small_window_step / 2:
        window_starts.pop()
    if len(window_starts) < 2:
        window_ends = [array.end]
    else:
        window_ends = [window_starts[k + 1] - 1 for k in range(len(window_starts) - 1)] + [array.end]
    if len(window_ends) != len(window_starts):
        window_ends = [array.end]
    if len(window_ends) > 1 and (window_ends[-1] - window_starts[-1]) < (small_window_step / 2):
        window_ends[-2] = window_ends[-1]
        window_ends.pop()
        window_starts.pop()

    extend = small_window_size - small_window_step
    window_ends = [min(e + extend, array.end) for e in window_ends]

    rem_starts = list(kmer_starts)
    rem_starts_2 = list(kmer_starts_2)
    rem_distances = list(distances)

    moving_top_distance = [0] * len(window_starts)
    for j, (ws, we) in enumerate(zip(window_starts, window_ends)):
        mask = [
            (ws <= s <= we) or (ws <= s2 <= we)
            for s, s2 in zip(rem_starts, rem_starts_2)
        ]
        if sum(mask) > small_window_min:
            in_distances = [d for d, m in zip(rem_distances, mask) if m]
            counts = Counter(in_distances)
            max_count = max(counts.values())
            moving_top_distance[j] = min(d for d, c in counts.items() if c == max_count)
        rem_starts = [s for s, m in zip(rem_starts, mask) if not m]
        rem_starts_2 = [s for s, m in zip(rem_starts_2, mask) if not m]
        rem_distances = [d for d, m in zip(rem_distances, mask) if not m]

    nonzero = [d for d in moving_top_distance if d != 0]
    if not nonzero:
        return 0, "", []

    counts_by_value = Counter(nonzero)
    sorted_by_name = sorted(counts_by_value.items(), key=lambda x: x[0])

    merged: list[tuple[str, int]] = [(str(sorted_by_name[0][0]), sorted_by_name[0][1])]
    last_name = sorted_by_name[0][0]
    for name_val, cnt in sorted_by_name[1:]:
        if name_val == last_name + 1:
            prev_name, prev_count = merged[-1]
            merged[-1] = (f"{prev_name},{name_val}", prev_count + cnt)
            last_name = name_val
        else:
            merged.append((str(name_val), cnt))
            last_name = name_val

    merged.sort(key=lambda x: -x[1])

    count_ns = min(5, len(merged))
    top_5_n = "".join(f"N_{merged[j][0]}_Count_{merged[j][1]}." for j in range(count_ns))

    first_name = merged[0][0]
    top_n_distances = [int(x) for x in first_name.split(",")]
    top_n = int(sum(top_n_distances) / len(top_n_distances))

    return top_n, top_5_n, top_n_distances


def test_find_top_repeat_distances_matches_reference_randomized() -> None:
    rng = random.Random(0)

    for _ in range(40):
        start = rng.randint(1, 200)
        end = start + rng.randint(1200, 5000)
        array = ArrayBreaks(start=start, end=end, seqID="seq", numID=1)

        collapsed: list[CollapsedKmer] = []
        for cluster_idx in range(rng.randint(1, 16)):
            locations = sorted(
                rng.sample(range(start, end + 1), rng.randint(0, 60))
            )
            collapsed.append(
                CollapsedKmer(
                    kmers=[f"k{cluster_idx}"],
                    count=len(locations),
                    locations=locations,
                    distances=[],
                )
            )

        got_input = _clone_collapsed(collapsed)
        ref_input = _clone_collapsed(collapsed)

        got = find_top_repeat_distances(got_input, array, max_repeat=1000, min_repeat=7)
        expected = _find_top_repeat_distances_reference(ref_input, array, max_repeat=1000, min_repeat=7)

        assert got == expected
        assert _snapshot(got_input) == _snapshot(ref_input)


def test_find_top_repeat_distances_matches_reference_overlap_case() -> None:
    array = ArrayBreaks(start=1, end=1600, seqID="seq", numID=1)
    collapsed = [
        CollapsedKmer(kmers=["a"], count=14, locations=[50, 240, 430, 620, 810, 1000, 1190], distances=[]),
        CollapsedKmer(kmers=["b"], count=14, locations=[80, 270, 460, 650, 840, 1030, 1220], distances=[]),
        CollapsedKmer(kmers=["c"], count=12, locations=[110, 300, 490, 680, 870, 1060], distances=[]),
    ]

    got_input = _clone_collapsed(collapsed)
    ref_input = _clone_collapsed(collapsed)

    got = find_top_repeat_distances(got_input, array, max_repeat=1000, min_repeat=7)
    expected = _find_top_repeat_distances_reference(ref_input, array, max_repeat=1000, min_repeat=7)

    assert got == expected
    assert _snapshot(got_input) == _snapshot(ref_input)
