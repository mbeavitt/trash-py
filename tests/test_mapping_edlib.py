"""Unit tests for the edlib-based repeat mapper (map_edlib).

Tests use synthetic tandem arrays with known ground-truth positions so we
can verify correctness independently of nhmmer.
"""
from __future__ import annotations

import random

import pytest

from trash_py.mapping import map_edlib
from trash_py.sequence import rev_comp_string


def _make_unit(length: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    return "".join(rng.choice("acgt") for _ in range(length))


def _mutate(seq: str, rate: float, seed: int = 0) -> str:
    rng = random.Random(seed)
    result = []
    for base in seq:
        if rng.random() < rate:
            result.append(rng.choice([b for b in "acgt" if b != base]))
        else:
            result.append(base)
    return "".join(result)


def _insert_random(seq: str, n_inserts: int, seed: int = 0) -> str:
    rng = random.Random(seed)
    s = list(seq)
    for _ in range(n_inserts):
        pos = rng.randint(0, len(s))
        s.insert(pos, rng.choice("acgt"))
    return "".join(s)


# ---------------------------------------------------------------------------
# Basic detection
# ---------------------------------------------------------------------------

def test_finds_exact_forward_repeats():
    unit = _make_unit(50)
    rows = map_edlib(unit, unit * 10, "chr1", 1, 1)
    fwd = [r for r in rows if r["strand"] == "+"]
    assert len(fwd) == 10


def test_finds_exact_reverse_repeats():
    unit = _make_unit(50)
    rev = rev_comp_string(unit)
    rows = map_edlib(unit, rev * 10, "chr1", 1, 1)
    rev_hits = [r for r in rows if r["strand"] == "-"]
    assert len(rev_hits) == 10


def test_no_false_positives_on_random_sequence():
    unit = _make_unit(50, seed=1)
    # Completely random sequence with no relationship to unit
    rng = random.Random(99)
    noise = "".join(rng.choice("acgt") for _ in range(500))
    rows = map_edlib(unit, noise, "chr1", 1, 1)
    assert len(rows) == 0


# ---------------------------------------------------------------------------
# Coordinate correctness
# ---------------------------------------------------------------------------

def test_forward_coordinates_exact():
    unit = _make_unit(50)
    array = unit * 3
    rows = sorted(
        [r for r in map_edlib(unit, array, "chr1", 1, 1) if r["strand"] == "+"],
        key=lambda r: r["start"],
    )
    assert len(rows) == 3
    for i, r in enumerate(rows):
        assert r["start"] == i * 50 + 1, f"unit {i}: start {r['start']} != {i*50+1}"
        assert r["end"] == (i + 1) * 50, f"unit {i}: end {r['end']} != {(i+1)*50}"


def test_start_offset_shifts_coordinates():
    unit = _make_unit(50)
    array = unit * 3
    rows_base = sorted(map_edlib(unit, array, "chr1", 1, 1), key=lambda r: r["start"])
    rows_off = sorted(map_edlib(unit, array, "chr1", 1, 1000), key=lambda r: r["start"])
    assert len(rows_base) == len(rows_off)
    for r0, r1 in zip(rows_base, rows_off):
        assert r1["start"] == r0["start"] + 999  # +1000 offset, -1 for 0→1 base
        assert r1["end"] == r0["end"] + 999


def test_start_offset_default_1_is_1based():
    unit = _make_unit(50)
    rows = map_edlib(unit, unit, "chr1", 1, 1)
    fwd = [r for r in rows if r["strand"] == "+"]
    assert len(fwd) == 1
    assert fwd[0]["start"] == 1
    assert fwd[0]["end"] == 50


# ---------------------------------------------------------------------------
# Metadata fields
# ---------------------------------------------------------------------------

def test_seqid_and_arrayid_propagated():
    unit = _make_unit(30)
    rows = map_edlib(unit, unit * 5, "my_contig", 42, 1)
    for r in rows:
        assert r["seqID"] == "my_contig"
        assert r["arrayID"] == 42


def test_eval_is_sentinel():
    unit = _make_unit(50)
    rows = map_edlib(unit, unit * 5, "chr1", 1, 1)
    for r in rows:
        assert r["eval"] == -1.0


def test_strand_values():
    unit = _make_unit(50)
    rev = rev_comp_string(unit)
    rows = map_edlib(unit, unit * 3 + rev * 3, "chr1", 1, 1)
    strands = {r["strand"] for r in rows}
    assert strands <= {"+", "-"}


# ---------------------------------------------------------------------------
# Tolerance: substitutions and indels
# ---------------------------------------------------------------------------

def test_tolerates_15pct_substitutions():
    unit = _make_unit(177)  # CEN178-like length
    units = [_mutate(unit, 0.15, seed=i) for i in range(20)]
    rows = map_edlib(unit, "".join(units), "chr1", 1, 1)
    fwd = [r for r in rows if r["strand"] == "+"]
    assert len(fwd) >= 18, f"found {len(fwd)}/20"


def test_tolerates_single_indels():
    unit = _make_unit(50)
    rng = random.Random(7)
    units = []
    for i in range(10):
        if rng.random() < 0.4:
            units.append(_insert_random(unit, 1, seed=i))
        else:
            units.append(unit)
    rows = map_edlib(unit, "".join(units), "chr1", 1, 1)
    fwd = [r for r in rows if r["strand"] == "+"]
    assert len(fwd) >= 8, f"found {len(fwd)}/10"


def test_mixed_strand_array():
    unit = _make_unit(50)
    rev = rev_comp_string(unit)
    array = unit * 5 + rev * 5
    rows = map_edlib(unit, array, "chr1", 1, 1)
    fwd = [r for r in rows if r["strand"] == "+"]
    rev_hits = [r for r in rows if r["strand"] == "-"]
    assert len(fwd) >= 4, f"found {len(fwd)} fwd"
    assert len(rev_hits) >= 4, f"found {len(rev_hits)} rev"


# ---------------------------------------------------------------------------
# Edge cases
# ---------------------------------------------------------------------------

def test_empty_array_returns_empty():
    unit = _make_unit(50)
    rows = map_edlib(unit, "", "chr1", 1, 1)
    assert rows == []


def test_array_shorter_than_unit():
    unit = _make_unit(50)
    rows = map_edlib(unit, unit[:20], "chr1", 1, 1)
    assert rows == []


def test_single_repeat():
    unit = _make_unit(50)
    rows = map_edlib(unit, unit, "chr1", 1, 1)
    fwd = [r for r in rows if r["strand"] == "+"]
    assert len(fwd) == 1
