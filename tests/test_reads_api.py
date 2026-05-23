"""Tests for the programmatic reads annotation API (trash_py.reads)."""
from __future__ import annotations

import pytest

from trash_py.reads import (
    RepeatArray,
    RepeatUnit,
    SequenceAnnotation,
    annotate,
    annotate_sequences,
)


# ── dataclass smoke tests ─────────────────────────────────────────────────────

def test_repeat_unit_fields():
    u = RepeatUnit(start=1, end=177, width=177, strand="+", array_id=1,
                   class_="177_1", score=45.2)
    assert u.end - u.start + 1 == u.width


def test_repeat_array_fields():
    a = RepeatArray(start=1, end=531, class_="177_1", n_repeats=3,
                    median_width=177, representative="ACGT", array_id=1)
    assert a.n_repeats == 3


def test_sequence_annotation_defaults():
    ann = SequenceAnnotation(name="read1", length=500)
    assert ann.arrays == []
    assert ann.repeats == []


# ── length guard ──────────────────────────────────────────────────────────────

def test_annotate_sequences_raises_on_short_input():
    """Sequences shorter than window_size must raise ValueError with a clear
    message rather than crashing inside window scoring."""
    short_seq = "ACGT" * 10  # 40 bp — well below the ~286 bp minimum
    with pytest.raises(ValueError, match="shorter than the"):
        annotate_sequences([("short", short_seq)])


def test_annotate_raises_on_short_input():
    short_seq = "ACGT" * 10
    with pytest.raises(ValueError, match="shorter than the"):
        annotate(short_seq, name="short")


# ── public API: non-repetitive sequence returns empty annotation ──────────────

def test_annotate_no_repeats():
    """A long non-repetitive sequence should return an empty annotation."""
    # 600 bp, no tandem structure — passes window scoring
    seq = "ACGTTTGCAA" * 30 + "TTGCAACGTA" * 30
    result = annotate(seq, name="nonrep")
    assert isinstance(result, SequenceAnnotation)
    assert result.name == "nonrep"
    assert result.length == len(seq)
    assert result.arrays == []
    assert result.repeats == []


def test_annotate_sequences_preserves_order():
    seqs = [
        ("a", "ACGT" * 150),
        ("b", "TTGC" * 150),
        ("c", "GCAA" * 150),
    ]
    results = annotate_sequences(seqs)
    assert [r.name for r in results] == ["a", "b", "c"]
    assert all(r.length == 600 for r in results)


def test_annotate_sequences_deduplicates_names():
    seqs = [("dup", "ACGT" * 150), ("dup", "TTGC" * 150)]
    results = annotate_sequences(seqs)
    assert len(results) == 2
    assert all(isinstance(r, SequenceAnnotation) for r in results)


# ── top-level package re-exports ─────────────────────────────────────────────

def test_package_exports():
    import trash_py

    assert hasattr(trash_py, "annotate")
    assert hasattr(trash_py, "annotate_sequences")
    assert hasattr(trash_py, "SequenceAnnotation")
    assert hasattr(trash_py, "RepeatArray")
    assert hasattr(trash_py, "RepeatUnit")
    assert hasattr(trash_py, "ARABIDOPSIS_CEN178")
    name, seq = trash_py.ARABIDOPSIS_CEN178
    assert name == "CEN178"
    assert len(seq) == 177
