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


# ── dataclass smoke tests (no external tools needed) ─────────────────────────

def test_repeat_unit_fields():
    u = RepeatUnit(start=1, end=178, width=178, strand="+", array_id=1,
                   class_="178_1", score=45.2)
    assert u.end - u.start + 1 == u.width


def test_repeat_array_fields():
    a = RepeatArray(start=1, end=534, class_="178_1", n_repeats=3,
                    median_width=178, representative="ACGT", array_id=1)
    assert a.n_repeats == 3


def test_sequence_annotation_defaults():
    ann = SequenceAnnotation(name="read1", length=500)
    assert ann.arrays == []
    assert ann.repeats == []


# ── public API: short sequence returns empty annotation ──────────────────────

def test_annotate_no_repeats(monkeypatch):
    """A non-repetitive sequence should return an empty annotation without
    calling any external tools (clustalo/nhmmer absent on CI is fine).
    We monkeypatch _require_external_tools so the test is self-contained.
    """
    import trash_py.reads as reads_mod

    monkeypatch.setattr(reads_mod, "_require_external_tools", lambda: None)

    # 600 bp random-ish, non-repetitive sequence — passes window scoring
    seq = ("ACGTTTGCAA" * 30 + "TTGCAACGTA" * 30)  # 600 bp, no tandem structure
    result = annotate(seq, name="nonrep")
    assert isinstance(result, SequenceAnnotation)
    assert result.name == "nonrep"
    assert result.length == len(seq)
    assert result.arrays == []
    assert result.repeats == []


def test_annotate_sequences_preserves_order(monkeypatch):
    import trash_py.reads as reads_mod

    monkeypatch.setattr(reads_mod, "_require_external_tools", lambda: None)

    # 600 bp each, no tandem structure
    seqs = [
        ("a", "ACGT" * 150),
        ("b", "TTGC" * 150),
        ("c", "GCAA" * 150),
    ]
    results = annotate_sequences(seqs)
    assert [r.name for r in results] == ["a", "b", "c"]
    assert all(r.length == 600 for r in results)


def test_annotate_sequences_deduplicates_names(monkeypatch):
    import trash_py.reads as reads_mod

    monkeypatch.setattr(reads_mod, "_require_external_tools", lambda: None)

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
