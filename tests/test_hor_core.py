"""Unit tests for the native HOR scanner (`_ext.find_hors`).

Every expected value here was produced by running the reference `HOR.V3.3`
binary on the same tiny alignment, so these lock the Python/C port to the
original tool's behaviour — including its quirks (see docs/HOR_source_bugs.md).
"""
from __future__ import annotations

from pathlib import Path

from trash_py import _ext


def _write_aln(path: Path, records: list[tuple[str, str]]) -> str:
    with path.open("w") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")
    return str(path)


def test_method1_parallel_hor(tmp_path: Path) -> None:
    # ABC ABC — a tandem (parallel) HOR of 3 units, exact match.
    aln = _write_aln(tmp_path / "a.fasta", [
        ("1_D1", "AAAA"), ("2_D1", "CCCC"), ("3_D1", "GGGG"),
        ("4_D1", "AAAA"), ("5_D1", "CCCC"), ("6_D1", "GGGG"),
    ])
    assert _ext.find_hors(aln, 1, 0, 3, 1) == [(1, 3, 4, 6, 1, 0)]


def test_method1_antiparallel_hor(tmp_path: Path) -> None:
    # ABC / reverse-CBA with opposite strands — a perpendicular (inverted) HOR.
    aln = _write_aln(tmp_path / "b.fasta", [
        ("1_D1", "AAAA"), ("2_D1", "CCCC"), ("3_D1", "GGGG"),
        ("4_D2", "GGGG"), ("5_D2", "CCCC"), ("6_D2", "AAAA"),
    ])
    assert _ext.find_hors(aln, 1, 0, 3, 1) == [(1, 3, 4, 6, 2, 0)]


def test_edge_close_undercounts_last_pair(tmp_path: Path) -> None:
    # Three matching pairs (each 1 SNV) whose run reaches the array boundary.
    # The original closes the HOR at the edge using the SNV total captured
    # *before* the final pair — so total_variant is 2, not 3 (upstream bug).
    aln = _write_aln(tmp_path / "c.fasta", [
        ("1_D1", "AAAA"), ("2_D1", "CCCC"), ("3_D1", "GGGG"),
        ("4_D1", "AAAT"), ("5_D1", "CCCT"), ("6_D1", "GGGT"),
    ])
    assert _ext.find_hors(aln, 1, 2, 3, 1) == [(1, 3, 4, 6, 1, 2)]


def test_mid_close_counts_all_pairs(tmp_path: Path) -> None:
    # Same three pairs, but a 7th repeat ends the run *before* the boundary,
    # so the close is a mid-loop close and total_variant is the full 3.
    aln = _write_aln(tmp_path / "cm.fasta", [
        ("1_D1", "AAAA"), ("2_D1", "CCCC"), ("3_D1", "GGGG"),
        ("4_D1", "AAAT"), ("5_D1", "CCCT"), ("6_D1", "GGGT"),
        ("7_D1", "TTTT"),
    ])
    assert _ext.find_hors(aln, 1, 2, 3, 1) == [(1, 3, 4, 6, 1, 3)]


def test_method2_split(tmp_path: Path) -> None:
    # Region A = repeats 1-3, region B = repeats 4-6; parallel match across split.
    aln = _write_aln(tmp_path / "d.fasta", [
        ("1_D1", "AAAA"), ("2_D1", "CCCC"), ("3_D1", "GGGG"),
        ("4_D1", "AAAA"), ("5_D1", "CCCC"), ("6_D1", "GGGG"),
    ])
    assert _ext.find_hors(aln, 3, 0, 2, 2) == [(1, 3, 4, 6, 1, 0)]


def test_wrapped_alignment_strips_newlines(tmp_path: Path) -> None:
    # 60-column wrapping must not change results: the scanner strips newlines
    # exactly like the binary. Same data as the parallel case, wrapped at 2.
    def wrap(s: str) -> str:
        return "\n".join(s[i:i + 2] for i in range(0, len(s), 2))
    recs = [("1_D1", "AAAA"), ("2_D1", "CCCC"), ("3_D1", "GGGG"),
            ("4_D1", "AAAA"), ("5_D1", "CCCC"), ("6_D1", "GGGG")]
    p = tmp_path / "w.fasta"
    with p.open("w") as f:
        for name, seq in recs:
            f.write(f">{name}\n{wrap(seq)}\n")
    assert _ext.find_hors(str(p), 1, 0, 3, 1) == [(1, 3, 4, 6, 1, 0)]


def test_no_hors_below_cutoff(tmp_path: Path) -> None:
    # A single matching pair can't form a HOR when cutoff = 3.
    aln = _write_aln(tmp_path / "n.fasta", [
        ("1_D1", "AAAA"), ("2_D1", "CCCC"), ("3_D1", "GGGG"), ("4_D1", "AAAA"),
    ])
    assert _ext.find_hors(aln, 1, 0, 3, 1) == []


def test_stream_matches_list(tmp_path: Path) -> None:
    """`find_hors_stream` (file, streaming) emits the same rows as `find_hors`
    (in-memory list), in the same order."""
    import numpy as np
    aln = _write_aln(tmp_path / "s.fasta", [
        ("1_D1", "AAAA"), ("2_D1", "CCCC"), ("3_D1", "GGGG"),
        ("4_D1", "AAAT"), ("5_D1", "CCCT"), ("6_D1", "GGGT"), ("7_D1", "TTTT"),
    ])
    lst = _ext.find_hors(aln, 1, 2, 3, 1)
    raw = tmp_path / "out.raw"
    n = _ext.find_hors_stream(aln, str(raw), 1, 2, 3, 1)
    rows = np.fromfile(raw, dtype=np.int32).reshape(-1, 6)
    assert n == len(lst)
    assert [tuple(int(x) for x in r) for r in rows] == [tuple(t) for t in lst]
