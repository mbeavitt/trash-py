from __future__ import annotations

from trash_py import arrays
from trash_py.arrays import ArrayBreaks, ArrayRow, CollapsedKmer


def test_split_and_check_arrays_reuses_consensus_score(monkeypatch) -> None:
    calls = {"seq_win_score_int": 0}

    def fake_seq_win_score_int(*args, **kwargs):
        calls["seq_win_score_int"] += 1
        return 0.0

    monkeypatch.setattr(
        arrays,
        "find_array_breaks",
        lambda *args, **kwargs: [ArrayBreaks(start=3, end=12, seqID="seq", numID=9)],
    )
    monkeypatch.setattr(
        arrays,
        "collapse_array_kmers",
        lambda *args, **kwargs: [CollapsedKmer(kmers=["a"], count=1, locations=[3, 9], distances=[])],
    )
    monkeypatch.setattr(arrays, "find_top_repeat_distances", lambda *args, **kwargs: (7, "N_7_Count_1.", [7]))
    monkeypatch.setattr(arrays, "build_consensus_repeat", lambda *args, **kwargs: (12.5, 7, "rep"))
    monkeypatch.setattr(arrays, "seq_win_score_int", fake_seq_win_score_int)

    rows = arrays.split_and_check_arrays(
        region_start=100,
        sequence="A" * 30,
        seqID="seq",
        numID=9,
    )

    assert rows == [
        ArrayRow(
            start=102,
            end=111,
            seqID="seq",
            numID=9,
            score=12.5,
            top_N=7,
            top_5_N="N_7_Count_1.",
            representative="rep",
        )
    ]
    assert calls["seq_win_score_int"] == 0


def test_split_and_check_arrays_falls_back_to_local_score_when_needed(monkeypatch) -> None:
    calls = {"seq_win_score_int": 0}

    def fake_seq_win_score_int(*args, **kwargs):
        calls["seq_win_score_int"] += 1
        return 25.0

    monkeypatch.setattr(
        arrays,
        "find_array_breaks",
        lambda *args, **kwargs: [ArrayBreaks(start=1, end=10, seqID="seq", numID=2)],
    )
    monkeypatch.setattr(
        arrays,
        "collapse_array_kmers",
        lambda *args, **kwargs: [CollapsedKmer(kmers=["a"], count=1, locations=[1, 5], distances=[])],
    )
    monkeypatch.setattr(arrays, "find_top_repeat_distances", lambda *args, **kwargs: (0, "", []))
    monkeypatch.setattr(arrays, "seq_win_score_int", fake_seq_win_score_int)

    rows = arrays.split_and_check_arrays(
        region_start=50,
        sequence="A" * 30,
        seqID="seq",
        numID=2,
    )

    assert rows == [
        ArrayRow(
            start=50,
            end=59,
            seqID="seq",
            numID=2,
            score=75.0,
            top_N=0,
            top_5_N="",
            representative="",
        )
    ]
    assert calls["seq_win_score_int"] == 1
