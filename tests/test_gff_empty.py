"""Regression test: exporting a GFF from an empty table must not crash.

``export_gff`` inspects ``rows[0]`` to tell an arrays table from a repeats
table. On a repeat-poor assembly every array can be filtered out before export,
leaving nothing to write; that previously raised ``IndexError`` and aborted the
run with five of the nine output files on disk. An empty table is a legitimate
result, so it now writes an empty GFF and the pipeline completes.
"""
from __future__ import annotations

from pathlib import Path

from trash_py.io_gff import export_gff


def test_export_gff_accepts_empty_rows(tmp_path: Path) -> None:
    out = tmp_path / "empty.gff"
    export_gff(
        [],
        out,
        seqid=1,
        source="TRASH",
        type_="Satellite_array",
        start=1,
        end=2,
        score=5,
        attributes=[9, 10, 14],
        attribute_names=["Name=", "Repeat_no=", "Repeat_median_width="],
    )
    assert out.exists()
    assert out.read_bytes() == b""


def test_export_gff_still_writes_rows(tmp_path: Path) -> None:
    # the guard must not change behaviour for a non-empty table
    out = tmp_path / "one.gff"
    export_gff(
        [{"start": 10, "end": 20, "seqID": "chr1", "class": "178_1"}],
        out,
        seqid="chr1",
        source="TRASH",
        type_="Satellite_DNA",
        start="10",
        end="20",
        attributes="Name=178_1",
    )
    body = out.read_bytes()
    assert body.endswith(b"\r")
    assert b"TRASH" in body and b"Satellite_DNA" in body
