"""External-tool precheck."""
from __future__ import annotations

import argparse
from pathlib import Path

import pytest

from trash_py import pipeline


def test_check_returns_missing(monkeypatch) -> None:
    monkeypatch.setattr(pipeline.shutil, "which", lambda name: None)
    missing = pipeline.check_external_tools()
    assert set(missing) == set(pipeline._EXTERNAL_TOOLS)


def test_check_returns_empty_when_all_present(monkeypatch) -> None:
    monkeypatch.setattr(pipeline.shutil, "which", lambda name: f"/fake/{name}")
    assert pipeline.check_external_tools() == []


def test_run_pipeline_exits_when_tools_missing(monkeypatch, tmp_path: Path, capsys) -> None:
    monkeypatch.setattr(pipeline.shutil, "which", lambda name: None)

    args = argparse.Namespace(
        fasta=tmp_path / "x.fasta",
        output=tmp_path / "out",
        max_rep_size=1000,
        min_rep_size=7,
        templates=None,
    )

    with pytest.raises(SystemExit) as excinfo:
        pipeline.run_pipeline(args)
    assert excinfo.value.code == 2

    err = capsys.readouterr().err
    assert "clustalo" in err
    assert "conda install" in err
