"""Parallel-mode parity: run trash-py with `-p 4` against the same fasta as
test_parity.py and apply all the stage comparators. Guards the parity
contract with `tests/golden/` for the parallel code path."""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from .compare import (
    compare_aregarrays,
    compare_arrays,
    compare_classarrays,
    compare_gff,
    compare_regarrays,
    compare_repeats,
    compare_repeats_with_seq,
)


REPO_ROOT = Path(__file__).resolve().parent.parent
SMALL_FASTA = REPO_ROOT / "tests" / "data" / "ath_Chr1_extraction_trc.fasta"
SMALL_GOLDEN = REPO_ROOT / "tests" / "golden" / "ath_Chr1_extraction_trc_r_ref"


@pytest.fixture(scope="module")
def python_output_parallel(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("python_output_small_parallel")
    result = subprocess.run(
        [sys.executable, "-m", "trash_py", "-f", str(SMALL_FASTA), "-o", str(out), "-p", "4"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"pipeline failed:\n{result.stderr}"
    return out


def _names(name: str) -> tuple[Path, Path]:
    fname = f"{SMALL_FASTA.name}_{name}"
    return SMALL_GOLDEN / fname, fname


def _check(diff) -> None:
    assert not diff, diff.report()


def test_regarrays_matches_parallel(python_output_parallel: Path) -> None:
    ref, fname = _names("regarrays.csv")
    _check(compare_regarrays(ref, python_output_parallel / fname))


def test_aregarrays_matches_parallel(python_output_parallel: Path) -> None:
    ref, fname = _names("aregarrays.csv")
    _check(compare_aregarrays(ref, python_output_parallel / fname))


def test_classarrays_matches_parallel(python_output_parallel: Path) -> None:
    ref, fname = _names("classarrays.csv")
    _check(compare_classarrays(ref, python_output_parallel / fname))


def test_no_repeats_arrays_matches_parallel(python_output_parallel: Path) -> None:
    ref, fname = _names("no_repeats_arrays.csv")
    _check(compare_regarrays(ref, python_output_parallel / fname))


def test_arrays_matches_parallel(python_output_parallel: Path) -> None:
    ref, fname = _names("arrays.csv")
    _check(compare_arrays(ref, python_output_parallel / fname))


def test_repeats_matches_parallel(python_output_parallel: Path) -> None:
    ref, fname = _names("repeats.csv")
    _check(compare_repeats(ref, python_output_parallel / fname))


def test_repeats_with_seq_matches_parallel(python_output_parallel: Path) -> None:
    ref, fname = _names("repeats_with_seq.csv")
    _check(compare_repeats_with_seq(ref, python_output_parallel / fname))


def test_arrays_gff_matches_parallel(python_output_parallel: Path) -> None:
    ref, fname = _names("arrays.gff")
    _check(compare_gff(ref, python_output_parallel / fname))


def test_repeats_gff_matches_parallel(python_output_parallel: Path) -> None:
    ref, fname = _names("repeats.gff")
    _check(compare_gff(ref, python_output_parallel / fname))
