"""End-to-end parity test against the upstream R reference output for
the Ath_CEN1.fasta input.

Slow: runs the full pipeline on the chromosome 1 centromere extraction.
Skipped automatically when either the input or the golden output is
missing.
"""
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
CEN1_FASTA = REPO_ROOT / "tests" / "data" / "Ath_CEN1.fasta.gz"
CEN1_GOLDEN = REPO_ROOT / "tests" / "golden" / "Ath_CEN1_r_ref"


pytestmark = pytest.mark.skipif(
    not CEN1_FASTA.exists() or not CEN1_GOLDEN.exists(),
    reason="Ath_CEN1 fasta or golden directory not available",
)


@pytest.fixture(scope="module")
def python_output(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("python_output_cen1")
    result = subprocess.run(
        [sys.executable, "-m", "trash_py", "-f", str(CEN1_FASTA), "-o", str(out), "-q"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, f"pipeline failed:\n{result.stderr}"
    return out


def _names(name: str) -> tuple[Path, Path]:
    fname = f"{CEN1_FASTA.name}_{name}"
    return CEN1_GOLDEN / fname, fname


def _check(diff) -> None:
    assert not diff, diff.report()


def test_regarrays_matches(python_output: Path) -> None:
    ref, fname = _names("regarrays.csv")
    _check(compare_regarrays(ref, python_output / fname))


def test_aregarrays_matches(python_output: Path) -> None:
    ref, fname = _names("aregarrays.csv")
    _check(compare_aregarrays(ref, python_output / fname))


def test_classarrays_matches(python_output: Path) -> None:
    ref, fname = _names("classarrays.csv")
    _check(compare_classarrays(ref, python_output / fname))


def test_no_repeats_arrays_matches(python_output: Path) -> None:
    ref, fname = _names("no_repeats_arrays.csv")
    _check(compare_regarrays(ref, python_output / fname))


def test_arrays_matches(python_output: Path) -> None:
    ref, fname = _names("arrays.csv")
    _check(compare_arrays(ref, python_output / fname))


def test_repeats_matches(python_output: Path) -> None:
    ref, fname = _names("repeats.csv")
    _check(compare_repeats(ref, python_output / fname))


def test_repeats_with_seq_matches(python_output: Path) -> None:
    ref, fname = _names("repeats_with_seq.csv")
    _check(compare_repeats_with_seq(ref, python_output / fname))


def test_arrays_gff_matches(python_output: Path) -> None:
    ref, fname = _names("arrays.gff")
    _check(compare_gff(ref, python_output / fname))


def test_repeats_gff_matches(python_output: Path) -> None:
    ref, fname = _names("repeats.gff")
    _check(compare_gff(ref, python_output / fname))
