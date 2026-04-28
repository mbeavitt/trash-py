"""End-to-end check for `-t templates.fasta` support."""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parent.parent
TEMPLATES = REPO_ROOT / "tests" / "data" / "test_templates.fasta"


def _run(fasta: Path, out: Path, templates: Path | None = None) -> subprocess.CompletedProcess:
    cmd = [sys.executable, "-m", "trash_py", "-f", str(fasta), "-o", str(out)]
    if templates is not None:
        cmd += ["-t", str(templates)]
    return subprocess.run(cmd, capture_output=True, text=True)


def _classes(csv_path: Path) -> set[str]:
    with csv_path.open() as f:
        return {r["class"] for r in csv.DictReader(f)}


@pytest.fixture(scope="module")
def small_fasta() -> Path:
    return REPO_ROOT / "tests" / "data" / "ath_Chr1_extraction_trc.fasta"


def test_pipeline_accepts_templates_flag(small_fasta: Path, tmp_path: Path) -> None:
    out = tmp_path / "with_templates"
    result = _run(small_fasta, out, templates=TEMPLATES)
    assert result.returncode == 0, f"pipeline failed:\n{result.stderr}"
    assert (out / f"{small_fasta.name}_classarrays.csv").exists()


def test_classes_use_template_names(small_fasta: Path, tmp_path: Path) -> None:
    """The small fasta has a ~159 bp satellite that should match the
    `CEN158` template; the resulting class column should reflect that
    rather than the auto-generated `<width>_<n>` name."""
    out = tmp_path / "with_templates"
    result = _run(small_fasta, out, templates=TEMPLATES)
    assert result.returncode == 0, result.stderr

    classes_with_t = _classes(out / f"{small_fasta.name}_classarrays.csv")
    assert "CEN158" in classes_with_t, (
        f"expected CEN158 in classes, got {classes_with_t}"
    )

    out2 = tmp_path / "no_templates"
    result2 = _run(small_fasta, out2)
    assert result2.returncode == 0, result2.stderr
    classes_no_t = _classes(out2 / f"{small_fasta.name}_classarrays.csv")
    assert classes_with_t != classes_no_t, (
        f"templates had no effect on classes: {classes_with_t}"
    )


def test_repeats_csv_uses_template_class(small_fasta: Path, tmp_path: Path) -> None:
    out = tmp_path / "with_templates"
    result = _run(small_fasta, out, templates=TEMPLATES)
    assert result.returncode == 0, result.stderr

    repeats = out / f"{small_fasta.name}_repeats.csv"
    assert repeats.exists()
    classes = _classes(repeats)
    assert "CEN158" in classes
