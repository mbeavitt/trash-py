"""Integration tests for the HOR pipeline (`trash_py.hor`).

The core test injects a committed alignment (so it does not depend on the local
MAFFT version) and requires the produced HOR table + repeats_with_hors table to
match committed golden files **byte-for-byte**. A second test runs the real
MAFFT path as a smoke test, and an opt-in test validates against the full
upstream reference output when it is available on disk.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

import trash_py.hor as hor
from trash_py.hor import HorArgs, run_hor_single
import trash_py._log as log

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "tests" / "data"
GOLD = REPO / "tests" / "golden" / "hor"
SAMPLE_REPEATS = DATA / "hor_sample_repeats.csv"
SAMPLE_ALN = DATA / "hor_sample_aligned.fasta"


def _inject_alignment(monkeypatch, aln: Path) -> None:
    def fake_align(repeats, out_dir, name):
        dst = Path(out_dir) / f"{name}temp.aligned.fasta"
        shutil.copy(aln, dst)
        return dst
    monkeypatch.setattr(hor, "align_repeats", fake_align)


def test_hor_tables_byte_identical(tmp_path: Path, monkeypatch) -> None:
    log.configure(quiet=True)
    _inject_alignment(monkeypatch, SAMPLE_ALN)
    args = HorArgs(repeats=SAMPLE_REPEATS, output_folder=tmp_path,
                   hor_class="178_1", make_plot=False)
    run_hor_single(args, "CP116282.1")

    for fname in ("HORs_178_1_CP116282.1.csv",
                  "repeats_with_hors_178_1_CP116282.1.csv"):
        got = (tmp_path / fname).read_bytes()
        ref = (GOLD / fname).read_bytes()
        assert got == ref, f"{fname} differs from golden"


def test_hor_plot_is_written(tmp_path: Path, monkeypatch) -> None:
    pytest.importorskip("matplotlib")
    log.configure(quiet=True)
    _inject_alignment(monkeypatch, SAMPLE_ALN)
    args = HorArgs(repeats=SAMPLE_REPEATS, output_folder=tmp_path,
                   hor_class="178_1", make_plot=True)
    run_hor_single(args, "CP116282.1")
    png = tmp_path / "HORs_lines_178_1_CP116282.1.png"
    assert png.exists() and png.stat().st_size > 0


@pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not installed")
def test_hor_subcommand_auto_selects_class(tmp_path: Path) -> None:
    """`trash-py hor <repeats> --chr-list <seq>` with no --class auto-picks the
    most abundant class and produces the expected (byte-identical) HOR table."""
    result = subprocess.run(
        [sys.executable, "-m", "trash_py", "hor", str(SAMPLE_REPEATS),
         "--chr-list", "CP116282.1", "-o", str(tmp_path), "--no-plot"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, result.stderr
    assert "auto-selected class '178_1'" in result.stdout
    got = (tmp_path / "HORs_178_1_CP116282.1.csv").read_bytes()
    ref = (GOLD / "HORs_178_1_CP116282.1.csv").read_bytes()
    assert got == ref


@pytest.mark.skipif(shutil.which("mafft") is None, reason="mafft not installed")
def test_full_mafft_path_runs(tmp_path: Path) -> None:
    """Smoke test of the real alignment path — HORs are found and tables written.
    (Byte-parity depends on the local MAFFT version, so we don't assert it here.)"""
    log.configure(quiet=True)
    args = HorArgs(repeats=SAMPLE_REPEATS, output_folder=tmp_path,
                   hor_class="178_1", make_plot=False)
    run_hor_single(args, "CP116282.1")
    assert (tmp_path / "HORs_178_1_CP116282.1.csv").exists()
    assert (tmp_path / "repeats_with_hors_178_1_CP116282.1.csv").exists()


# Opt-in full-reference parity: set HOR_REF_DIR to the upstream hor_output dir
# (and HOR_REF_REPEATS to the repeats_with_seq.csv). Validates a full chromosome
# byte-for-byte, modulo R's scientific-notation contractions which we keep long.
@pytest.mark.skipif(
    not os.environ.get("HOR_REF_DIR") or shutil.which("mafft") is None,
    reason="HOR_REF_DIR not set or mafft missing",
)
def test_full_reference_parity(tmp_path: Path) -> None:
    ref_dir = Path(os.environ["HOR_REF_DIR"])
    repeats = Path(os.environ["HOR_REF_REPEATS"])
    chrom = os.environ.get("HOR_REF_CHR", "CP116282.1")
    log.configure(quiet=True)
    args = HorArgs(repeats=repeats, output_folder=tmp_path,
                   hor_class="178_1", make_plot=False)
    run_hor_single(args, chrom)
    got = (tmp_path / f"HORs_178_1_{chrom}.csv").read_text().splitlines()
    ref = (ref_dir / f"HORs_178_1_{chrom}.csv").read_text().splitlines()
    assert len(got) == len(ref)
    diffs = [(g, r) for g, r in zip(got, ref) if g != r and _canon(g) != _canon(r)]
    assert not diffs, f"{len(diffs)} genuinely differing lines, e.g. {diffs[:3]}"


def _canon(line: str) -> str:
    """Normalise R's scientific contractions (e.g. 2e+05) to plain integers so
    the deliberate long-number difference doesn't count as a mismatch."""
    out = []
    for tok in line.split(","):
        t = tok.strip('"')
        try:
            f = float(t)
            out.append(str(int(f)) if f.is_integer() else repr(f))
        except ValueError:
            out.append(tok)
    return ",".join(out)
