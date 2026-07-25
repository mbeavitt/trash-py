"""Regression test: the genome CLI must not crash on sub-window scaffolds.

Genome assemblies routinely carry short unplaced scaffolds. A sequence shorter
than one scoring window cannot host a scored tandem array, so the pipeline skips
it with a warning rather than crashing inside ``seq_win_score_int`` with
"window extends outside sequence bounds" (which previously aborted the whole
run). The reads API, by contrast, deliberately raises on short input -- that
contract is covered in ``test_reads_api.py``.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent
SMALL_FASTA = REPO_ROOT / "tests" / "data" / "ath_Chr1_extraction_trc.fasta"


def _read_fasta(path: Path) -> str:
    return path.read_text()


def test_cli_skips_short_scaffold(tmp_path: Path) -> None:
    # A real repeat-bearing chromosome plus a sub-window scaffold (800 bp,
    # below the ~1,111 bp default window). The run must complete (exit 0) and
    # still annotate the real chromosome.
    mixed = tmp_path / "mixed.fasta"
    mixed.write_text(_read_fasta(SMALL_FASTA).rstrip("\n") + "\n"
                     ">tiny_scaffold\n" + "ACGT" * 200 + "\n")
    out = tmp_path / "out"
    result = subprocess.run(
        [sys.executable, "-m", "trash_py", "-f", str(mixed), "-o", str(out)],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"pipeline crashed:\n{result.stderr}"
    combined = result.stdout + result.stderr
    assert "tiny_scaffold" in combined and "skipped" in combined
    # The real chromosome was still processed: array outputs exist and the
    # crashing message never appears.
    assert "window extends outside sequence bounds" not in combined
    assert (out / f"{mixed.name}_arrays.csv").exists()
