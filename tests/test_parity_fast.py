"""End-to-end parity test for the ``--fast`` flag against the upstream R
reference output.

The ``--fast`` path uses edlib for repeat mapping and a custom
in-process consensus instead of nhmmer + clustalo. It does not need to
hit the strict rotation-equivalence bar the default path does (small
SNP differences in the consensus can rotate the canonical anchor to a
different position), but the user's quality bar of ≥99.9 % on repeat
positions and scoring should hold. Specifically we assert:

* per-array boundaries match exactly (count, start, end, top_N, class)
* the representative is within 5 % Hamming after best rotation
  (same tolerance as ``compare_aregarrays``)
* every R-reported repeat is recovered within ±10 bp at ≥99.9 % recall

We re-use the small Ath Chr1 extraction fixture, which exercises the
two stages (Stage-1 cluster consensus + Stage-2 rescore consensus) on
realistic CEN178-like input.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import pytest

from .compare import compare_aregarrays, compare_regarrays


REPO_ROOT = Path(__file__).resolve().parent.parent
SMALL_FASTA = REPO_ROOT / "tests" / "data" / "ath_Chr1_extraction_trc.fasta"
SMALL_GOLDEN = REPO_ROOT / "tests" / "golden" / "ath_Chr1_extraction_trc_r_ref"


@pytest.fixture(scope="module")
def fast_output(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("python_output_fast")
    result = subprocess.run(
        [
            sys.executable, "-m", "trash_py",
            "-f", str(SMALL_FASTA),
            "-o", str(out),
            "--fast", "-q",
        ],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"pipeline failed:\n{result.stderr}"
    return out


def _names(name: str) -> tuple[Path, Path]:
    fname = f"{SMALL_FASTA.name}_{name}"
    return SMALL_GOLDEN / fname, fname


def test_regarrays_matches_exactly(fast_output: Path) -> None:
    # Stage 5 output is upstream of any consensus work.
    ref, fname = _names("regarrays.csv")
    diff = compare_regarrays(ref, fast_output / fname)
    assert not diff, diff.report()


def test_aregarrays_representative_within_5pct(fast_output: Path) -> None:
    # Same 5 % rotation-Hamming tolerance the default-path parity test
    # uses. The custom consensus is allowed to drift a handful of SNPs
    # vs the R reference but must stay below that bar.
    ref, fname = _names("aregarrays.csv")
    diff = compare_aregarrays(ref, fast_output / fname)
    assert not diff, diff.report()


def _read_repeats(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def test_repeat_positions_recall_99_9pct(fast_output: Path) -> None:
    """≥99.9 % of R-reported repeats are recovered within ±10 bp.

    This is the user's quality bar for the ``--fast`` mode. Matching is
    done per-repeat against the nearest candidate of the same strand,
    not via bucketed midpoints — adjacent buckets would otherwise
    spuriously miss positions sitting 1 bp on the wrong side of a
    bucket boundary.
    """
    ref_path, fname = _names("repeats.csv")
    got_path = fast_output / fname
    assert got_path.exists(), f"missing output: {got_path}"

    pos_tol = 10
    ref_rows = _read_repeats(ref_path)
    got_rows = _read_repeats(got_path)

    def midpoint(r: dict[str, str]) -> int:
        return (int(r["start"]) + int(r["end"])) // 2

    # Index got midpoints by (seqID, strand)
    got_by_key: dict[tuple[str, str], list[int]] = {}
    for g in got_rows:
        got_by_key.setdefault((g["seqID"], g["strand"]), []).append(midpoint(g))
    for v in got_by_key.values():
        v.sort()

    from bisect import bisect_left

    matched = 0
    for r in ref_rows:
        mids = got_by_key.get((r["seqID"], r["strand"]), [])
        if not mids:
            continue
        m = midpoint(r)
        i = bisect_left(mids, m)
        # Check neighbours on either side of insertion point
        best = min(
            (abs(mids[j] - m) for j in (i - 1, i) if 0 <= j < len(mids)),
            default=pos_tol + 1,
        )
        if best <= pos_tol:
            matched += 1

    recall = matched / max(1, len(ref_rows))
    assert recall >= 0.999, (
        f"position recall {recall:.4%} below 99.9 % "
        f"(matched {matched}/{len(ref_rows)})"
    )


def test_repeat_count_within_0_1pct(fast_output: Path) -> None:
    """Total repeat count must be within 0.1 % of the R reference."""
    ref_path, fname = _names("repeats.csv")
    got_path = fast_output / fname
    r_n = sum(1 for _ in _read_repeats(ref_path))
    g_n = sum(1 for _ in _read_repeats(got_path))
    drift = abs(r_n - g_n) / max(1, r_n)
    assert drift <= 0.001, f"repeat count drift {drift:.4%} (ref={r_n} got={g_n})"
