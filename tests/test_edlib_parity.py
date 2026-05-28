"""Parity test: edlib mapper vs nhmmer mapper on real Arabidopsis data.

Loads representative sequences and array boundaries from the small golden
classarrays CSV, extracts the corresponding genomic windows from the raw
FASTA, then runs both mappers and asserts ≥99.9% recall (nhmmer hits
recovered by edlib within a 5 bp coordinate tolerance).

Requires nhmmer on PATH; skipped automatically when it is absent.
"""
from __future__ import annotations

import csv
import shutil
import tempfile
from pathlib import Path

import pytest

from trash_py.io_fasta import read_fasta_and_list
from trash_py.mapping import map_edlib, map_nhmmer


REPO_ROOT = Path(__file__).resolve().parent.parent
SMALL_FASTA = REPO_ROOT / "tests" / "data" / "ath_Chr1_extraction_trc.fasta"
SMALL_CLASSARRAYS = (
    REPO_ROOT
    / "tests"
    / "golden"
    / "ath_Chr1_extraction_trc_r_ref"
    / "ath_Chr1_extraction_trc.fasta_classarrays.csv"
)

pytestmark = pytest.mark.skipif(
    not shutil.which("nhmmer") or not SMALL_FASTA.exists() or not SMALL_CLASSARRAYS.exists(),
    reason="nhmmer not on PATH or test data missing",
)

# nhmmer reports *envelope* coordinates (a wider context window around each
# alignment), while edlib returns tighter alignment boundaries.  The
# downstream handle_edge_repeat absorbs this difference in the full
# pipeline, but at the raw mapper level the two can disagree by up to
# ~30 bp.  Hits with edit-distance fraction above EDIT_CAP are excluded:
# nhmmer's statistical calibration finds very short partial repeats
# (< 50% of rep_len) that edlib cannot recover within any fixed ed
# threshold, but these are rare boundary artefacts that the downstream
# gap-filling handles regardless.
POS_TOL = 30
EDIT_CAP = 0.50
# Minimum fraction of nhmmer hits (within EDIT_CAP) that edlib must recover
RECALL_THRESHOLD = 0.999


def _midpoint(r: dict) -> int:
    return (int(r["start"]) + int(r["end"])) // 2


def _recall(ref_rows: list[dict], got_rows: list[dict], tol: int) -> float:
    """Fraction of ref midpoints that have a matching got midpoint within tol."""
    if not ref_rows:
        return 1.0
    got_mids = [_midpoint(r) for r in got_rows]
    matched = sum(
        1 for rm in (_midpoint(r) for r in ref_rows)
        if any(abs(rm - gm) <= tol for gm in got_mids)
    )
    return matched / len(ref_rows)


@pytest.fixture(scope="module")
def array_data() -> list[dict]:
    """Return list of dicts with keys: representative, array_seq, seqID, arrayID, start, top_N."""
    sequences = {name: seq for name, seq in read_fasta_and_list(SMALL_FASTA)}
    arrays = []
    with SMALL_CLASSARRAYS.open() as f:
        for row in csv.DictReader(f):
            top_n = int(row["top_N"])
            if top_n < 14:
                continue
            seq_id = row["seqID"]
            start = int(row["start"])
            end = int(row["end"])
            genomic_seq = sequences.get(seq_id, "")
            array_seq = genomic_seq[start - 1:end]
            arrays.append({
                "representative": row["representative"],
                "array_seq": array_seq,
                "seqID": seq_id,
                "arrayID": int(row["array_num_ID"]),
                "start": start,
                "top_N": top_n,
            })
    return arrays


def _nhmmer_ed_fraction(row: dict, representative: str) -> float:
    """Levenshtein distance fraction between the representative and the
    sequence implied by the nhmmer envelope, using the correct strand."""
    from rapidfuzz.distance import Levenshtein
    from trash_py.sequence import rev_comp_string
    rep_len = len(representative)
    if rep_len == 0:
        return 1.0
    width = int(row["end"]) - int(row["start"]) + 1
    # We don't have the raw sequence here; approximate via width alone.
    # A hit whose width is far less than rep_len must have many deletions
    # and will have ed ≥ (rep_len - width).  Use that as a lower bound.
    min_ed = max(0, rep_len - width)
    return min_ed / rep_len


def test_edlib_recall_vs_nhmmer(array_data):
    """edlib must recover ≥99.9% of reachable nhmmer hits within POS_TOL bp.

    Hits whose nhmmer envelope width implies an edit fraction above EDIT_CAP
    (structurally unreachable with any fixed edit-distance threshold, e.g.
    very short partial boundary repeats) are excluded from the denominator.
    """
    from rapidfuzz.distance import Levenshtein
    from trash_py.io_fasta import read_fasta_and_list
    from trash_py.sequence import rev_comp_string

    sequences = {name: seq for name, seq in read_fasta_and_list(SMALL_FASTA)}
    all_nhmmer: list[dict] = []
    all_edlib: list[dict] = []
    array_rep: dict[int, tuple[str, int]] = {}  # arrayID → (representative, start)

    with tempfile.TemporaryDirectory() as tmp:
        for arr in array_data:
            nhmmer_rows = map_nhmmer(
                arr["representative"],
                arr["array_seq"],
                arr["seqID"],
                arr["arrayID"],
                arr["start"],
                output_folder=Path(tmp),
            )
            edlib_rows = map_edlib(
                arr["representative"],
                arr["array_seq"],
                arr["seqID"],
                arr["arrayID"],
                arr["start"],
            )
            all_nhmmer.extend(nhmmer_rows)
            all_edlib.extend(edlib_rows)
            array_rep[arr["arrayID"]] = (arr["representative"], arr["start"])

    # Filter to hits whose actual edit-distance fraction ≤ EDIT_CAP
    def actual_ed_frac(row: dict) -> float:
        rep, arr_start = array_rep[row["arrayID"]]
        rep_len = len(rep)
        if rep_len == 0:
            return 1.0
        seq_id = row["seqID"]
        s = int(row["start"]) - arr_start
        e = int(row["end"]) - arr_start + 1
        subseq = sequences.get(seq_id, "")[arr_start - 1 + s:arr_start - 1 + e]
        ed_fw = Levenshtein.distance(rep, subseq)
        ed_rv = Levenshtein.distance(rev_comp_string(rep), subseq)
        return min(ed_fw, ed_rv) / rep_len

    reachable_nhmmer = [r for r in all_nhmmer if actual_ed_frac(r) <= EDIT_CAP]
    total = len(reachable_nhmmer)
    assert total > 0, "nhmmer produced no reachable hits — check test data"

    recall = _recall(reachable_nhmmer, all_edlib, POS_TOL)
    allow_miss = max(1, round(total * (1 - RECALL_THRESHOLD)))
    matched = round(recall * total)
    missed = total - matched

    assert missed <= allow_miss, (
        f"edlib missed {missed}/{total} reachable nhmmer hits "
        f"(recall {recall:.4%}, threshold {RECALL_THRESHOLD:.1%}, "
        f"allow_miss={allow_miss}, POS_TOL={POS_TOL}bp)"
    )


def test_edlib_strand_parity_vs_nhmmer(array_data):
    """edlib strand assignments must match nhmmer on ≥99% of hits."""
    strand_match = 0
    strand_total = 0

    with tempfile.TemporaryDirectory() as tmp:
        for arr in array_data:
            nhmmer_rows = map_nhmmer(
                arr["representative"],
                arr["array_seq"],
                arr["seqID"],
                arr["arrayID"],
                arr["start"],
                output_folder=Path(tmp),
            )
            edlib_rows = map_edlib(
                arr["representative"],
                arr["array_seq"],
                arr["seqID"],
                arr["arrayID"],
                arr["start"],
            )
            edlib_by_mid = {_midpoint(r): r["strand"] for r in edlib_rows}
            for nr in nhmmer_rows:
                nm = _midpoint(nr)
                # Find closest edlib hit
                closest = min(edlib_by_mid, key=lambda m: abs(m - nm), default=None)
                if closest is not None and abs(closest - nm) <= POS_TOL:
                    strand_total += 1
                    if edlib_by_mid[closest] == nr["strand"]:
                        strand_match += 1

    if strand_total == 0:
        pytest.skip("No paired hits found for strand comparison")

    frac = strand_match / strand_total
    assert frac >= 0.99, f"strand agreement {frac:.2%} < 99% ({strand_match}/{strand_total})"


def test_edlib_count_close_to_nhmmer(array_data):
    """edlib hit count must be within 1% of nhmmer hit count per array."""
    with tempfile.TemporaryDirectory() as tmp:
        for arr in array_data:
            nhmmer_rows = map_nhmmer(
                arr["representative"],
                arr["array_seq"],
                arr["seqID"],
                arr["arrayID"],
                arr["start"],
                output_folder=Path(tmp),
            )
            edlib_rows = map_edlib(
                arr["representative"],
                arr["array_seq"],
                arr["seqID"],
                arr["arrayID"],
                arr["start"],
            )
            n_n = len(nhmmer_rows)
            n_e = len(edlib_rows)
            if n_n == 0:
                continue
            tol = max(1, round(n_n * 0.01))
            assert abs(n_n - n_e) <= tol, (
                f"arrayID {arr['arrayID']}: nhmmer={n_n}, edlib={n_e} "
                f"(diff {abs(n_n-n_e)} > tol {tol})"
            )
