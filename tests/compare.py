"""Comparators for R reference vs Python output.

Design notes:

Stage 5 (`*_regarrays.csv`) is byte-deterministic in the R reference, so
we hash-compare it.

Everything from stage 6 on is RNG-sensitive in the R reference (parallel
`sample()`/`runif()` calls don't inherit the parent seed). Array
boundaries (start/end/seqID/numID/top_N) are still stable, but the
chosen representative is a **random circular rotation** of the "true"
consensus, and downstream nhmmer hits depend on which rotation won. So
for stage-6+ files we:

* compare array boundaries exactly
* compare the representative up to circular rotation (canonicalised to
  the lexicographically smallest rotation)
* compare floats with a small absolute tolerance
* compare derived counts (repeats_number, etc.) with a small integer
  tolerance

Per-repeat position sets are compared loosely: we require the Python
output to recover a large fraction of the R-reported repeat positions
within a small bp tolerance.
"""
from __future__ import annotations

import csv
import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


# ---- low-level helpers -----------------------------------------------------

def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open() as f:
        return list(csv.DictReader(f))


def _int_r(x: str) -> int:
    """Parse an integer-valued field that R may have formatted in scientific
    notation (e.g. R's default `scipen=0` writes 900000 as "9e+05")."""
    try:
        return int(x)
    except ValueError:
        return int(float(x))


def canonical_rotation(seq: str) -> str:
    """Lex-smallest circular rotation of `seq`."""
    if not seq:
        return seq
    doubled = seq + seq
    n = len(seq)
    return min(doubled[i:i + n] for i in range(n))


def rotations_equivalent(a: str, b: str) -> bool:
    if len(a) != len(b):
        return False
    return canonical_rotation(a) == canonical_rotation(b)


def float_close(a: str, b: str, *, abs_tol: float = 1e-6, rel_tol: float = 1e-6) -> bool:
    fa, fb = float(a), float(b)
    return abs(fa - fb) <= max(abs_tol, rel_tol * max(abs(fa), abs(fb)))


# ---- per-file comparators --------------------------------------------------

@dataclass
class Diff:
    path: str
    messages: list[str]

    def __bool__(self) -> bool:
        return bool(self.messages)

    def report(self) -> str:
        if not self.messages:
            return f"{self.path}: OK"
        return f"{self.path}:\n  " + "\n  ".join(self.messages)


def _key_by(rows: list[dict[str, str]], keys: tuple[str, ...]) -> dict[tuple, dict[str, str]]:
    return {tuple(r[k] for k in keys): r for r in rows}


def compare_regarrays(ref: Path, got: Path) -> Diff:
    """Stage 5 output: deterministic, hash-compare."""
    msgs: list[str] = []
    if not got.exists():
        msgs.append(f"missing: {got}")
    elif sha256_of(ref) != sha256_of(got):
        msgs.append("sha256 mismatch")
    return Diff(ref.name, msgs)


def compare_aregarrays(ref: Path, got: Path, *, rep_hamming_tol_frac: float = 0.05) -> Diff:
    """Stage 6 output: array boundaries exact; representative up-to-rotation."""
    msgs: list[str] = []
    if not got.exists():
        return Diff(ref.name, [f"missing: {got}"])

    r_rows = read_csv(ref)
    g_rows = read_csv(got)
    if len(r_rows) != len(g_rows):
        msgs.append(f"row count: ref={len(r_rows)} got={len(g_rows)}")
        return Diff(ref.name, msgs)

    r_keyed = _key_by(r_rows, ("seqID", "start", "end"))
    g_keyed = _key_by(g_rows, ("seqID", "start", "end"))
    missing_keys = sorted(set(r_keyed) - set(g_keyed))
    if missing_keys:
        msgs.append(f"missing array boundaries: {missing_keys[:5]}{' …' if len(missing_keys) > 5 else ''}")
    for key, r_row in r_keyed.items():
        g_row = g_keyed.get(key)
        if g_row is None:
            continue
        if r_row["numID"] != g_row["numID"]:
            msgs.append(f"{key}: numID {r_row['numID']} vs {g_row['numID']}")
        if r_row["top_N"] != g_row["top_N"]:
            msgs.append(f"{key}: top_N {r_row['top_N']} vs {g_row['top_N']}")
        if not float_close(r_row["score"], g_row["score"], abs_tol=1e-4):
            msgs.append(f"{key}: score {r_row['score']} vs {g_row['score']}")
        # top_5_N: set of "N_x_Count_y" clauses separated by "."
        r_clauses = set(filter(None, r_row["top_5_N"].split(".")))
        g_clauses = set(filter(None, g_row["top_5_N"].split(".")))
        if r_clauses != g_clauses:
            msgs.append(f"{key}: top_5_N clauses differ")
        # Representative: rotation-equivalent, tolerate a small hamming distance
        rr = r_row["representative"]
        gr = g_row["representative"]
        if len(rr) != len(gr):
            msgs.append(f"{key}: representative len {len(rr)} vs {len(gr)}")
        elif rr and gr:
            crr = canonical_rotation(rr)
            cgr = canonical_rotation(gr)
            if crr != cgr:
                # Compute min hamming across all rotations of gr vs rr
                min_d = min(
                    sum(x != y for x, y in zip(rr, (gr * 2)[i:i + len(gr)]))
                    for i in range(len(gr))
                )
                if min_d > rep_hamming_tol_frac * len(rr):
                    msgs.append(
                        f"{key}: representative differs: min hamming={min_d}/{len(rr)} "
                        f"(tol={int(rep_hamming_tol_frac * len(rr))})"
                    )
    return Diff(ref.name, msgs)


def compare_arrays(ref: Path, got: Path, *, count_tol: int = 5, width_tol: int = 5, score_tol: int = 10) -> Diff:
    """Stage 11 final arrays.csv: boundaries exact + class + fuzzy aggregates."""
    msgs: list[str] = []
    if not got.exists():
        return Diff(ref.name, [f"missing: {got}"])
    r_rows = read_csv(ref)
    g_rows = read_csv(got)
    if len(r_rows) != len(g_rows):
        msgs.append(f"row count: ref={len(r_rows)} got={len(g_rows)}")
        return Diff(ref.name, msgs)
    r_keyed = _key_by(r_rows, ("seqID", "start", "end"))
    g_keyed = _key_by(g_rows, ("seqID", "start", "end"))
    for key, r_row in r_keyed.items():
        g_row = g_keyed.get(key)
        if g_row is None:
            msgs.append(f"missing boundary: {key}")
            continue
        for col in ("numID", "top_N", "class"):
            if r_row[col] != g_row[col]:
                msgs.append(f"{key}: {col} {r_row[col]} vs {g_row[col]}")
        if not rotations_equivalent(r_row["representative"], g_row["representative"]):
            msgs.append(f"{key}: representative not rotation-equivalent")
        for col, tol in (("repeats_number", count_tol), ("median_repeat_width", width_tol), ("median_score", score_tol)):
            if abs(int(r_row[col]) - int(g_row[col])) > tol:
                msgs.append(f"{key}: {col} {r_row[col]} vs {g_row[col]} (tol {tol})")
    return Diff(ref.name, msgs)


def compare_classarrays(ref: Path, got: Path) -> Diff:
    """Stage 8 classarrays.csv: boundaries + class exact."""
    msgs: list[str] = []
    if not got.exists():
        return Diff(ref.name, [f"missing: {got}"])
    r_rows = read_csv(ref)
    g_rows = read_csv(got)
    if len(r_rows) != len(g_rows):
        msgs.append(f"row count: ref={len(r_rows)} got={len(g_rows)}")
        return Diff(ref.name, msgs)
    r_keyed = _key_by(r_rows, ("seqID", "start", "end"))
    g_keyed = _key_by(g_rows, ("seqID", "start", "end"))
    for key, r_row in r_keyed.items():
        g_row = g_keyed.get(key)
        if g_row is None:
            msgs.append(f"missing boundary: {key}")
            continue
        for col in ("numID", "top_N", "class"):
            if r_row[col] != g_row[col]:
                msgs.append(f"{key}: {col} {r_row[col]} vs {g_row[col]}")
        if not rotations_equivalent(r_row["representative"], g_row["representative"]):
            msgs.append(f"{key}: representative not rotation-equivalent")
    return Diff(ref.name, msgs)


def compare_repeats(ref: Path, got: Path, *, count_tol_frac: float = 0.05, pos_tol: int = 10) -> Diff:
    """Stage 9/13 repeats.csv: count within frac; ≥95% of R positions recovered within pos_tol."""
    msgs: list[str] = []
    if not got.exists():
        return Diff(ref.name, [f"missing: {got}"])
    r_rows = read_csv(ref)
    g_rows = read_csv(got)

    r_n, g_n = len(r_rows), len(g_rows)
    if abs(r_n - g_n) > max(2, int(count_tol_frac * r_n)):
        msgs.append(f"row count: ref={r_n} got={g_n} (tol {count_tol_frac*100:.0f}%)")

    # Overlap of (seqID, strand, mid-point rounded to pos_tol granularity)
    def midpoints(rows: list[dict[str, str]]) -> set[tuple[str, str, int]]:
        return {
            (r["seqID"], r["strand"], (_int_r(r["start"]) + _int_r(r["end"])) // 2 // pos_tol * pos_tol)
            for r in rows
        }

    r_mid = midpoints(r_rows)
    g_mid = midpoints(g_rows)
    recall = len(r_mid & g_mid) / max(1, len(r_mid))
    if recall < 0.95:
        msgs.append(f"position recall {recall:.2%} (<95%) (shared {len(r_mid & g_mid)}/{len(r_mid)})")

    # Classes should be a superset / subset match
    r_classes = {r["class"] for r in r_rows}
    g_classes = {r["class"] for r in g_rows}
    if r_classes != g_classes:
        msgs.append(f"class set: ref={sorted(r_classes)} got={sorted(g_classes)}")
    return Diff(ref.name, msgs)


def compare_repeats_with_seq(ref: Path, got: Path, *, pos_tol: int = 10) -> Diff:
    """Same as compare_repeats + sanity-check that the extracted sequence
    column length matches (start, end) width on shared positions.
    """
    base = compare_repeats(ref, got, pos_tol=pos_tol)
    if not got.exists():
        return base
    msgs = list(base.messages)
    g_rows = read_csv(got)
    for r in g_rows:
        expected_len = _int_r(r["end"]) - _int_r(r["start"]) + 1
        if len(r["sequence"]) != expected_len:
            msgs.append(f"seq len mismatch at arrayID={r['arrayID']} start={r['start']}")
            break
    return Diff(ref.name, msgs)


def compare_gff(ref: Path, got: Path, *, count_tol_frac: float = 0.05, pos_tol: int = 10) -> Diff:
    """Very forgiving GFF comparator: same seqid+type+strand population,
    similar count, majority of positions overlap within tolerance.
    """
    msgs: list[str] = []
    if not got.exists():
        return Diff(ref.name, [f"missing: {got}"])

    def parse(path: Path) -> list[tuple[str, str, int, int, str]]:
        recs = []
        for line in path.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 7:
                continue
            recs.append((parts[0], parts[2], int(parts[3]), int(parts[4]), parts[6]))
        return recs

    r_recs, g_recs = parse(ref), parse(got)
    if abs(len(r_recs) - len(g_recs)) > max(2, int(count_tol_frac * max(1, len(r_recs)))):
        msgs.append(f"record count: ref={len(r_recs)} got={len(g_recs)}")
    r_buckets = {(r[0], r[1], r[4], r[2] // pos_tol * pos_tol) for r in r_recs}
    g_buckets = {(g[0], g[1], g[4], g[2] // pos_tol * pos_tol) for g in g_recs}
    recall = len(r_buckets & g_buckets) / max(1, len(r_buckets))
    if recall < 0.95:
        msgs.append(f"position recall {recall:.2%} (<95%)")
    return Diff(ref.name, msgs)


# ---- top-level dispatcher --------------------------------------------------

ComparatorFn = callable

FILE_COMPARATORS: dict[str, ComparatorFn] = {
    "_regarrays.csv": compare_regarrays,
    "_aregarrays.csv": compare_aregarrays,
    "_classarrays.csv": compare_classarrays,
    "_arrays.csv": compare_arrays,
    "_no_repeats_arrays.csv": compare_regarrays,  # tiny/empty — exact compare is fine
    "_repeats.csv": compare_repeats,
    "_repeats_with_seq.csv": compare_repeats_with_seq,
    "_arrays.gff": compare_gff,
    "_repeats.gff": compare_gff,
}


def compare_all(fasta_name: str, ref_dir: Path, got_dir: Path) -> list[Diff]:
    out: list[Diff] = []
    for suffix, fn in FILE_COMPARATORS.items():
        ref = ref_dir / f"{fasta_name}{suffix}"
        got = got_dir / f"{fasta_name}{suffix}"
        if not ref.exists():
            out.append(Diff(ref.name, [f"reference missing: {ref}"]))
            continue
        out.append(fn(ref, got))
    return out
