"""Map each array's representative onto the array's genomic span and
emit a list of repeat occurrences.

Pipeline (per array):
* `map_repeats` — `nhmmer` if the representative is long enough
  (top_N ≥ 14), else an exact-length Hamming-tolerant scan.
* `clean_overlaps` — drop overlapping hits and fill obvious tandem gaps.
* `rescore_repeats` — recompute the representative via clustal
  alignment, fill remaining gaps, rescore each repeat, and merge
  adjacent split repeats.
* `handle_edge_repeat` — tweak array-edge boundaries to maximise
  alignment with the representative.
"""
from __future__ import annotations

import math
import shutil
import subprocess
from pathlib import Path
from typing import Any

from rapidfuzz.distance import Levenshtein

from .arrays import _clustalo_align, _consensus_N
from .sequence import rev_comp_string


_COMP_ONLY = str.maketrans("acgtn", "tgcan")


RepeatRow = dict[str, Any]

NHMMER_ARGS = ("--popen", "0.1", "--pextend", "0.8", "--dna")


def _fasta_bytes(header: str, sequence: str, width: int = 60) -> str:
    wrapped = "\n".join(sequence[i:i + width] for i in range(0, len(sequence), width))
    return f">{header}\n{wrapped}\n"


def _parse_nhmmer_tblout(path: Path, seqID: str, arrayID: int) -> list[RepeatRow]:
    """Port of `read_and_format_nhmmer.R`.

    nhmmer --tblout columns (DNA search): target, accession, query,
    accession, hmmfrom, hmm to, alifrom, ali to, envfrom, env to, sq len,
    strand, E-value, score, bias, description of target.
    """
    rows: list[RepeatRow] = []
    if not path.exists():
        return rows
    for line in path.read_text().splitlines():
        if "#" in line:
            continue
        stripped = line.strip()
        if not stripped:
            continue
        fields = stripped.split()
        if len(fields) < 14:
            continue
        alifrom = int(fields[6])
        ali_to = int(fields[7])
        envfrom = int(fields[8])
        env_to = int(fields[9])
        strand = fields[11]
        evalue = float(fields[12])
        score = float(fields[13])

        # R swaps the "from"/"to" pair for "-" strand so start < end.
        if strand == "-":
            alifrom, ali_to = ali_to, alifrom
            envfrom, env_to = env_to, envfrom

        rows.append({
            "seqID": seqID,
            "arrayID": arrayID,
            "start": envfrom,
            "end": env_to,
            "strand": strand,
            "score": score,
            "eval": evalue,
        })
    return rows


def map_nhmmer(
    representative: str,
    array_sequence: str,
    seqID: str,
    arrayID: int,
    start_offset: int,
    output_folder: Path,
    nhmmer_exe: str = "nhmmer",
) -> list[RepeatRow]:
    """Run nhmmer against the array sequence and parse --tblout."""
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)
    rep_file = output_folder / f"Array_{arrayID}_repeat.fasta"
    seq_file = output_folder / f"Array_{arrayID}_sequence.fasta"
    tbl_file = output_folder / f"nhmmer_{arrayID}_output.txt"

    rep_file.write_text(_fasta_bytes("reference_repeat", representative))
    seq_file.write_text(_fasta_bytes(seqID, array_sequence))
    try:
        subprocess.run(
            [nhmmer_exe, *NHMMER_ARGS, "--tblout", str(tbl_file), str(rep_file), str(seq_file)],
            capture_output=True,
            check=True,
        )
        rows = _parse_nhmmer_tblout(tbl_file, seqID, arrayID)
    finally:
        for p in (rep_file, seq_file, tbl_file):
            if p.exists():
                p.unlink()

    for r in rows:
        r["start"] += start_offset - 1
        r["end"] += start_offset - 1
    return rows


def _match_pattern_mismatch(pattern: str, subject: str, max_mismatch: int) -> list[tuple[int, int]]:
    """Hamming-tolerant exact-length match — mirrors Biostrings::matchPattern
    with default `fixed=TRUE` (no IUPAC expansion, no indels)."""
    n = len(pattern)
    L = len(subject)
    hits: list[tuple[int, int]] = []
    if n == 0 or L < n:
        return hits
    for i in range(L - n + 1):
        mm = 0
        window = subject[i:i + n]
        for a, b in zip(pattern, window):
            if a != b:
                mm += 1
                if mm > max_mismatch:
                    break
        if mm <= max_mismatch:
            hits.append((i + 1, i + n))  # 1-based inclusive
    return hits


def map_default(
    representative: str,
    array_sequence: str,
    seqID: str,
    arrayID: int,
    start_offset: int,
) -> list[RepeatRow]:
    """Hamming-tolerant mapping for short (top_N < 14) repeats."""
    max_mismatch = math.ceil(len(representative) / 10)
    rows: list[RepeatRow] = []

    for s, e in _match_pattern_mismatch(representative, array_sequence, max_mismatch):
        rows.append({
            "seqID": seqID, "arrayID": arrayID,
            "start": s + start_offset - 1,
            "end": e + start_offset - 1,
            "strand": "+", "score": 0.0, "eval": -1.0,
        })

    rev = rev_comp_string(representative)
    for s, e in _match_pattern_mismatch(rev, array_sequence, max_mismatch):
        rows.append({
            "seqID": seqID, "arrayID": arrayID,
            "start": s + start_offset - 1,
            "end": e + start_offset - 1,
            "strand": "-", "score": 0.0, "eval": -1.0,
        })

    return rows


def map_repeats(
    representative: str,
    top_N: int,
    array_sequence: str,
    seqID: str,
    arrayID: int,
    start_offset: int,
    output_folder: Path | None = None,
    nhmmer_exe: str | None = None,
) -> list[RepeatRow]:
    """Dispatch to nhmmer or default matcher based on `top_N`."""
    if top_N >= 14:
        if output_folder is None:
            raise ValueError("output_folder required for nhmmer path")
        return map_nhmmer(
            representative, array_sequence, seqID, arrayID, start_offset,
            output_folder, nhmmer_exe or shutil.which("nhmmer") or "nhmmer",
        )
    return map_default(representative, array_sequence, seqID, arrayID, start_offset)


def _argmax_first(values: list[float]) -> int:
    best_idx = -1
    best_val = -math.inf
    for i, v in enumerate(values):
        if v > best_val:
            best_val = v
            best_idx = i
    return best_idx


def _argmin_first(values: list[float]) -> int:
    best_idx = -1
    best_val = math.inf
    for i, v in enumerate(values):
        if v < best_val:
            best_val = v
            best_idx = i
    return best_idx


def handle_overlaps(rows: list[RepeatRow], overlap_threshold: float = 0.1) -> list[RepeatRow]:
    """Port `src/repeats/handle_overlaps.R`.

    Returns rows with columns: seqID, arrayID, start, end, strand, score,
    eval, width, class. Mutates score/eval (to 0/-1) for partial overlaps
    that are sliced rather than dropped.
    """
    if len(rows) <= 1:
        return [dict(r) for r in rows]

    table = sorted((dict(r) for r in rows), key=lambda r: int(r["start"]))
    for r in table:
        r["start"] = int(r["start"])
        r["end"] = int(r["end"])
        r["width"] = int(r["width"])
        r["score"] = float(r["score"])
        r["eval"] = float(r["eval"])

    def recompute_overlap():
        for idx in range(len(table) - 1):
            ov = table[idx]["end"] - table[idx + 1]["start"] + 1
            table[idx]["overlap_with_next"] = max(ov, 0)
        table[-1]["overlap_with_next"] = 0

    recompute_overlap()

    while sum(r["overlap_with_next"] for r in table) > 0:
        i = _argmax_first([r["overlap_with_next"] for r in table])
        overlap = table[i]["overlap_with_next"]
        min_width = min(table[i]["width"], table[i + 1]["width"])
        overlap_fraction = overlap / min_width
        if overlap_fraction >= overlap_threshold:
            if -1 in (table[i]["eval"], table[i + 1]["eval"]):
                # Default-matcher path: eval is -1 sentinel; R drops lower
                # score via `which.min` — ties resolve to the earlier row.
                drop = i if table[i]["score"] <= table[i + 1]["score"] else i + 1
                table.pop(drop)
            else:
                # nhmmer path: higher eval is worse; `which.max` picks the
                # earlier row on tie.
                drop = i if table[i]["eval"] >= table[i + 1]["eval"] else i + 1
                table.pop(drop)
        else:
            table[i + 1]["start"] += overlap // 2
            table[i + 1]["score"] = 0.0
            table[i + 1]["eval"] = -1.0
            table[i]["end"] -= -(-overlap // 2)  # ceil
            table[i]["score"] = 0.0
            table[i]["eval"] = -1.0

        for r in table:
            r["overlap_with_next"] = 0
        if len(table) > 1:
            recompute_overlap()

    for r in table:
        r.pop("overlap_with_next", None)

    return [
        {k: r[k] for k in ("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class")}
        for r in table
    ]


def handle_gaps(rows: list[RepeatRow], representative_len: int) -> list[RepeatRow]:
    """Port `src/repeats/handle_gaps.R`.

    Extends repeats across small gaps (< rep_len / 2); drops orphan
    repeats separated by larger gaps.
    """
    if not rows:
        return []

    table = sorted((dict(r) for r in rows), key=lambda r: int(r["start"]))
    for r in table:
        r["start"] = int(r["start"])
        r["end"] = int(r["end"])

    gap = [0] * len(table)
    for i in range(len(table) - 1):
        gap[i] = max(table[i + 1]["start"] - table[i]["end"] - 1, 0)
    gap[-1] = 9_999_999

    to_remove: list[int] = []
    for i in range(len(table)):
        if gap[i] <= 0:
            continue
        if gap[i] < representative_len / 2:
            # Extend the upstream repeat (by strand).
            if table[i]["strand"] == "+":
                table[i]["end"] += gap[i]
            else:
                if i + 1 < len(table):
                    table[i + 1]["start"] -= gap[i]
            gap[i] = 0
        else:
            if i == 0 or gap[i - 1] != 0:
                to_remove.append(i)

    if to_remove:
        keep = [r for idx, r in enumerate(table) if idx not in set(to_remove)]
        table = keep

    return [
        {k: r[k] for k in ("seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class")}
        for r in table
    ]


def clean_overlaps(
    rows: list[RepeatRow],
    representative_len: int,
    arr_class: str,
) -> list[RepeatRow]:
    """Apply the per-array mapping cleanup block from main.R:386-395.

    Adds width + class, runs handle_overlaps (0.1), returns empty if < 3
    rows survive, then runs handle_gaps.
    """
    working: list[RepeatRow] = []
    for r in rows:
        w = int(r["end"]) - int(r["start"]) + 1
        working.append({**r, "width": w, "class": arr_class})

    working = handle_overlaps(working, overlap_threshold=0.1)
    if len(working) < 3:
        return working
    return handle_gaps(working, representative_len=representative_len)


def _evenly_spaced_ids(sample_ids: list[int], k: int) -> list[int]:
    """R: `sample_ids[round(seq(1, n, length.out=k))]` (deterministic stand-in
    for `sample`, shared with the stage 6 replacement)."""
    n = len(sample_ids)
    if n <= k:
        return list(sample_ids)
    if k == 1:
        return [sample_ids[0]]
    step = (n - 1) / (k - 1)
    picks = []
    for i in range(k):
        # R `round()` uses banker's rounding (round-half-to-even). Matches
        # what stage 6 uses; see `src/array/split_and_check_arrays.R:339`.
        x = 1 + i * step
        frac = x - math.floor(x)
        if frac == 0.5:
            base = math.floor(x)
            idx = base if base % 2 == 0 else base + 1
        else:
            idx = round(x)
        picks.append(sample_ids[idx - 1])
    return picks


def fill_gaps(
    rows: list[RepeatRow],
    array_sequence: str,
    array_start: int,
) -> list[RepeatRow]:
    """Port `src/repeats/fill_gaps.R`.

    Walks the sorted repeat list, and for any gap whose length sits within
    50% of the mean repeat width, re-checks the gap sequence (fw + rv) for
    similarity to the current representative. If edit distance / rep_len
    is below 0.8, appends a new repeat row for that gap.
    """
    gap_range = 0.5
    max_score = 0.8
    if len(rows) < 2:
        return [dict(r) for r in rows]

    table = [dict(r) for r in rows]

    widths = [int(r["end"]) - int(r["start"]) + 1 for r in table]
    mean_rep_len = round(sum(widths) / len(widths))
    lo = math.floor(mean_rep_len * (1 - gap_range))
    hi = math.ceil(mean_rep_len * (1 + gap_range))

    i = 1
    while i < len(table):
        gap_len = int(table[i]["start"]) - int(table[i - 1]["end"]) - 1
        if lo <= gap_len <= hi:
            new_start = int(table[i - 1]["end"]) + 1
            new_end = int(table[i]["start"]) - 1
            check_fw = array_sequence[new_start - array_start:new_end - array_start + 1]
            check_rv = rev_comp_string(check_fw)
            rep_i = table[i]["representative"]
            rep_len = len(rep_i)
            score_fw = Levenshtein.distance(rep_i, check_fw) / rep_len
            score_rv = Levenshtein.distance(rep_i, check_rv) / rep_len
            if min(score_fw, score_rv) < max_score:
                strand = "+" if score_fw <= score_rv else "-"
                table.append({
                    "seqID": table[i]["seqID"],
                    "arrayID": table[i]["arrayID"],
                    "start": new_start,
                    "end": new_end,
                    "strand": strand,
                    "score": -1.0,
                    "eval": -1.0,
                    "width": new_end - new_start + 1,
                    "class": table[i]["class"],
                    "representative": rep_i,
                    "score_template": min(score_fw, score_rv),
                })
                table.sort(key=lambda r: int(r["start"]))
                i += 1
        i += 1

    return table


def _adist(pattern: str, target: str) -> int:
    return Levenshtein.distance(pattern, target)


def rescore_repeats(
    rows: list[RepeatRow],
    arr_representative: str,
    arr_class: str,
    top_N: int,
    array_sequence: str,
    array_start: int,
    templates: dict[str, str] | None = None,
    clustalo_exe: str = "clustalo",
) -> list[RepeatRow]:
    """Ports main.R:401-473: recalculate representative via MSA, fill
    gaps, adist-rescore, then merge adjacent high-score short pairs.
    """
    max_repeats_to_align = 15
    min_repeats_to_recalculate = 10
    score_min_to_merge = 30.0
    size_max_to_merge = 1.0
    templates = templates or {}

    def _seq(start: int, end: int) -> str:
        return array_sequence[start - array_start:end - array_start + 1]

    table = [dict(r) for r in rows]
    for r in table:
        r["representative"] = arr_representative
        r["score_template"] = -1.0

    sample_ids = [i for i, r in enumerate(table) if r["strand"] != "."]
    if len(sample_ids) >= min_repeats_to_recalculate:
        if len(sample_ids) > max_repeats_to_align:
            sample_ids = _evenly_spaced_ids(sample_ids, max_repeats_to_align)
        sample_seqs = [_seq(int(table[j]["start"]), int(table[j]["end"])) for j in sample_ids]
        strands = [table[j]["strand"] for j in sample_ids]
        sample_seqs = [
            rev_comp_string(s) if st == "-" else s
            for s, st in zip(sample_seqs, strands)
        ]
        alignment = _clustalo_align(sample_seqs, clustalo_exe)
        consensus = _consensus_N(alignment, top_N)
        if consensus:
            for r in table:
                r["representative"] = consensus

    table = fill_gaps(table, array_sequence, array_start)

    # Re-extract sequences after fill_gaps (new rows may be present, and
    # order may have changed).
    rep0 = table[0]["representative"]
    rep0_len = len(rep0)
    rep0_rc = rev_comp_string(rep0)
    seqs = [_seq(int(r["start"]), int(r["end"])) for r in table]
    for idx, r in enumerate(table):
        if r["strand"] == "+":
            r["score"] = _adist(rep0, seqs[idx]) / rep0_len * 100
        elif r["strand"] == "-":
            r["score"] = _adist(rep0_rc, seqs[idx]) / rep0_len * 100

    template = templates.get(arr_class)
    if template is not None:
        tpl_rc = rev_comp_string(template)
        tpl_len = len(template)
        for idx, r in enumerate(table):
            if r["strand"] == "+":
                r["score_template"] = _adist(template, seqs[idx]) / tpl_len * 100
            elif r["strand"] == "-":
                r["score_template"] = _adist(tpl_rc, seqs[idx]) / tpl_len * 100

    # Merge adjacent split repeats (same strand, both > min_score, combined
    # width < rep_len). `i_r` iterates in place on a shrinking table.
    i_r = 0
    while i_r < len(table) - 1:
        a, b = table[i_r], table[i_r + 1]
        if (
            a["score"] > score_min_to_merge
            and b["score"] > score_min_to_merge
            and (int(a["width"]) + int(b["width"]) < size_max_to_merge * rep0_len)
        ):
            if a["strand"] == b["strand"] == "+":
                new_score = _adist(rep0, seqs[i_r] + seqs[i_r + 1]) / rep0_len * 100
                if new_score < min(a["score"], b["score"]):
                    a["end"] = int(b["end"])
                    table.pop(i_r + 1)
                    seqs.pop(i_r + 1)
                    a["width"] = int(a["end"]) - int(a["start"]) + 1
                    seqs[i_r] = _seq(int(a["start"]), int(a["end"]))
                    a["score"] = new_score
                    if arr_class in templates:
                        a["score_template"] = _adist(templates[arr_class], seqs[i_r]) / len(templates[arr_class]) * 100
            elif a["strand"] == b["strand"] == "-":
                new_score = _adist(rep0_rc, seqs[i_r] + seqs[i_r + 1]) / rep0_len * 100
                if new_score < min(a["score"], b["score"]):
                    a["end"] = int(b["end"])
                    table.pop(i_r + 1)
                    seqs.pop(i_r + 1)
                    a["width"] = int(a["end"]) - int(a["start"]) + 1
                    seqs[i_r] = _seq(int(a["start"]), int(a["end"]))
                    a["score"] = new_score
                    if arr_class in templates:
                        a["score_template"] = _adist(rev_comp_string(templates[arr_class]), seqs[i_r]) / len(templates[arr_class]) * 100
        i_r += 1

    return table


def find_edge_best_start_end(sequence_potential: str, representative: str) -> int:
    """Port `src/repeats/find_edge_best_start_end.R`.

    Scans a candidate edge window for a plateau point — the index at
    which kmer-match score stops accumulating new unique values relative
    to the rest of the window.
    """
    kmer = 12
    L_s = len(sequence_potential)
    L_r = len(representative)
    if L_s <= kmer or L_r <= kmer:
        return 0

    # R: `1:(length(x) - kmer)` — yields `length - kmer` kmers, NOT
    # `length - kmer + 1`. Off-by-one kept for parity.
    rep_kmers = {representative[i:i + kmer] for i in range(L_r - kmer)}
    pot_kmers = [sequence_potential[i:i + kmer] for i in range(L_s - kmer)]
    found = [i + 1 for i, km in enumerate(pot_kmers) if km in rep_kmers]

    # `scores[i]` = count of matched kmers whose 1-based index <= i.
    scores = [0] * L_s
    # Two-pointer since `found` is already sorted ascending.
    fj = 0
    for i in range(1, L_s + 1):
        while fj < len(found) and found[fj] <= i:
            fj += 1
        scores[i - 1] = fj

    division_score = [0.0] * L_s
    for i in range(1, L_s + 1):
        left_unique = len(set(scores[0:i]))
        right_unique = len(set(scores[i - 1:L_s]))
        numer = left_unique / i
        right_len = L_s - i + 1
        denom = right_unique / right_len if right_len > 0 else 0
        division_score[i - 1] = 0.0 if denom == 0 else numer / denom

    if L_s > 1:
        step = (1.0 - 0.25) / (L_s - 1)
        weights = [0.25 + j * step for j in range(L_s)]
    else:
        weights = [0.25]
    division_score = [d * w for d, w in zip(division_score, weights)]

    best_idx = 0
    best_val = -math.inf
    for i, v in enumerate(division_score):
        if v > best_val:
            best_val = v
            best_idx = i
    new_end = best_idx + 1 + kmer

    hits = sum(1 for km in pot_kmers[:new_end] if km in rep_kmers)
    if hits / new_end > 0.15:
        return new_end
    return 0


def _clamp(x: int, lo: int, hi: int) -> int:
    if x < lo:
        return lo
    if x > hi:
        return hi
    return x


def handle_edge_repeat(
    rows: list[RepeatRow],
    sequence_vector: str,
    sequence_vector_start: int,
    template_sequence: str = "",
) -> list[RepeatRow]:
    """Port `src/repeats/handle_edge_repeat.R`.

    Adjusts the length of each edge repeat (array boundaries and internal
    gap boundaries) to maximise kmer-alignment against the representative,
    and appends new repeats just outside the adjusted edge when there's
    enough similarity.

    Mirrors R's per-branch arithmetic quirks — notably the R code indexes
    `sequence_vector` with *absolute* positions (no `- svs + 1` shift)
    when scoring newly-added repeats' `score` column, but uses relative
    indexing for their `score_template`. Both are reproduced verbatim.
    """
    potential_seq_width = 1.2
    if len(rows) < 2:
        return [dict(r) for r in rows]

    svs = sequence_vector_start + 1  # R: `sequence_vector_start <- sequence_vector_start + 1`

    table: list[RepeatRow] = []
    for r in rows:
        nr = dict(r)
        nr["start"] = int(nr["start"])
        nr["end"] = int(nr["end"])
        nr["score"] = float(nr["score"])
        nr["eval"] = float(nr["eval"])
        nr["score_template"] = float(nr["score_template"])
        nr["width"] = int(nr["width"])
        table.append(nr)

    rep_fw = table[0]["representative"]
    rep_len = len(rep_fw)
    rep_rev = rep_fw[::-1]
    rep_comp = rep_fw.translate(_COMP_ONLY)
    rep_rev_comp = rep_comp[::-1]

    mean_rep_size = sum(int(r["width"]) for r in table) / len(table)
    window = round(mean_rep_size * potential_seq_width)
    L = len(sequence_vector)

    left_edges: list[int] = [table[0]["start"]]
    left_edges_rep_id: list[int] = [0]
    right_edges: list[int] = []
    right_edges_rep_id: list[int] = []
    for i in range(1, len(table)):
        if (table[i]["start"] - table[i - 1]["end"]) > mean_rep_size:
            left_edges.append(table[i]["start"])
            left_edges_rep_id.append(i)
            right_edges.append(table[i - 1]["end"])
            right_edges_rep_id.append(i - 1)
    right_edges.append(table[-1]["end"])
    right_edges_rep_id.append(len(table) - 1)

    if len(left_edges) != len(right_edges):
        raise RuntimeError("handle_edge_repeat: wrong edge definition")

    def abs_slice(start: int, end: int) -> str:
        """1-based inclusive slice into sequence_vector."""
        if start > end or start < 1 or end < 1:
            return ""
        return sequence_vector[start - 1:end]

    def score_row(rep: str, seq: str) -> float:
        if rep == "" or seq == "":
            return 0.0
        return _adist(rep, seq) / len(rep) * 100

    for i_e in range(len(left_edges)):
        lrid = left_edges_rep_id[i_e]
        rrid = right_edges_rep_id[i_e]

        # --- Left edge -------------------------------------------------
        strand_l = table[lrid]["strand"]
        pe1 = table[lrid]["end"] - svs + 1
        ps1 = pe1 - window
        ps1 = _clamp(ps1, 1, max(L, 1))
        pe1 = _clamp(pe1, 1, max(L, 1))
        sp = abs_slice(ps1, pe1)[::-1]
        rep_for_find = rep_rev if strand_l == "+" else rep_comp
        new_offset = find_edge_best_start_end(sp, rep_for_find)
        table[lrid]["start"] = table[lrid]["end"] - new_offset
        table[lrid]["eval"] = -1.0
        subseq = abs_slice(table[lrid]["start"] - svs + 1, pe1)
        table[lrid]["score"] = score_row(table[lrid]["representative"], subseq)
        if template_sequence:
            table[lrid]["score_template"] = score_row(template_sequence, subseq)
        left_edges[i_e] = table[lrid]["start"]

        if (table[lrid]["end"] - table[lrid]["start"] + 1) >= (mean_rep_size / 2):
            ps3 = table[lrid]["start"] - svs + 1 - window
            pe3 = table[lrid]["start"] - svs
            ps3 = _clamp(ps3, 1, max(L, 1))
            pe3 = _clamp(pe3, 1, max(L, 1))
            sp3 = abs_slice(ps3, pe3)[::-1]
            best_new_start = find_edge_best_start_end(sp3, rep_for_find)
            if best_new_start > 0:
                new_row = dict(table[lrid])
                new_row["start"] = table[lrid]["start"] - best_new_start - 1
                new_row["end"] = table[lrid]["start"] - 1
                # R line 79: absolute-position indexing (apparent bug — matched for parity).
                sub_new_abs = sequence_vector[max(new_row["start"], 1) - 1:new_row["end"]] if new_row["start"] <= L else ""
                new_row["score"] = score_row(table[lrid]["representative"], sub_new_abs)
                if template_sequence:
                    sub_new_rel = abs_slice(new_row["start"] - svs + 1, new_row["end"] - svs + 1)
                    new_row["score_template"] = score_row(template_sequence, sub_new_rel)
                table.append(new_row)
                left_edges[i_e] = new_row["start"]

        # --- Right edge ------------------------------------------------
        strand_r = table[rrid]["strand"]
        ps2 = table[rrid]["start"] - svs + 1
        pe2 = table[rrid]["start"] - svs + window
        ps2 = _clamp(ps2, 1, max(L, 1))
        pe2 = _clamp(pe2, 1, max(L, 1))
        sp2 = abs_slice(ps2, pe2)
        if strand_r == "+":
            rep_for_find_r = rep_fw
        else:
            rep_for_find_r = rep_rev_comp  # rev(array_representative_rev)
        new_offset_r = find_edge_best_start_end(sp2, rep_for_find_r)
        table[rrid]["end"] = table[rrid]["start"] + new_offset_r
        table[rrid]["eval"] = -1.0
        subseq_r = abs_slice(ps2, table[rrid]["end"] - svs + 1)
        if strand_r == "+":
            table[rrid]["score"] = score_row(table[rrid]["representative"], subseq_r)
            if template_sequence:
                table[rrid]["score_template"] = score_row(template_sequence, subseq_r)
        else:
            # R uses rev_comp_string(rep) for the minus strand score.
            table[rrid]["score"] = (
                _adist(rev_comp_string(table[rrid]["representative"]), subseq_r)
                / len(table[rrid]["representative"]) * 100 if subseq_r else 0.0
            )
            if template_sequence:
                table[rrid]["score_template"] = (
                    _adist(rev_comp_string(template_sequence), subseq_r)
                    / len(template_sequence) * 100 if subseq_r else 0.0
                )
        right_edges[i_e] = table[rrid]["end"]

        if (table[rrid]["end"] - table[rrid]["start"] + 1) >= (mean_rep_size / 2):
            ps4 = table[rrid]["end"] - svs + 2
            pe4 = table[rrid]["end"] - svs + 1 + window
            ps4 = _clamp(ps4, 1, max(L, 1))
            pe4 = _clamp(pe4, 1, max(L, 1))
            sp4 = abs_slice(ps4, pe4)
            best_new_end = find_edge_best_start_end(sp4, rep_for_find_r)
            # R line 140 uses `> 0` for plus, line 171 uses `> 1` for minus —
            # the asymmetry is kept.
            min_ok = 0 if strand_r == "+" else 1
            if best_new_end > min_ok:
                new_row = dict(table[rrid])
                new_row["start"] = table[rrid]["end"] + 1
                new_row["end"] = table[rrid]["end"] + best_new_end + 1
                # R lines 144 / 175: score uses *relative* indexing on the
                # right side (unlike the left-side bug on line 79). Matches
                # R's inconsistency verbatim.
                sub_r_new = abs_slice(new_row["start"] - svs + 1, new_row["end"] - svs + 1)
                if strand_r == "+":
                    new_row["score"] = score_row(table[rrid]["representative"], sub_r_new)
                    if template_sequence:
                        new_row["score_template"] = score_row(template_sequence, sub_r_new)
                else:
                    new_row["score"] = (
                        _adist(rev_comp_string(table[rrid]["representative"]), sub_r_new)
                        / len(table[rrid]["representative"]) * 100 if sub_r_new else 0.0
                    )
                    if template_sequence:
                        new_row["score_template"] = (
                            _adist(rev_comp_string(template_sequence), sub_r_new)
                            / len(template_sequence) * 100 if sub_r_new else 0.0
                        )
                table.append(new_row)
                right_edges[i_e] = new_row["end"]

    for r in table:
        r["width"] = int(r["end"]) - int(r["start"]) + 1
    table = [r for r in table if r["width"] > 2]

    # R re-selects the final column order (line 201) — done here too.
    col_order = ["seqID", "arrayID", "start", "end", "strand", "score", "eval", "width", "class", "representative", "score_template"]
    return [{k: r[k] for k in col_order} for r in table]
