"""Final-stage tabulation: per-array aggregates, repeat-table cleanup,
and (optional) sequence-column extraction for the repeats table.
"""
from __future__ import annotations

import math
import statistics
from typing import Any

from .sequence import rev_comp_string


ARRAYS_COLUMNS = [
    "start", "end", "seqID", "numID", "score", "top_N", "top_5_N",
    "representative", "class", "array_num_ID",
    "repeats_number", "median_repeat_width", "median_score",
]

REPEATS_COLUMNS = [
    "seqID", "arrayID", "start", "end", "strand", "score", "eval",
    "width", "class", "score_template",
]

REPEATS_WITH_SEQ_COLUMNS = REPEATS_COLUMNS + ["sequence"]


def summarise_arrays(
    arrays: list[dict[str, Any]],
    repeats: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    """Reassign each array's representative from its first repeat, append
    aggregate columns (`repeats_number`, `median_repeat_width`,
    `median_score`), and drop arrays that ended up with zero repeats.

    Note: matches an upstream R quirk where the lookup key is the array's
    1-based row index rather than `array_num_ID`. After stage 8's
    sort+renumber the two are equal in practice, so this is consistent
    on real inputs but is the literal behaviour to match for parity.
    """
    out = [dict(a) for a in arrays]
    first_repeat_by_one_based: dict[int, dict] = {}
    for r in repeats:
        key = int(r["arrayID"])
        if key not in first_repeat_by_one_based:
            first_repeat_by_one_based[key] = r

    for i, arr in enumerate(out):
        anid = int(arr["array_num_ID"])
        has_any = any(int(r["arrayID"]) == anid for r in repeats)
        if has_any:
            lookup = first_repeat_by_one_based.get(i + 1)
            if lookup is not None:
                arr["representative"] = lookup["representative"]

    for arr in out:
        anid = int(arr["array_num_ID"])
        widths = [int(r["width"]) for r in repeats if int(r["arrayID"]) == anid]
        scores = [float(r["score"]) for r in repeats if int(r["arrayID"]) == anid]
        arr["repeats_number"] = len(widths)
        arr["median_repeat_width"] = math.ceil(statistics.median(widths)) if widths else 0
        arr["median_score"] = math.ceil(statistics.median(scores)) if scores else -1

    return [a for a in out if a["repeats_number"] != 0]


def strip_repeats_columns(repeats: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Drop the `representative` column before the repeats table is written."""
    return [{k: r[k] for k in REPEATS_COLUMNS} for r in repeats]


def append_sequence_column(
    repeats: list[dict[str, Any]],
    fasta_by_seqID: dict[str, str],
) -> list[dict[str, Any]]:
    """Append the genomic-sequence slice for each repeat row (rev-complemented
    on the minus strand). Returns a fresh list."""
    out: list[dict[str, Any]] = []
    for r in repeats:
        start = int(r["start"])
        end = int(r["end"])
        seq = fasta_by_seqID[r["seqID"]][start - 1:end]
        if r["strand"] == "-":
            seq = rev_comp_string(seq)
        out.append({**r, "sequence": seq})
    return out
