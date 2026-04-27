"""GFF writer matching upstream TRASH output.

`export_gff` writes tab-separated records with a CR (`\\r`) record
separator — that's what the R reference produces, and downstream
tooling depends on the exact bytes.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any

from .summarise import ARRAYS_COLUMNS, REPEATS_COLUMNS


def _format_attr_value(value: Any) -> str:
    if isinstance(value, float) and value.is_integer():
        return str(int(value))
    return str(value)


def export_gff(
    rows: list[dict[str, Any]],
    path: Path,
    *,
    seqid: str | int,
    source: str,
    type_: str,
    start: str | int,
    end: str | int,
    score: str | int = ".",
    strand: str | int = ".",
    phase: str = ".",
    attributes: list[int] | str = ".",
    attribute_names: list[str] | None = None,
) -> None:
    """Write `rows` to `path` as GFF3.

    Each parameter is either a literal string (used verbatim on every
    row) or a 1-based column index (looked up per-row from
    `ARRAYS_COLUMNS` or `REPEATS_COLUMNS` depending on row shape).
    """
    attribute_names = attribute_names or []
    ordered_keys = ARRAYS_COLUMNS if "top_N" in rows[0] else REPEATS_COLUMNS

    def _col(spec: str | int, default: str = ".") -> list[str]:
        if isinstance(spec, int):
            key = ordered_keys[spec - 1]
            return [_format_attr_value(r[key]) for r in rows]
        return [str(spec)] * len(rows)

    seqid_col = _col(seqid)
    source_col = _col(source)
    type_col = _col(type_)
    start_col = _col(start)
    end_col = _col(end)
    score_col = _col(score)
    strand_col = _col(strand)
    phase_col = _col(phase)

    if isinstance(attributes, list):
        attr_col: list[str] = []
        for r in rows:
            parts = [
                f"{name}{_format_attr_value(r[ordered_keys[col_idx - 1]])}"
                for name, col_idx in zip(attribute_names, attributes)
            ]
            attr_col.append(";".join(parts))
    else:
        attr_col = [str(attributes)] * len(rows)

    with Path(path).open("w", newline="") as f:
        for i in range(len(rows)):
            fields = [
                seqid_col[i], source_col[i], type_col[i],
                start_col[i], end_col[i], score_col[i],
                strand_col[i], phase_col[i], attr_col[i],
            ]
            f.write("\t".join(fields))
            f.write("\r")
