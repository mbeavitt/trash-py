"""CSV writer matching R's `write.csv(x, row.names=FALSE)` byte-for-byte.

Rules observed on macOS R 4.4.3:
* Headers are wrapped in double-quotes.
* Character values are wrapped in double-quotes.
* Numeric values are unquoted, formatted with 15 significant digits
  (`%.15g` — R's default when `getOption("digits")` is 15, which is
  what R's internal `as.character.numeric()` uses for CSV export).
* Integer-valued floats are written without a decimal point
  (e.g., `2223`, not `2223.0`).
* Fields separated by `,`. Rows terminated by `\\n`. File ends with `\\n`.
"""
from __future__ import annotations

from pathlib import Path
from typing import Any, Sequence


def _format(value: Any) -> str:
    if isinstance(value, bool):
        return '"TRUE"' if value else '"FALSE"'
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if value.is_integer():
            return str(int(value))
        return f"{value:.15g}"
    return f'"{value}"'


def write_csv_r_style(path: Path, columns: Sequence[str], rows: Sequence[Sequence[Any]]) -> None:
    with Path(path).open("w", newline="") as f:
        f.write(",".join(f'"{c}"' for c in columns))
        f.write("\n")
        for row in rows:
            f.write(",".join(_format(v) for v in row))
            f.write("\n")
