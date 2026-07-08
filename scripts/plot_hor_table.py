#!/usr/bin/env python
"""Plot a single-chromosome HOR table as a self-similarity dot-plot.

Takes a `HORs_<class>_<seq>.csv` file (the table written by `trash-py hor`) and
renders the same dot-plot the tool produces — handy for re-plotting or tweaking
without re-running detection.

Usage:
    python scripts/plot_hor_table.py HORs_178_1_Chr3.csv
    python scripts/plot_hor_table.py HORs_178_1_Chr3.csv -o chr3.png --label Chr3 -t 25
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

# Run straight from a source checkout without installing.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))
from trash_py.hor_plot import plot_hors  # noqa: E402

# Column layout of the HOR table (after dropping R's leading row-name column).
# 0 start_A 1 end_A 2 start_B 3 end_B 4 direction 5 total_variant
# 6 start.A.bp 7 start.B.bp 8 end.A.bp 9 end.B.bp 10 chrA 11 chrB
# 12 block.size.in.units 13 block.A.size.bp 14 block.B.size.bp 15 SNV_per_kbp
_STR_COLS = {10, 11}  # chrA, chrB stay as strings; everything else is numeric


def load_table(path: Path) -> list[list]:
    rows: list[list] = []
    with path.open(newline="") as f:
        reader = csv.reader(f)
        next(reader, None)  # header
        for raw in reader:
            cells = raw[1:]  # drop R's row-name column
            if len(cells) < 16:
                continue
            rows.append([c if i in _STR_COLS else float(c)
                         for i, c in enumerate(cells[:16])])
    return rows


def derive_class(filename: str, seq_id: str) -> str:
    """Recover the repeat class from a `HORs_<class>_<seq>.csv` name."""
    stem = filename[:-4] if filename.endswith(".csv") else filename
    if stem.startswith("HORs_"):
        stem = stem[len("HORs_"):]
    if seq_id and stem.endswith(f"_{seq_id}"):
        stem = stem[: -(len(seq_id) + 1)]
    return stem or "?"


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("table", type=Path, help="HOR table CSV (HORs_<class>_<seq>.csv)")
    p.add_argument("-o", "--output", type=Path, default=None,
                   help="output PNG (default: alongside the table)")
    p.add_argument("--label", default=None,
                   help="axis/title name for the sequence (default: chrA from the table)")
    p.add_argument("-c", "--class", dest="hor_class", default=None,
                   help="repeat class for the title (default: parsed from the filename)")
    p.add_argument("-t", "--threshold", type=int, default=25,
                   help="HOR divergence threshold used (sets the colour scale; default 25)")
    a = p.parse_args(argv)

    if not a.table.exists():
        print(f"table not found: {a.table}", file=sys.stderr)
        return 2
    rows = load_table(a.table)
    if not rows:
        print("no HOR rows in table", file=sys.stderr)
        return 1

    seq_id = str(rows[0][10])
    label = a.label or seq_id
    hor_class = a.hor_class or derive_class(a.table.name, seq_id)
    out = a.output or a.table.with_name(f"HORs_dotplot_{hor_class}_{seq_id}.png")

    plot_hors(rows, out, label, hor_class, a.threshold)
    print(f"wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
