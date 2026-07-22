#!/usr/bin/env python
"""Plot a single-chromosome HOR table as a self-similarity dot-plot.

Takes a `HORs_<class>_<seq>.csv` file (the table written by `trash-py hor`/`HORT.R`) and
renders a dot-plot-style view of higher order repeat pairs

Usage:
    python scripts/plot_hor_table.py HORs_178_1_Chr3.csv
    python scripts/plot_hor_table.py HORs_178_1_Chr3.csv -o chr3.png --label Chr3 -t 25
"""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize

# Dark "stained-glass" palette.
BG, INK, GRID = "#0d0d12", "#c8c8d0", "#2a2a33"


def read_table(path: Path):
    """Read a HOR table, returning each HOR's block-A and block-B centres (in
    Mbp), its divergence (SNV/kbp), and the sequence name."""
    xa, yb, snv = [], [], []
    seq_id = ""
    with open(path, newline="") as f:
        for row in csv.DictReader(f):
            xa.append((float(row["start.A.bp"]) + float(row["end.A.bp"])) / 2e6)
            yb.append((float(row["start.B.bp"]) + float(row["end.B.bp"])) / 2e6)
            snv.append(float(row["SNV_per_kbp"]))
            seq_id = seq_id or row["chrA"]
    return np.array(xa), np.array(yb), np.array(snv), seq_id


def plot_dotplot(xa, yb, snv, out_path: Path, label: str, hor_class: str, threshold: int):
    """One translucent dot per HOR at its paired-block centres (block A on x,
    block B on y), mirrored across y = x so the plot is the symmetric
    self-comparison. Tandem HORs sit near the diagonal, inversions on the
    anti-diagonal; colour is divergence."""
    xs = np.concatenate([xa, yb])
    ys = np.concatenate([yb, xa])

    # Conserved HORs (low divergence) glow bright; divergent ones recede.
    vmax = max(threshold * 10, 1e-9)
    shade = np.clip(np.concatenate([snv, snv]) / vmax, 0.0, 1.0)
    colours = plt.get_cmap("magma")(1.0 - shade)
    colours[:, 3] = 0.5  # alpha, so overlapping dots blend

    lo, hi = min(xs.min(), ys.min()), max(xs.max(), ys.max())
    pad = 0.02 * (hi - lo) or 1.0
    lo, hi = lo - pad, hi + pad
    size = 1.0 if len(xa) > 200_000 else 3.0 if len(xa) > 20_000 else 8.0

    fig, ax = plt.subplots(figsize=(9, 9), dpi=260)
    fig.patch.set_facecolor(BG)
    ax.set_facecolor(BG)
    ax.plot([lo, hi], [lo, hi], color=GRID, lw=0.8, zorder=0)  # main diagonal
    ax.scatter(xs, ys, s=size, c=colours, linewidths=0, rasterized=True)
    ax.set(xlim=(lo, hi), ylim=(lo, hi), aspect="equal",
           xlabel=f"HOR block A position — {label} (Mbp)",
           ylabel=f"HOR block B position — {label} (Mbp)",
           title=f"HOR self-similarity — {hor_class}, {label}")

    for spine in ax.spines.values():
        spine.set_color(GRID)
    ax.tick_params(colors=INK, labelsize=8)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.title.set_color(INK)
    ax.grid(True, color=GRID, lw=0.4, alpha=0.5)

    bar = fig.colorbar(plt.cm.ScalarMappable(Normalize(0, vmax), "magma_r"),
                       ax=ax, fraction=0.046, pad=0.02)
    bar.set_label("SNV per kbp", color=INK)
    bar.ax.tick_params(colors=INK, labelsize=8)
    bar.outline.set_edgecolor(GRID)

    fig.tight_layout()
    fig.savefig(out_path, facecolor=BG)
    plt.close(fig)
    return out_path


def class_from_filename(name: str, seq_id: str) -> str:
    """Recover the repeat class from a `HORs_<class>_<seq>.csv` filename."""
    stem = name[:-4] if name.endswith(".csv") else name
    stem = stem[len("HORs_"):] if stem.startswith("HORs_") else stem
    if seq_id and stem.endswith(f"_{seq_id}"):
        stem = stem[: -len(seq_id) - 1]
    return stem or "?"


def main() -> None:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("table", type=Path, help="HOR table CSV (HORs_<class>_<seq>.csv)")
    p.add_argument("-o", "--output", type=Path,
                   help="output PNG (default: alongside the table)")
    p.add_argument("--label",
                   help="sequence name for the axes/title (default: chrA from the table)")
    p.add_argument("-c", "--class", dest="hor_class",
                   help="repeat class for the title (default: parsed from the filename)")
    p.add_argument("-t", "--threshold", type=int, default=4,
                   help="HOR threshold used, which sets the colour scale (default 4)")
    args = p.parse_args()

    xa, yb, snv, seq_id = read_table(args.table)
    if len(xa) == 0:
        p.error(f"no HOR rows in {args.table}")

    label = args.label or seq_id or args.table.stem
    hor_class = args.hor_class or class_from_filename(args.table.name, seq_id)
    out = args.output or args.table.with_name(f"HORs_dotplot_{hor_class}_{seq_id or label}.png")

    plot_dotplot(xa, yb, snv, out, label, hor_class, args.threshold)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
