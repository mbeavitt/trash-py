"""HOR line plots — a best-effort matplotlib rendering of the upstream R plot.

The original R plots were acknowledged to be rough; we don't aim for pixel
parity, just a readable dot-plot of HOR blocks coloured by divergence, with the
same axis-compression (empty regions squeezed out) the R version used.
"""
from __future__ import annotations

from pathlib import Path


def _lwd(ax_len: float) -> float:
    for thresh, lw in ((1e7, 1), (5e6, 2), (2.5e6, 3), (1.25e6, 4)):
        if ax_len >= thresh:
            return lw
    return 5


def plot_hors(annotated: list[list], start_adjusted: list[float], bins,
              out_png: Path, chrA: str, hor_class: str, threshold: int) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    unit = 1_000_000.0
    ax_len = sum(bins.sizes)
    if ax_len <= 0:
        return

    # Adjust HOR bp coords into the compressed axis (same as R).
    def adj(v: float) -> float:
        for bs, be, corr in zip(bins.starts, bins.ends, bins.corrections):
            if bs <= v < be:
                return v - corr
        return v

    snv = [row[15] for row in annotated]
    smax = max(snv) if snv else 1.0
    red_at = max(threshold * 10, 1e-9)
    order = sorted(range(len(annotated)), key=lambda i: snv[i], reverse=True)

    fig, ax = plt.subplots(figsize=(10, 10), dpi=110)
    segs, colors = [], []
    cmap = _green_yellow_red()
    for i in order:
        row = annotated[i]
        xa, xb = adj(row[6]) / unit, adj(row[8]) / unit   # start.A.bp, end.A.bp
        ya, yb = adj(row[7]) / unit, adj(row[9]) / unit   # start.B.bp, end.B.bp
        segs.append([(xa, ya), (xb, yb)])
        colors.append(cmap(min(snv[i] / red_at, 1.0)))
    ax.add_collection(LineCollection(segs, colors=colors, linewidths=_lwd(ax_len)))

    xs = [s / unit for s in start_adjusted]
    ax.plot(xs, [0] * len(xs), ",", color="#22222290")
    ax.plot([0] * len(xs), xs, ",", color="#22222290")

    lim = ax_len / unit
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel(f"{chrA}, Mbp")
    ax.set_ylabel(f"{chrA}, Mbp")
    ax.set_title(f"HORs {hor_class} {chrA}")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def _green_yellow_red():
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("gyr", ["green", "yellow", "red"])
