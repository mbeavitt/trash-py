"""HOR self-similarity dot-plot — a "stained-glass" rendering.

One translucent dot per HOR, placed at the centre of its paired blocks (block A
on x, block B on y) and mirrored across y=x, so tandem HORs trace the main
diagonal and inverted-repeat HORs the anti-diagonal. Dots are coloured by
divergence — a *magnitude*, so a single perceptual sequential ramp (magma), not
the old green→yellow→red rainbow — with alpha so dense regions blend like stained
glass on a dark ground.

Best-effort: we don't aim for pixel parity with the original R plot.
"""
from __future__ import annotations

from pathlib import Path

# Dark surface + recessive ink (kept local so the module is self-contained).
_BG = "#0d0d12"
_INK = "#c8c8d0"
_GRID = "#2a2a33"


def plot_hors(annotated: list[list], out_png: Path, chrA: str, hor_class: str,
              threshold: int, unit_val: float = 1_000_000.0) -> None:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not annotated:
        return

    a = np.asarray(
        [(r[6], r[7], r[8], r[9], r[15]) for r in annotated], dtype=float
    )
    # One dot per HOR at the centre of its paired blocks (A on x, B on y), plus
    # the mirror across y = x for the symmetric self-dot-plot.
    xa = (a[:, 0] + a[:, 2]) / 2.0 / unit_val   # block A centre
    yb = (a[:, 1] + a[:, 3]) / 2.0 / unit_val   # block B centre
    snv = a[:, 4]
    xs = np.concatenate([xa, yb])
    ys = np.concatenate([yb, xa])

    # Colour by divergence (low = bright/hot, high = dim), magma on dark ground.
    vmax = max(threshold * 10, 1e-9)
    t = np.clip(np.concatenate([snv, snv]) / vmax, 0.0, 1.0)
    colors = plt.get_cmap("magma")(1.0 - t)     # conserved HORs glow, divergent recede
    colors[:, 3] = 0.5                           # alpha -> stained-glass blend

    # Zoom to the occupied region (HORs only form inside the satellite array, so
    # most of the chromosome is empty) — a square window with a small margin.
    lo, hi = float(min(xs.min(), ys.min())), float(max(xs.max(), ys.max()))
    pad = 0.02 * (hi - lo) if hi > lo else 1.0
    lo, hi = lo - pad, hi + pad
    size = 1.0 if len(annotated) > 200_000 else (3.0 if len(annotated) > 20_000 else 8.0)

    fig, ax = plt.subplots(figsize=(9, 9), dpi=130)
    fig.patch.set_facecolor(_BG)
    ax.set_facecolor(_BG)

    ax.plot([lo, hi], [lo, hi], color=_GRID, lw=0.8, zorder=0)  # main diagonal
    sc = ax.scatter(xs, ys, s=size, c=colors, marker="o", linewidths=0)
    sc.set_rasterized(True)                     # keep the PNG small/fast with ~10^6 dots
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_aspect("equal")
    ax.set_xlabel(f"{chrA}, Mbp", color=_INK)
    ax.set_ylabel(f"{chrA}, Mbp", color=_INK)
    ax.set_title(f"HOR self-similarity — {hor_class}, {chrA}", color=_INK)
    for spine in ax.spines.values():
        spine.set_color(_GRID)
    ax.tick_params(colors=_INK, labelsize=8)
    ax.grid(True, color=_GRID, lw=0.4, alpha=0.5)

    # Colorbar keyed to divergence.
    import matplotlib as mpl
    sm = mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(0, vmax), cmap="magma_r")
    cb = fig.colorbar(sm, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label("SNV per kbp", color=_INK)
    cb.ax.yaxis.set_tick_params(color=_INK, labelsize=8)
    plt.setp(cb.ax.get_yticklabels(), color=_INK)
    cb.outline.set_edgecolor(_GRID)

    fig.tight_layout()
    fig.savefig(out_png, facecolor=_BG)
    plt.close(fig)
