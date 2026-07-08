"""HOR self-similarity dot-plot — a "stained-glass" rendering.

Each HOR block is one translucent diagonal segment on a chromosome-vs-chromosome
plane: parallel (tandem) HORs lie along the main diagonal, antiparallel
(inverted) HORs along the anti-diagonal, and the plane is mirrored across y=x so
the plot is the familiar symmetric self-dot-plot. Segments are coloured by
divergence — a *magnitude*, so a single perceptual sequential ramp (magma), not
the old green→yellow→red rainbow — and drawn with alpha so overlapping panes
blend like stained glass on a dark ground.

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
    from matplotlib.collections import LineCollection

    if not annotated:
        return

    a = np.asarray(
        [(r[6], r[7], r[8], r[9], r[4], r[15]) for r in annotated], dtype=float
    )
    sA, sB, eA, eB, direction, snv = (a[:, 0], a[:, 1], a[:, 2], a[:, 3],
                                      a[:, 4], a[:, 5])
    sA, sB, eA, eB = (x / unit_val for x in (sA, sB, eA, eB))

    # Parallel HORs run start->end on both axes; antiparallel are drawn along the
    # anti-diagonal (flip the B ends) so inverted repeats read as the other slope.
    para = direction == 1
    p0x, p0y = sA, np.where(para, sB, eB)
    p1x, p1y = eA, np.where(para, eB, sB)

    # Segments + their mirror across y = x -> symmetric self-dot-plot.
    seg = np.stack([np.column_stack([p0x, p0y]), np.column_stack([p1x, p1y])], axis=1)
    seg_m = seg[:, :, ::-1]
    segments = np.concatenate([seg, seg_m], axis=0)

    # Colour by divergence (low = bright/hot, high = dim), magma on dark ground.
    vmax = max(threshold * 10, 1e-9)
    t = np.clip(np.concatenate([snv, snv]) / vmax, 0.0, 1.0)
    cmap = plt.get_cmap("magma")
    colors = cmap(1.0 - t)                 # conserved HORs glow, divergent recede
    colors[:, 3] = 0.55                     # alpha -> stained-glass blend

    # Zoom to the occupied region (HORs only form inside the satellite array, so
    # most of the chromosome is empty) — a square window with a small margin.
    lo = float(min(sA.min(), sB.min(), eA.min(), eB.min()))
    hi = float(max(sA.max(), sB.max(), eA.max(), eB.max()))
    pad = 0.02 * (hi - lo) if hi > lo else 1.0
    lo, hi = lo - pad, hi + pad
    lw = 0.6 if len(annotated) > 200_000 else (1.0 if len(annotated) > 20_000 else 1.6)

    fig, ax = plt.subplots(figsize=(9, 9), dpi=130)
    fig.patch.set_facecolor(_BG)
    ax.set_facecolor(_BG)

    lc = LineCollection(segments, colors=colors, linewidths=lw, capstyle="round")
    lc.set_rasterized(True)                 # keep the PNG small/fast with ~10^6 panes
    ax.add_collection(lc)

    ax.plot([lo, hi], [lo, hi], color=_GRID, lw=0.8, zorder=0)  # main diagonal
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
