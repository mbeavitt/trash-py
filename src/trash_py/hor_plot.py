"""HOR self-similarity dot-plot — a "stained-glass" rendering.

One translucent dot per matched repeat *unit*: within each HOR the paired
monomers (A_u ↔ B_u) become individual points, so a HOR reads as a run of dots
along a diagonal — tandem HORs on the main diagonal, inverted-repeat HORs on the
anti-diagonal — and the plane is mirrored across y=x for the symmetric self-plot.
Dots are coloured by divergence — a *magnitude*, so a single perceptual
sequential ramp (magma), not the old green→yellow→red rainbow — with alpha so
dense regions blend like stained glass on a dark ground.

Best-effort: we don't aim for pixel parity with the original R plot.
"""
from __future__ import annotations

from pathlib import Path

# Dark surface + recessive ink (kept local so the module is self-contained).
_BG = "#0d0d12"
_INK = "#c8c8d0"
_GRID = "#2a2a33"


def plot_hors(annotated: list[list], out_png: Path, chrA: str, hor_class: str,
              threshold: int, repeats, unit_val: float = 1_000_000.0) -> None:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if not annotated:
        return

    # Genomic midpoint (bp) of each repeat, indexed by its 1-based unit number.
    mid = np.array([(r.start + r.end) / 2.0 for r in repeats], dtype=float)

    a = np.asarray(
        [(r[0], r[1], r[2], r[3], r[4], r[15]) for r in annotated], dtype=np.int64
    )
    start_A, end_A, start_B, end_B, direction = (a[:, i] for i in range(5))
    snv = a[:, 5].astype(float)
    ks = (end_A - start_A + 1).astype(np.int64)     # units per HOR

    # Ragged expansion: one row per (HOR, unit). `off` is the within-HOR unit u.
    n = int(ks.sum())
    off = np.arange(n) - np.repeat(ks.cumsum() - ks, ks)
    rA = np.repeat(start_A, ks) + off                # A monomer index (1-based)
    dirr = np.repeat(direction, ks)
    # parallel pairs A_u with B_u (start_B+u); antiparallel with B reversed (end_B-u)
    rB = np.where(dirr == 1, np.repeat(start_B, ks) + off, np.repeat(end_B, ks) - off)
    snv = np.repeat(snv, ks)

    # Satellite arrays can yield tens of millions of unit-pairs (Chr1 178_1:
    # ~70M) — far too many to scatter. Subsample to a representative cloud.
    # (A proper fix is a 2D density raster; see TODO in the branch.)
    cap = 1_500_000
    if n > cap:
        rng = np.random.default_rng(0)
        pick = rng.choice(n, size=cap, replace=False)
        rA, rB, snv = rA[pick], rB[pick], snv[pick]
        n = cap

    xa = mid[rA - 1] / unit_val
    yb = mid[rB - 1] / unit_val
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
    size = 0.5 if n > 1_000_000 else (1.5 if n > 100_000 else 6.0)

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
