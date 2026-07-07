"""Plot-region binning for the HOR module.

Faithful port of `calc_plot_regs.R`, which leans on R's `hist(..., breaks = n)`
— and therefore on R's `pretty()` algorithm — to choose bin edges. We reproduce
`pretty()` (R's `src/appl/pretty.c::R_pretty`) closely enough to match R's bin
edges, so the `start.adjusted` column of `repeats_with_hors_*.csv` lines up with
the reference output.
"""
from __future__ import annotations

import math
import sys
from dataclasses import dataclass

_DBL_EPSILON = sys.float_info.epsilon        # 2.220446049250313e-16
_DBL_MIN = sys.float_info.min                # 2.2250738585072014e-308
_DBL_MAX = sys.float_info.max
_ROUNDING_EPS = 1e-10


def r_pretty(lo: float, hi: float, ndiv: int, min_n: int = 1,
             shrink_sml: float = 0.75, h: float = 1.5, h5: float = 2.75) -> list[float]:
    """Return break points the way R's ``pretty(c(lo, hi), n = ndiv, min.n)`` does."""
    lo = float(lo)
    hi = float(hi)
    dx = hi - lo
    if dx == 0 and hi == 0:
        cell = 1.0
        i_small = True
    else:
        cell = max(abs(lo), abs(hi))
        U = 1 + (1.0 / (1.0 + h) if (h5 >= 1.5 * h + 0.5) else 1.5 / (1.0 + h5))
        i_small = dx < cell * U * max(1, ndiv) * _DBL_EPSILON * 3

    if i_small:
        if cell > 10:
            cell = 9 + cell / 10
        cell *= shrink_sml
        if min_n > 1:
            cell /= min_n
    else:
        cell = dx
        if ndiv > 1:
            cell /= ndiv

    if cell < 20 * _DBL_MIN:
        cell = 20 * _DBL_MIN
    elif cell * 10 > _DBL_MAX:
        cell = 0.1 * _DBL_MAX

    base = 10.0 ** math.floor(math.log10(cell))
    unit = base
    if (2 * base) - cell < h * (cell - unit):
        unit = 2 * base
        if (5 * base) - cell < h5 * (cell - unit):
            unit = 5 * base
            if (10 * base) - cell < h * (cell - unit):
                unit = 10 * base

    ns = math.floor(lo / unit + _ROUNDING_EPS)
    nu = math.ceil(hi / unit - _ROUNDING_EPS)
    while ns * unit > lo + _ROUNDING_EPS * unit:
        ns -= 1
    while nu * unit < hi - _ROUNDING_EPS * unit:
        nu += 1

    k = int(0.5 + nu - ns)
    if k < min_n:
        k = min_n - k
        if ns >= 0:
            nu += k // 2
            ns -= k // 2 + k % 2
        else:
            ns -= k // 2
            nu += k // 2 + k % 2
        ndiv_out = min_n
    else:
        ndiv_out = k

    lo_out = ns * unit
    hi_out = nu * unit
    n = ndiv_out + 1
    # seq.int(lo_out, hi_out, length.out = n)
    if n == 1:
        return [lo_out]
    step = (hi_out - lo_out) / (n - 1)
    return [lo_out + i * step for i in range(n)]


@dataclass
class Bins:
    starts: list[float]
    ends: list[float]
    counts: list[int]
    sizes: list[float]
    corrections: list[float]
    start_adjusted: list[float]
    end_adjusted: list[float]


def _bin_counts(x: list[int], breaks: list[float]) -> list[int]:
    """Mirror R's ``hist`` counting: right-closed intervals, first bin closed on
    the left (``right = TRUE, include.lowest = TRUE``)."""
    nb = len(breaks) - 1
    counts = [0] * nb
    lo0 = breaks[0]
    for v in x:
        # binary search for right-closed bin
        left, right = 0, nb
        while left < right:
            mid = (left + right) // 2
            if v <= breaks[mid + 1]:
                right = mid
            else:
                left = mid + 1
        b = left
        # include.lowest: value equal to the very first edge belongs to bin 0
        if v == lo0:
            b = 0
        counts[b] += 1
    return counts


def calc_plot_regs(starts: list[int], window: int) -> Bins:
    """Port of ``calc_plot_regs``: pretty bins, drop empties, merge touching
    bins, then compute the cumulative left-shift (`correction`) that squeezes
    out the empty gaps for plotting."""
    smin, smax = min(starts), max(starts)
    nbreaks = math.ceil((smax - smin) / window)
    breaks = r_pretty(smin, smax, nbreaks, min_n=1)
    counts = _bin_counts(starts, breaks)

    b_starts = list(breaks[:-1])
    b_ends = list(breaks[1:])
    rows = [[bs, be, bc] for bs, be, bc in zip(b_starts, b_ends, counts) if bc != 0]

    # merge bins whose edges touch (ends[i] == starts[i+1])
    i = 0
    while i < len(rows) - 1:
        if rows[i][1] == rows[i + 1][0]:
            rows[i][1] = rows[i + 1][1]
            rows[i][2] += rows[i + 1][2]
            del rows[i + 1]
            i -= 1
        i += 1

    starts_f = [r[0] for r in rows]
    ends_f = [r[1] for r in rows]
    counts_f = [r[2] for r in rows]
    sizes = [e - s for s, e in zip(starts_f, ends_f)]

    n = len(rows)
    dist_to_next = [0.0] * n
    correction = [starts_f[0]] * n
    i = 0
    while i < n - 1:
        dist_to_next[i] = starts_f[i + 1] - ends_f[i]
        correction[i + 1] = correction[i + 1] + sum(dist_to_next[: i + 1])
        i += 1

    start_adj = [s - c for s, c in zip(starts_f, correction)]
    end_adj = [e - c for e, c in zip(ends_f, correction)]
    return Bins(starts_f, ends_f, counts_f, sizes, correction, start_adj, end_adj)
