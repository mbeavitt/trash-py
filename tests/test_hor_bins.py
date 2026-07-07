"""Tests for the R-`pretty()`/`hist()` port used to compute plot bins and the
`start.adjusted` column. Every expected value was cross-checked against R
(`pretty()` and `calc_plot_regs.R`) directly.
"""
from __future__ import annotations

from trash_py.hor_bins import r_pretty, calc_plot_regs


def test_r_pretty_matches_r():
    # Verified against R: pretty(c(0,100), n=1) and pretty(c(3,137), n=1)
    assert r_pretty(0, 100, 1) == [0.0, 100.0]
    assert r_pretty(3, 137, 1) == [0.0, 100.0, 200.0]


def test_calc_plot_regs_single_merged_region():
    # starts spread across three touching 100k bins -> merged into one region.
    # Verified against R's calc_plot_regs.
    starts = [10, 20, 30, 120_000, 130_000, 240_000]
    bins = calc_plot_regs(starts, 100_000)
    assert bins.start_adjusted == [0.0]
    assert bins.end_adjusted == [300_000.0]
    assert bins.corrections == [0.0]
    assert bins.sizes == [300_000.0]


def test_calc_plot_regs_compression():
    # Three clusters separated by empty gaps: the gaps are squeezed out and the
    # cumulative correction grows. Verified against R's calc_plot_regs.
    starts = [12_600_050] * 5 + [12_950_000] * 3 + [15_800_000] * 10
    bins = calc_plot_regs(starts, 100_000)
    assert bins.start_adjusted == [0.0, 100_000.0, 200_000.0]
    assert bins.corrections == [12_600_000.0, 12_800_000.0, 15_500_000.0]
    assert bins.sizes == [100_000.0, 100_000.0, 100_000.0]
