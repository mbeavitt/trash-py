"""Port of `src/auxi/genomic_bins_starts.R`."""
from __future__ import annotations

import math


def genomic_bins_starts(start: int = 1, end: int = 0, bin_number: int = 0, bin_size: int = 0) -> list[int]:
    if bin_number > 0 and bin_size > 0:
        raise ValueError("Use either bin_number or bin_size, not both")
    if bin_number == 0 and bin_size == 0:
        raise ValueError("Use either bin_number or bin_size")
    if end < start:
        raise ValueError("end < start")
    if bin_size >= (end - start):
        return [start]

    if bin_number > 0:
        # Dead branch in current call sites; preserved for parity.
        raise NotImplementedError("bin_number path not needed by current stages")

    # R: seq(start, end - bin_size, bin_size) — inclusive both ends, step=bin_size.
    last = end - bin_size
    n = (last - start) // bin_size + 1
    start_positions = [start + i * bin_size for i in range(n)]
    # Drop the tail if the final window would be less than half a bin wide.
    if (end - start_positions[-1]) < (bin_size / 2):
        start_positions = start_positions[:-1]
    return start_positions
