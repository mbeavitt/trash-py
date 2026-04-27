"""DNA-string utilities used across stages."""
from __future__ import annotations


_COMP_TABLE = str.maketrans("acgtn", "tgcan")


def rev_comp_string(seq: str) -> str:
    """Reverse-complement a lowercase DNA string. Unknown bases pass through."""
    return seq.translate(_COMP_TABLE)[::-1]
