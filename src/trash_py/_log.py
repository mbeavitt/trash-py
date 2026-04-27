"""Tiny logging helper. Stdout, optional quiet, sectioned output."""
from __future__ import annotations

import sys
import time


_QUIET = False
_T0 = time.monotonic()


def configure(*, quiet: bool = False) -> None:
    """Reset elapsed-time origin and set quiet mode."""
    global _QUIET, _T0
    _QUIET = quiet
    _T0 = time.monotonic()


def _emit(line: str = "") -> None:
    if _QUIET:
        return
    print(line, file=sys.stdout, flush=True)


def header(text: str) -> None:
    """Top-of-output banner. No leading blank line, no indent."""
    _emit(text)


def info(msg: str = "") -> None:
    """A line of body text. Plain, no indent."""
    _emit(msg)


def detail(msg: str) -> None:
    """An indented detail line under a section."""
    _emit(f"  {msg}")


def section(title: str) -> None:
    """Section heading — blank line, then `==> title`, then blank line."""
    _emit("")
    _emit(f"==> {title}")


def elapsed_marker() -> None:
    """Indented bracketed cumulative-elapsed marker, used at end of a section."""
    _emit(f"  [{format_elapsed(elapsed_since_start())}]")


def warn(msg: str) -> None:
    print(f"warning: {msg}", file=sys.stderr)


def format_elapsed(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f}s"
    if seconds < 3600:
        m = int(seconds // 60)
        s = int(seconds - 60 * m)
        return f"{m}m{s:02d}s"
    h = int(seconds // 3600)
    m = int((seconds - 3600 * h) // 60)
    return f"{h}h{m:02d}m"


def elapsed_since_start() -> float:
    return time.monotonic() - _T0
