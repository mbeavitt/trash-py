"""Tiny logging helper. Stdout, optional quiet, sectioned output."""
from __future__ import annotations

import subprocess
import sys
import time
from dataclasses import dataclass


_QUIET = False
_T0 = time.monotonic()


@dataclass
class _ToolStat:
    count: int = 0
    total: float = 0.0


_EXTERNAL: dict[str, _ToolStat] = {}


def configure(*, quiet: bool = False) -> None:
    """Reset elapsed-time origin, external-tool counters, and quiet mode."""
    global _QUIET, _T0
    _QUIET = quiet
    _T0 = time.monotonic()
    _EXTERNAL.clear()


def run_external(tool: str, *args, **kwargs) -> subprocess.CompletedProcess:
    """Wrap `subprocess.run` and record per-tool wall-time + call count.

    Emits a one-time `running {tool}...` line on the first call since the
    last summary, so the user knows we've entered that subprocess phase.
    """
    stat = _EXTERNAL.setdefault(tool, _ToolStat())
    if stat.count == 0:
        detail(f"running {tool}...")
    t0 = time.monotonic()
    try:
        return subprocess.run(*args, **kwargs)
    finally:
        stat.count += 1
        stat.total += time.monotonic() - t0


def tool_summary(tool: str) -> None:
    """Emit a detail line summarising calls/time recorded for `tool` since the
    previous call, then clear the counter. No-op if no calls were recorded
    since the last summary."""
    stat = _EXTERNAL.pop(tool, None)
    if stat is None or stat.count == 0:
        return
    plural = "" if stat.count == 1 else "s"
    detail(f"{tool}: {stat.count:,} call{plural}, {format_elapsed(stat.total)}")


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
