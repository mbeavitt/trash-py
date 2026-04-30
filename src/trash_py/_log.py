"""Tiny logging helper. Stdout, optional quiet, sectioned output."""
from __future__ import annotations

import subprocess
import sys
import time
from dataclasses import dataclass


_QUIET = False
_WORKER = False
_T0 = time.monotonic()


@dataclass
class _ToolStat:
    count: int = 0
    # Wall-clock span across all calls (possibly across workers): the time
    # the user actually waited, not summed CPU. Stored as `time.time()`
    # values so they're comparable across processes.
    first_start: float = float("inf")
    last_end: float = float("-inf")


_EXTERNAL: dict[str, _ToolStat] = {}
_ANNOUNCED: set[str] = set()


def configure(*, quiet: bool = False) -> None:
    """Reset elapsed-time origin, external-tool counters, and quiet mode."""
    global _QUIET, _T0
    _QUIET = quiet
    _T0 = time.monotonic()
    _EXTERNAL.clear()
    _ANNOUNCED.clear()


def set_worker_mode(on: bool) -> None:
    """In worker mode, `run_external` skips the one-shot `running {tool}...`
    line so parallel workers don't all print it. Stats are still recorded
    locally in the worker process and drained via `pop_stats`."""
    global _WORKER
    _WORKER = on


def announce_tool(tool: str) -> None:
    """Emit the one-shot `running {tool}...` line from the parent process
    before submitting parallel work. No-op if already announced this round."""
    if tool in _ANNOUNCED:
        return
    _ANNOUNCED.add(tool)
    detail(f"running {tool}...")


def pop_stats() -> dict[str, _ToolStat]:
    """Snapshot and clear the per-process external-tool counters. Called by
    a worker before returning its result so the parent can aggregate."""
    snap = dict(_EXTERNAL)
    _EXTERNAL.clear()
    return snap


def merge_stats(snapshot: dict[str, _ToolStat]) -> None:
    """Add a worker's per-tool counts/times into the parent's `_EXTERNAL`."""
    for tool, stat in snapshot.items():
        local = _EXTERNAL.setdefault(tool, _ToolStat())
        local.count += stat.count
        local.first_start = min(local.first_start, stat.first_start)
        local.last_end = max(local.last_end, stat.last_end)


def run_external(tool: str, *args, **kwargs) -> subprocess.CompletedProcess:
    """Wrap `subprocess.run` and record per-tool wall-time + call count.

    Emits a one-time `running {tool}...` line on the first call since the
    last summary, so the user knows we've entered that subprocess phase.
    Suppressed in worker mode — the parent announces before dispatch.
    """
    stat = _EXTERNAL.setdefault(tool, _ToolStat())
    if stat.count == 0 and not _WORKER and tool not in _ANNOUNCED:
        _ANNOUNCED.add(tool)
        detail(f"running {tool}...")
    t0 = time.time()
    try:
        return subprocess.run(*args, **kwargs)
    finally:
        t1 = time.time()
        stat.count += 1
        stat.first_start = min(stat.first_start, t0)
        stat.last_end = max(stat.last_end, t1)


def tool_summary(tool: str) -> None:
    """Emit a detail line summarising calls/time recorded for `tool` since the
    previous call, then clear the counter. No-op if no calls were recorded
    since the last summary."""
    tool_summary_group(tool)


def tool_summary_group(*tools: str) -> None:
    """Combine summaries for tools dispatched together (e.g. both running
    inside one parallel block) into a single line: per-tool call counts,
    then one wall-clock span covering all calls of all listed tools."""
    items: list[tuple[str, _ToolStat]] = []
    first_start = float("inf")
    last_end = float("-inf")
    for tool in tools:
        stat = _EXTERNAL.pop(tool, None)
        _ANNOUNCED.discard(tool)
        if stat is None or stat.count == 0:
            continue
        items.append((tool, stat))
        first_start = min(first_start, stat.first_start)
        last_end = max(last_end, stat.last_end)
    if not items:
        return
    parts = [
        f"{tool}: {stat.count:,} call{'' if stat.count == 1 else 's'}"
        for tool, stat in items
    ]
    wall = max(0.0, last_end - first_start)
    detail(f"{', '.join(parts)}, {format_elapsed(wall)}")


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
