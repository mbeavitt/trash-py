"""Resumable pipeline state for runs that outlive a scheduler's time limit.

A large genome can take longer than the wall-clock allocation an HPC
queue will give it (CSD3's 12 h cap, say). This module lets a run

* stream its per-stage results to a checkpoint directory as they are
  produced,
* notice that the job is about to be killed — either a signal the
  scheduler sends ahead of the limit (`sbatch --signal=B:USR1@600`) or a
  deadline read from `SLURM_JOB_END_TIME`/`--time-limit` — and stop at
  the next task boundary with everything finished so far on disk,
* pick up where the previous run stopped when the same command is run
  again against the same checkpoint directory.

The unit of progress is one *item*: a sequence in the window-scoring
stage, a region in array identification, an array in repeat mapping.
Items are recorded strictly in order, so a stage's state is just "the
first N items are done, and here are the rows they produced".

A run that stops to be resumed exits `EXIT_CHECKPOINTED` (75,
`EX_TEMPFAIL`) so a batch script can tell "requeue me" from a real
failure.
"""
from __future__ import annotations

import json
import os
import re
import signal
import time
from pathlib import Path
from typing import Any, Iterable

from . import _log as log


# Bumped when the on-disk layout changes in a way older state can't satisfy.
STATE_VERSION = 1

# Exit code for "stopped early, state saved, run me again". 75 is
# EX_TEMPFAIL from sysexits.h: a temporary failure, retry later.
EXIT_CHECKPOINTED = 75

DEFAULT_INTERVAL = 300.0
DEFAULT_MARGIN = 600.0
DEFAULT_SIGNALS = "USR1,TERM"

STATE_FILE = "state.json"


class Checkpointed(Exception):
    """Raised once state is safely on disk, to unwind out of the pipeline.

    Caught in `cli.main`, which turns it into `EXIT_CHECKPOINTED`.
    """

    def __init__(self, reason: str, stage: str, done: int, total: int, directory: Path):
        self.reason = reason
        self.stage = stage
        self.done = done
        self.total = total
        self.directory = directory
        super().__init__(
            f"stopped during {stage} ({done:,}/{total:,} done) after {reason}"
        )


_UNITS = {"d": 86400.0, "h": 3600.0, "m": 60.0, "s": 1.0}
_COMPOUND = re.compile(
    r"^(?:(\d+(?:\.\d+)?)d)?(?:(\d+(?:\.\d+)?)h)?"
    r"(?:(\d+(?:\.\d+)?)m)?(?:(\d+(?:\.\d+)?)s)?$"
)


def parse_duration(text: str) -> float:
    """Seconds from `11h`, `1h30m`, `90m`, `[dd-]hh:mm:ss`, `mm:ss`, or a
    bare number of seconds. Raises ValueError on anything else."""
    raw = str(text).strip().lower()
    if not raw:
        raise ValueError("empty duration")

    if ":" in raw:
        days = 0.0
        rest = raw
        if "-" in rest:
            day_part, rest = rest.split("-", 1)
            days = float(day_part)
        try:
            parts = [float(p) for p in rest.split(":")]
        except ValueError:
            raise ValueError(f"bad duration: {text}") from None
        if len(parts) == 2:
            hours, minutes, seconds = 0.0, parts[0], parts[1]
        elif len(parts) == 3:
            hours, minutes, seconds = parts
        else:
            raise ValueError(f"bad duration: {text}")
        return days * 86400.0 + hours * 3600.0 + minutes * 60.0 + seconds

    try:
        return float(raw)
    except ValueError:
        pass

    match = _COMPOUND.match(raw)
    if not match or not any(match.groups()):
        raise ValueError(f"bad duration: {text}")
    total = 0.0
    for value, unit in zip(match.groups(), "dhms"):
        if value is not None:
            total += float(value) * _UNITS[unit]
    return total


def parse_signals(spec: str) -> list[signal.Signals]:
    """`"USR1,TERM"` -> [SIGUSR1, SIGTERM]. Unknown names raise ValueError."""
    out: list[signal.Signals] = []
    for name in str(spec).split(","):
        name = name.strip().upper()
        if not name:
            continue
        if not name.startswith("SIG"):
            name = f"SIG{name}"
        try:
            sig = getattr(signal, name)
        except AttributeError:
            raise ValueError(f"unknown signal: {name}") from None
        if not isinstance(sig, signal.Signals):
            raise ValueError(f"unknown signal: {name}")
        if sig not in out:
            out.append(sig)
    return out


def resolve_deadline(
    time_limit: str | None, margin: float, now: float | None = None
) -> tuple[float | None, str]:
    """Absolute unix time at which to stop, and where it came from.

    `--time-limit` is a budget measured from now; otherwise Slurm's
    `SLURM_JOB_END_TIME` (unix seconds, exported into every job) gives the
    real end of the allocation. `margin` is subtracted from either, so the
    run stops early enough to save state and exit cleanly.
    """
    now = time.time() if now is None else now
    if time_limit:
        return now + parse_duration(time_limit) - margin, "--time-limit"
    env_end = os.environ.get("SLURM_JOB_END_TIME", "").strip()
    if env_end:
        try:
            return float(env_end) - margin, "SLURM_JOB_END_TIME"
        except ValueError:
            log.warn(f"ignoring unparseable SLURM_JOB_END_TIME={env_end!r}")
    return None, ""


class Stage:
    """One resumable pipeline stage.

    `resume_index` is how many items a previous run already finished, so
    a stage restarts by skipping that many; `rows` holds every row the
    stage has produced, replayed ones first.
    """

    def __init__(
        self,
        checkpoint: "Checkpoint",
        name: str,
        total: int,
        resume_index: int = 0,
        rows: list | None = None,
    ) -> None:
        self._checkpoint = checkpoint
        self.name = name
        self.total = total
        self.resume_index = resume_index
        self.rows: list = [] if rows is None else rows
        self.items_done = resume_index
        self.done = False

    @property
    def pending(self) -> int:
        return max(0, self.total - self.resume_index)

    def skip(self, tasks: list) -> list:
        """The slice of `tasks` this run still has to do."""
        return tasks[self.resume_index:]

    def record(self, rows: Iterable) -> None:
        """Record one item's output rows. May raise `Checkpointed`."""
        rows = list(rows)
        self.rows.extend(rows)
        self.items_done += 1
        self._checkpoint._after_record(self, rows)

    def finish(self) -> None:
        self.done = True
        self._checkpoint._after_stage(self)

    def __enter__(self) -> "Stage":
        return self

    def __exit__(self, exc_type, exc, tb) -> bool:
        if exc_type is None:
            self.finish()
        return False


class _NullCheckpoint:
    """Stand-in used when checkpointing is off: stages just accumulate."""

    enabled = False
    directory: Path | None = None
    signal_numbers: tuple[int, ...] = ()

    def stage(self, name: str, total: int) -> Stage:
        return Stage(self, name, total)

    def install_signal_handlers(self) -> None:  # pragma: no cover - trivial
        pass

    def describe(self) -> None:  # pragma: no cover - trivial
        pass

    def complete(self) -> None:  # pragma: no cover - trivial
        pass

    def close(self) -> None:  # pragma: no cover - trivial
        pass

    def _after_record(self, stage: Stage, rows: list) -> None:
        pass

    def _after_stage(self, stage: Stage) -> None:
        pass


NO_CHECKPOINT = _NullCheckpoint()


class Checkpoint:
    """Streaming, resumable state in `directory`.

    One `<stage>.jsonl` per stage holds the rows, appended as they are
    produced; `state.json` records how many items and rows of each file
    are committed. `state.json` is written last and atomically, so it is
    the source of truth — trailing rows from a run that died mid-flight
    are truncated away on load.
    """

    enabled = True

    def __init__(
        self,
        directory: Path,
        fingerprint: dict[str, Any],
        *,
        interval: float = DEFAULT_INTERVAL,
        deadline: float | None = None,
        deadline_source: str = "",
        signals: Iterable[signal.Signals] = (),
        restart: bool = False,
        keep: bool = False,
    ) -> None:
        self.directory = Path(directory)
        self.fingerprint = fingerprint
        self.interval = interval
        self.deadline = deadline
        self.deadline_source = deadline_source
        self.signals = list(signals)
        self.keep = keep
        self._stop_reason: str | None = None
        self._last_save = time.monotonic()
        self._files: dict[str, Any] = {}
        self._stages: dict[str, Stage] = {}
        self._resumed = False

        self._guard_directory()
        if restart and self.directory.exists():
            self._wipe()
        self.directory.mkdir(parents=True, exist_ok=True)
        self._state = self._load_state()

    # -- setup ---------------------------------------------------------

    def _guard_directory(self) -> None:
        """Refuse a directory that holds anything but our own state — we
        delete this directory when the run finishes, and `--checkpoint out/`
        is an easy thing to type."""
        if not self.directory.is_dir():
            return
        strays = [
            entry.name for entry in self.directory.iterdir()
            if entry.name != STATE_FILE
            and entry.name != STATE_FILE + ".tmp"
            and entry.suffix != ".jsonl"
        ]
        if strays:
            raise SystemExit(
                f"refusing to use {self.directory} as a checkpoint directory: "
                f"it holds files that are not checkpoint state "
                f"({', '.join(sorted(strays)[:3])}).\n"
                f"Point --checkpoint at a directory of its own."
            )

    def _wipe(self) -> None:
        """Remove our own files, then the directory itself. Deliberately not
        `rmtree`: anything we did not write stays put."""
        for entry in sorted(self.directory.iterdir()):
            if entry.name in (STATE_FILE, STATE_FILE + ".tmp") or entry.suffix == ".jsonl":
                entry.unlink()
        try:
            self.directory.rmdir()
        except OSError as exc:  # pragma: no cover - guarded against above
            log.warn(f"left {self.directory} in place: {exc}")

    def _state_path(self) -> Path:
        return self.directory / STATE_FILE

    def _load_state(self) -> dict[str, Any]:
        path = self._state_path()
        if not path.exists():
            return {
                "version": STATE_VERSION,
                "fingerprint": self.fingerprint,
                "stages": {},
            }
        try:
            state = json.loads(path.read_text())
        except (OSError, ValueError) as exc:
            raise SystemExit(
                f"checkpoint at {self.directory} is unreadable ({exc}).\n"
                f"Delete it, or re-run with --restart to start over."
            ) from None
        if state.get("version") != STATE_VERSION:
            raise SystemExit(
                f"checkpoint at {self.directory} was written by an "
                f"incompatible version of trash-py.\n"
                f"Re-run with --restart to start over."
            )
        stored = state.get("fingerprint", {})
        if stored != self.fingerprint:
            differing = sorted(
                k for k in set(stored) | set(self.fingerprint)
                if stored.get(k) != self.fingerprint.get(k)
            )
            raise SystemExit(
                f"checkpoint at {self.directory} was written for a different "
                f"run (differs in: {', '.join(differing)}).\n"
                f"Point --checkpoint somewhere else, or re-run with --restart."
            )
        self._resumed = bool(state.get("stages"))
        return state

    @property
    def signal_numbers(self) -> tuple[int, ...]:
        """The signals this run treats as "stop and save" — worker processes
        ignore exactly these, so the parent can drain the pool and save."""
        return tuple(int(sig) for sig in self.signals)

    def install_signal_handlers(self) -> None:
        """Ask for a graceful stop when the scheduler warns us."""
        for sig in self.signals:
            try:
                signal.signal(sig, self._handle_signal)
            except (OSError, ValueError) as exc:  # pragma: no cover - platform
                log.warn(f"cannot handle {sig.name}: {exc}")

    def _handle_signal(self, signum: int, frame: Any) -> None:
        # Async-signal context: only set a flag. The stop happens at the
        # next item boundary, where it is safe to write files.
        if self._stop_reason is None:
            try:
                name = signal.Signals(signum).name
            except ValueError:  # pragma: no cover - defensive
                name = str(signum)
            self._stop_reason = f"SIG{name.removeprefix('SIG')}"

    def describe(self) -> None:
        log.detail(f"checkpoint: {self.directory}/" + (" (resuming)" if self._resumed else ""))
        if self.deadline is not None:
            left = max(0.0, self.deadline - time.time())
            log.detail(
                f"stopping in {log.format_elapsed(left)} if unfinished "
                f"(from {self.deadline_source})"
            )
        if self.signals:
            log.detail(
                "stop signals: " + ", ".join(s.name for s in self.signals)
            )

    # -- stages --------------------------------------------------------

    def stage(self, name: str, total: int) -> Stage:
        entry = self._state["stages"].get(name, {"items": 0, "rows": 0, "done": False})
        items = int(entry.get("items", 0))
        n_rows = int(entry.get("rows", 0))
        if items > total:
            raise SystemExit(
                f"checkpoint at {self.directory} has {items:,} {name} items but "
                f"this run only has {total:,}; re-run with --restart."
            )
        rows = self._replay(name, n_rows)
        stage = Stage(self, name, total, resume_index=items, rows=rows)
        stage.done = bool(entry.get("done", False)) and items >= total
        self._stages[name] = stage
        if items:
            log.detail(
                f"resuming {name}: {items:,}/{total:,} already done "
                f"({len(rows):,} rows restored)"
            )
        return stage

    def _path(self, name: str) -> Path:
        return self.directory / f"{name}.jsonl"

    def _replay(self, name: str, n_rows: int) -> list:
        """Read the first `n_rows` committed rows and drop anything after
        them (a previous run may have died between two saves)."""
        path = self._path(name)
        rows: list = []
        if path.exists():
            # readline() rather than iteration: we need an exact byte offset
            # to truncate at, and the file iterator's read-ahead hides it.
            with path.open("r+", encoding="utf-8", newline="\n") as fh:
                offset = 0
                while len(rows) < n_rows:
                    line = fh.readline()
                    if not line:
                        break
                    offset += len(line.encode("utf-8"))
                    stripped = line.strip()
                    if stripped:
                        rows.append(json.loads(stripped))
                fh.seek(offset)
                fh.truncate()
        if len(rows) < n_rows:
            raise SystemExit(
                f"checkpoint file {path} is short ({len(rows):,} of "
                f"{n_rows:,} rows); re-run with --restart."
            )
        self._files[name] = self._path(name).open("a", encoding="utf-8")
        return rows

    def _after_record(self, stage: Stage, rows: list) -> None:
        fh = self._files.get(stage.name)
        if fh is not None:
            for row in rows:
                fh.write(json.dumps(row, separators=(",", ":")) + "\n")

        if self._stop_reason is None and self._past_deadline():
            self._stop_reason = "time limit"

        if self._stop_reason is not None:
            self.save(stage_note=stage)
            log.info()
            log.warn(
                f"stopping early ({self._stop_reason}): "
                f"{stage.name} {stage.items_done:,}/{stage.total:,} saved to "
                f"{self.directory}"
            )
            raise Checkpointed(
                self._stop_reason, stage.name, stage.items_done, stage.total,
                self.directory,
            )

        if time.monotonic() - self._last_save >= self.interval:
            self.save(stage_note=stage)

    def _past_deadline(self) -> bool:
        return self.deadline is not None and time.time() >= self.deadline

    def _after_stage(self, stage: Stage) -> None:
        self.save(stage_note=None)

    # -- persistence ---------------------------------------------------

    def save(self, stage_note: Stage | None = None) -> None:
        """Flush every stage file, then atomically publish `state.json`."""
        for fh in self._files.values():
            fh.flush()
            os.fsync(fh.fileno())

        stages = dict(self._state.get("stages", {}))
        for name, stage in self._stages.items():
            stages[name] = {
                "items": stage.items_done,
                "rows": len(stage.rows),
                "total": stage.total,
                "done": stage.done,
            }
        self._state["stages"] = stages
        self._state["version"] = STATE_VERSION
        self._state["fingerprint"] = self.fingerprint
        self._state["updated"] = time.time()

        tmp = self.directory / (STATE_FILE + ".tmp")
        tmp.write_text(json.dumps(self._state, indent=1))
        os.replace(tmp, self._state_path())
        self._last_save = time.monotonic()

        if stage_note is not None and self._stop_reason is None:
            log.detail(
                f"checkpoint saved: {stage_note.name} "
                f"{stage_note.items_done:,}/{stage_note.total:,}"
            )

    def complete(self) -> None:
        """The whole run finished: drop the state unless asked to keep it."""
        self.close()
        if self.keep:
            log.detail(f"checkpoint kept: {self.directory}/")
            return
        if self.directory.exists():
            self._wipe()

    def close(self) -> None:
        for fh in self._files.values():
            try:
                fh.close()
            except OSError:  # pragma: no cover - defensive
                pass
        self._files.clear()


def fingerprint_for(args: Any, stem: str) -> dict[str, Any]:
    """Identity of a run: resuming is only safe if these all match."""
    from . import __version__

    def file_stamp(path: Any) -> list[Any]:
        p = Path(path)
        try:
            st = p.stat()
        except OSError:
            return [str(p), None, None]
        return [str(p.resolve()), st.st_size, int(st.st_mtime)]

    fp: dict[str, Any] = {
        "trash_py": __version__,
        "state": STATE_VERSION,
        "fasta": file_stamp(args.fasta),
        "max_rep_size": int(args.max_rep_size),
        "min_rep_size": int(args.min_rep_size),
        "name": stem,
    }
    templates = getattr(args, "templates", None)
    fp["templates"] = file_stamp(templates) if templates is not None else None
    return fp


def from_args(args: Any, stem: str) -> _NullCheckpoint | Checkpoint:
    """Build the checkpoint described by the CLI flags (or the null one)."""
    directory = getattr(args, "checkpoint", None)
    if not directory:
        return NO_CHECKPOINT

    margin = getattr(args, "time_margin", DEFAULT_MARGIN)
    deadline, source = resolve_deadline(
        getattr(args, "time_limit", None), margin
    )
    return Checkpoint(
        Path(directory),
        fingerprint_for(args, stem),
        interval=float(getattr(args, "checkpoint_interval", DEFAULT_INTERVAL)),
        deadline=deadline,
        deadline_source=source,
        signals=parse_signals(getattr(args, "checkpoint_signal", DEFAULT_SIGNALS)),
        restart=bool(getattr(args, "restart", False)),
        keep=bool(getattr(args, "keep_checkpoint", False)),
    )
