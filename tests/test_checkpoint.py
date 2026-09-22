"""Checkpoint/resume: a run stopped part-way must produce exactly the
output an uninterrupted run produces.

The stop is forced with `--time-limit 0`, which puts the deadline in the
past: every run then stops after committing its first item, so driving the
CLI in a loop walks the pipeline forward one item at a time through every
stage boundary.
"""
from __future__ import annotations

import json
import signal
import subprocess
import sys
from pathlib import Path

import pytest

from trash_py import checkpoint as ckpt_mod


REPO_ROOT = Path(__file__).resolve().parent.parent
SMALL_FASTA = REPO_ROOT / "tests" / "data" / "ath_Chr1_extraction_trc.fasta"

OUTPUT_FILES = [
    "regarrays.csv", "aregarrays.csv", "classarrays.csv",
    "no_repeats_arrays.csv", "arrays.csv", "arrays.gff",
    "repeats.csv", "repeats.gff", "repeats_with_seq.csv",
]


def run_cli(*extra: str, output: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, "-m", "trash_py", "-f", str(SMALL_FASTA),
         "-o", str(output), *extra],
        capture_output=True,
        text=True,
    )


@pytest.fixture(scope="module")
def uninterrupted(tmp_path_factory) -> Path:
    out = tmp_path_factory.mktemp("checkpoint_reference")
    result = run_cli(output=out)
    assert result.returncode == 0, result.stderr
    return out


def assert_same_output(produced: Path, reference: Path) -> None:
    for name in OUTPUT_FILES:
        fname = f"{SMALL_FASTA.name}_{name}"
        assert (produced / fname).read_bytes() == (reference / fname).read_bytes(), (
            f"{fname} differs from the uninterrupted run"
        )


def drive_to_completion(
    out: Path, ckpt: Path, *extra: str, max_runs: int = 40
) -> int:
    """Re-run the same command until it completes; returns the run count."""
    for run in range(1, max_runs + 1):
        result = run_cli(
            "--checkpoint", str(ckpt), "--time-limit", "0", *extra, output=out
        )
        if result.returncode == 0:
            return run
        assert result.returncode == ckpt_mod.EXIT_CHECKPOINTED, (
            f"run {run} failed:\n{result.stderr}"
        )
        assert (ckpt / "state.json").exists()
    pytest.fail(f"still unfinished after {max_runs} runs")


@pytest.mark.parametrize("processes", ["1", "4"])
def test_resumed_run_matches_uninterrupted(
    tmp_path: Path, uninterrupted: Path, processes: str
) -> None:
    out = tmp_path / "out"
    runs = drive_to_completion(out, tmp_path / "ckpt", "-p", processes)
    assert runs > 1, "the run was never actually interrupted"
    assert_same_output(out, uninterrupted)


def test_checkpoint_is_removed_once_the_run_finishes(tmp_path: Path) -> None:
    ckpt = tmp_path / "ckpt"
    drive_to_completion(tmp_path / "out", ckpt)
    assert not ckpt.exists()


def test_keep_checkpoint_leaves_the_state_behind(tmp_path: Path) -> None:
    ckpt = tmp_path / "ckpt"
    drive_to_completion(tmp_path / "out", ckpt, "--keep-checkpoint")
    state = json.loads((ckpt / "state.json").read_text())
    assert all(stage["done"] for stage in state["stages"].values())


def test_default_checkpoint_directory_sits_under_the_output_dir(
    tmp_path: Path,
) -> None:
    out = tmp_path / "out"
    result = run_cli("--checkpoint", "--time-limit", "0", "--keep-checkpoint",
                     output=out)
    assert result.returncode == ckpt_mod.EXIT_CHECKPOINTED
    assert (out / f"{SMALL_FASTA.name}.checkpoint" / "state.json").exists()


def test_resume_refuses_a_checkpoint_from_a_different_run(tmp_path: Path) -> None:
    out = tmp_path / "out"
    ckpt = tmp_path / "ckpt"
    first = run_cli("--checkpoint", str(ckpt), "--time-limit", "0", output=out)
    assert first.returncode == ckpt_mod.EXIT_CHECKPOINTED

    # Same checkpoint, different repeat-size parameters: results from the
    # first run are not valid input for the second.
    clash = run_cli("--checkpoint", str(ckpt), "-m", "500", output=out)
    assert clash.returncode != 0
    assert "--restart" in clash.stderr
    assert "max_rep_size" in clash.stderr

    # --restart throws the stale state away and runs the whole thing.
    fixed = run_cli("--checkpoint", str(ckpt), "-m", "500", "--restart",
                    output=out)
    assert fixed.returncode == 0, fixed.stderr
    assert "resuming" not in fixed.stdout


def test_stdin_input_cannot_be_checkpointed(tmp_path: Path) -> None:
    result = subprocess.run(
        [sys.executable, "-m", "trash_py", "-f", "-", "-o", str(tmp_path / "out"),
         "--checkpoint", str(tmp_path / "ckpt")],
        input=SMALL_FASTA.read_text(),
        capture_output=True,
        text=True,
    )
    assert result.returncode == 2
    assert "--checkpoint" in result.stderr


def test_a_signal_stops_the_run_at_the_next_item(tmp_path: Path) -> None:
    ckpt = ckpt_mod.Checkpoint(tmp_path / "ckpt", {"fp": 1})
    stage = ckpt.stage("things", total=10)

    stage.record([{"a": 1}])            # no stop requested yet
    assert stage.items_done == 1

    ckpt._handle_signal(signal.SIGUSR1, None)
    with pytest.raises(ckpt_mod.Checkpointed) as stop:
        stage.record([{"a": 2}])
    assert stop.value.reason == "SIGUSR1"
    assert stop.value.stage == "things"

    # Both items — the one before the signal and the one that tripped the
    # stop — are committed, and a fresh Checkpoint sees them.
    ckpt.close()
    resumed = ckpt_mod.Checkpoint(tmp_path / "ckpt", {"fp": 1})
    stage = resumed.stage("things", total=10)
    assert stage.resume_index == 2
    assert stage.rows == [{"a": 1}, {"a": 2}]


def test_rows_written_after_the_last_save_are_discarded(tmp_path: Path) -> None:
    """A run killed outright (SIGKILL, node failure) can leave rows past the
    last committed point. state.json is the source of truth, so they go."""
    ckpt = ckpt_mod.Checkpoint(tmp_path / "ckpt", {"fp": 1})
    stage = ckpt.stage("things", total=10)
    stage.record([{"a": 1}])
    stage.finish()
    ckpt.close()

    with (tmp_path / "ckpt" / "things.jsonl").open("a") as fh:
        fh.write(json.dumps({"a": "uncommitted"}) + "\n")

    resumed = ckpt_mod.Checkpoint(tmp_path / "ckpt", {"fp": 1})
    stage = resumed.stage("things", total=10)
    assert stage.rows == [{"a": 1}]
    assert (tmp_path / "ckpt" / "things.jsonl").read_text().count("\n") == 1


@pytest.mark.parametrize(
    "text,seconds",
    [
        ("11h", 39600.0),
        ("1h30m", 5400.0),
        ("690m", 41400.0),
        ("11:30:00", 41400.0),
        ("30:00", 1800.0),
        ("1-00:00:00", 86400.0),
        ("600", 600.0),
        ("600s", 600.0),
    ],
)
def test_parse_duration(text: str, seconds: float) -> None:
    assert ckpt_mod.parse_duration(text) == seconds


@pytest.mark.parametrize("text", ["", "soon", "12x", "1:2:3:4"])
def test_parse_duration_rejects_nonsense(text: str) -> None:
    with pytest.raises(ValueError):
        ckpt_mod.parse_duration(text)


def test_parse_signals() -> None:
    assert ckpt_mod.parse_signals("USR1,TERM") == [signal.SIGUSR1, signal.SIGTERM]
    assert ckpt_mod.parse_signals("SIGUSR2") == [signal.SIGUSR2]
    with pytest.raises(ValueError):
        ckpt_mod.parse_signals("NOTASIGNAL")


def test_deadline_comes_from_slurm_when_no_limit_is_given(monkeypatch) -> None:
    monkeypatch.setenv("SLURM_JOB_END_TIME", "1000000")
    deadline, source = ckpt_mod.resolve_deadline(None, margin=600.0)
    assert (deadline, source) == (999400.0, "SLURM_JOB_END_TIME")

    # An explicit budget wins, and is measured from now.
    deadline, source = ckpt_mod.resolve_deadline("2h", margin=600.0, now=0.0)
    assert (deadline, source) == (6600.0, "--time-limit")

    monkeypatch.delenv("SLURM_JOB_END_TIME")
    assert ckpt_mod.resolve_deadline(None, margin=600.0) == (None, "")


def test_a_directory_with_other_files_is_refused(tmp_path: Path) -> None:
    """The checkpoint directory is deleted on success, so it must never be
    a directory holding anything else."""
    ckpt = tmp_path / "not-a-checkpoint"
    ckpt.mkdir()
    (ckpt / "precious.csv").write_text("keep me\n")

    result = run_cli("--checkpoint", str(ckpt), output=tmp_path / "out")
    assert result.returncode != 0
    assert "refusing to use" in result.stderr
    assert (ckpt / "precious.csv").exists()
