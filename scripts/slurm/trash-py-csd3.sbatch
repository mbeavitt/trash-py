#!/bin/bash
#! Run trash-py on CSD3 (Slurm) so that a genome needing more than one
#! 12-hour allocation finishes across several of them.
#!
#! How it works:
#!   * --checkpoint makes trash-py stream its progress to a directory and
#!     resume from it when the same command runs again.
#!   * --signal=B:USR1@900 asks Slurm to signal this script 15 minutes
#!     before the wall-clock limit. The trap below relays that to trash-py,
#!     which finishes the tasks in flight, commits its state and exits 75.
#!   * exit 75 means "not finished, state saved", so we requeue the job;
#!     any other code is a real success or a real failure.
#!
#! Even with no signal wiring at all, trash-py reads SLURM_JOB_END_TIME and
#! stops ~10 minutes (--time-margin) before the allocation ends. The signal
#! is what makes it react to a limit that changes under it (a requeue onto a
#! shorter partition, say).

#SBATCH --job-name=trash-py
#SBATCH --account=CHANGEME-SL2-CPU
#SBATCH --partition=icelake
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --time=12:00:00
#SBATCH --signal=B:USR1@900
#SBATCH --requeue
#SBATCH --open-mode=append
#SBATCH --output=logs/trash-py-%j.log

set -u

#! CSD3 module environment (adjust to your setup).
. /etc/profile.d/modules.sh
module purge
module load rhel8/default-icl

#! trash-py, clustalo, nhmmer (and mafft for HORs) must all be on PATH.
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda activate trash-py

GENOME=/path/to/genome.fasta
OUTDIR=/path/to/output
CKPT="$OUTDIR/checkpoint"

mkdir -p "$OUTDIR" logs

#! ---------------------------------------------------------------- relay --
#! Slurm signals the batch shell (B:), not the program it started, so the
#! trap has to pass the signal on. Bash only runs a trap between commands,
#! which is why trash-py runs in the background with an explicit `wait`
#! rather than in the foreground.
trash_pid=""
relay_signal() {
    echo "[$(date -Is)] wall-clock limit approaching — asking trash-py to checkpoint"
    [ -n "$trash_pid" ] && kill -USR1 "$trash_pid"
}
trap relay_signal USR1

echo "[$(date -Is)] starting trash-py (job $SLURM_JOB_ID, attempt ${SLURM_RESTART_COUNT:-0})"

trash-py \
    -f "$GENOME" \
    -o "$OUTDIR" \
    -p "${SLURM_CPUS_PER_TASK:-1}" \
    --checkpoint "$CKPT" \
    --time-margin 10m &
trash_pid=$!

#! A trapped signal interrupts `wait` (returning >128), so wait again for
#! the real exit status.
wait "$trash_pid"; status=$?
while [ "$status" -gt 128 ]; do
    wait "$trash_pid"; status=$?
done

#! ------------------------------------------------------------- requeue --
if [ "$status" -eq 75 ]; then
    echo "[$(date -Is)] checkpointed in $CKPT — requeueing job $SLURM_JOB_ID"
    #! Same job ID, same script, and trash-py picks the checkpoint back up.
    #! Needs --requeue above; --open-mode=append keeps the earlier log.
    scontrol requeue "$SLURM_JOB_ID"
    #! Alternatively, submit a fresh job instead of requeueing this one:
    #!     sbatch "$0"
    exit 0
fi

echo "[$(date -Is)] trash-py exited $status"
exit "$status"
