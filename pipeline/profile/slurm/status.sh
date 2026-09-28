#!/bin/bash
# Query SLURM job status for Snakemake cluster-generic executor
# Returns: running, success, or failed

jobid="$1"

if [ -z "$jobid" ]; then
    echo "failed"
    exit 0
fi

# First check squeue (works for PENDING and RUNNING jobs)
squeue_state=$(squeue -j "$jobid" -h -o "%T" 2>/dev/null)

if [ -n "$squeue_state" ]; then
    # Job is still in the queue (PENDING, RUNNING, etc.)
    echo "running"
    exit 0
fi

# Job not in squeue, check sacct for final state
# Sleep briefly to allow sacct to catch up
sleep 2
state=$(sacct -j "$jobid" --format=State --noheader | head -1 | awk '{print $1}')

case "$state" in
    COMPLETED)
        echo "success"
        ;;
    RUNNING|PENDING|CONFIGURING|COMPLETING|REQUEUED|SUSPENDED)
        echo "running"
        ;;
    "")
        # No info yet, assume still running
        echo "running"
        ;;
    *)
        echo "failed"
        ;;
esac
