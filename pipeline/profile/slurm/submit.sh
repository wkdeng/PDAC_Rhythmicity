#!/bin/bash
# Submit wrapper for Snakemake cluster-generic executor
# Filters sbatch output to return only the job ID

output=$(sbatch "$@" 2>/dev/null)
# Extract just the numeric job ID from sbatch --parsable output
echo "$output" | grep -oE '[0-9]+' | tail -1
