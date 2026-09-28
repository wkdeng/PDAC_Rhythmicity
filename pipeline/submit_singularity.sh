#!/bin/bash
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

module load singularity

snakemake --executor cluster-generic --cluster-generic-submit-cmd \
    "${SCRIPT_DIR}/qsub_wrapper.sh" \
    --jobs 50 --snakefile Snakefile_rmats_singularity.smk --rerun-incomplete --configfile config_rmats_singularity.yaml "$@"
