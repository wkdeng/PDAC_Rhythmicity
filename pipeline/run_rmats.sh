#!/bin/bash
###############################################################################
# ChronoTherapy rMATS/DaPars Pipeline - TACC LS6 Launcher
# For pre-aligned BAM cohorts (e.g., TCGA genomic)
#
# Usage:
#   ./run_rmats.sh              # Full pipeline
#   ./run_rmats.sh -n           # Dry run only
#   ./run_rmats.sh --target X   # Run specific target
###############################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Defaults
DRY_RUN=false
TARGET=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        -n|--dry-run)
            DRY_RUN=true
            shift
            ;;
        --target)
            TARGET="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [-n|--dry-run] [--target <snakemake_target>]"
            echo ""
            echo "Options:"
            echo "  -n, --dry-run     Show what would be done without executing"
            echo "  --target TARGET   Run a specific Snakemake target"
            echo "  -h, --help        Show this help message"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Create logs directory
mkdir -p logs

# Redirect Snakemake's temp directory to local /tmp to avoid Lustre metadata
# latency issues (missing tmp files in .snakemake/)
export TMPDIR="/tmp/${USER}_snakemake_$$"
mkdir -p "$TMPDIR"
trap 'rm -rf "$TMPDIR"' EXIT

echo "================================================"
echo "ChronoTherapy rMATS/DaPars Pipeline (TACC LS6)"
echo "================================================"
echo "Working directory: $SCRIPT_DIR"
echo "Snakefile: Snakefile_rmats.smk"
echo "Config: config_rmats.yaml"
echo "Profile: profile/slurm"
echo ""

if $DRY_RUN; then
    echo "=== DRY RUN ==="
    snakemake -s Snakefile_rmats.smk --configfile config_rmats.yaml --profile profile/slurm -n $TARGET
else
    # Dry run first for verification
    echo "Running dry-run check..."
    snakemake -s Snakefile_rmats.smk --configfile config_rmats.yaml --profile profile/slurm -n $TARGET
    echo ""
    echo "Dry run passed. Press Enter to submit jobs, or Ctrl+C to cancel..."
    read

    echo "Submitting jobs to SLURM..."
    snakemake -s Snakefile_rmats.smk --retries 2 --configfile config_rmats.yaml --profile profile/slurm --rerun-incomplete $TARGET

    echo ""
    echo "================================================"
    echo "Pipeline completed!"
    echo "================================================"
    echo "Check job status: squeue -u \$USER"
    echo "View logs in: data/logs/"
fi
