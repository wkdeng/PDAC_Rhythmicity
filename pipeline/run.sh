#!/bin/bash
###############################################################################
# ChronoTherapy RNA-seq Pipeline - TACC LS6 Launcher
#
# Usage:
#   ./run.sh              # Full pipeline
#   ./run.sh -n           # Dry run only
#   ./run.sh --target X   # Run specific target
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

echo "=========================================="
echo "ChronoTherapy RNA-seq Pipeline (TACC LS6)"
echo "=========================================="
echo "Working directory: $SCRIPT_DIR"
echo "Profile: profile/slurm"
echo ""

if $DRY_RUN; then
    echo "=== DRY RUN ==="
    snakemake --profile profile/slurm -n $TARGET
else
    # Dry run first for verification
    echo "Running dry-run check..."
    snakemake --profile profile/slurm -n $TARGET
    echo ""
    echo "Dry run passed. Press Enter to submit jobs, or Ctrl+C to cancel..."
    read

    echo "Submitting jobs to SLURM..."
    snakemake --profile profile/slurm --rerun-incomplete $TARGET

    echo ""
    echo "=========================================="
    echo "Pipeline completed!"
    echo "=========================================="
    echo "Check job status: squeue -u \$USER"
    echo "View logs in: logs/"
fi
