#!/bin/bash
# Wrapper around qsub that excludes down nodes by picking a healthy one
# Usage: called by snakemake as the cluster submit command

EXCLUDE="node03 node06 node11 node12 node13 node14 node15 node16 node17 node18 node19 node20 \
node21 node22 node23 node24 node25 node26 node27 node28 node29 node30"

# Get a random free node not in the exclude list
GOOD_NODE=$(pbsnodes -a 2>/dev/null | awk '
    /^[a-z]/ { node=$1 }
    /state = free/ || /state = job-exclusive/ { print node }
' | grep -vFf <(printf '%s\n' $EXCLUDE) | shuf | head -1)

if [ -z "$GOOD_NODE" ]; then
    echo "No available nodes found" >&2
    exit 1
fi

mkdir -p logs

qsub -V -j oe -o logs/ -e logs/ \
    -l select=1:ncpus=20:mem=200gb:vnode=${GOOD_NODE} \
    "$@"
