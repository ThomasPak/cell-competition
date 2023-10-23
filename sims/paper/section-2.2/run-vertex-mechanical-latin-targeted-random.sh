#!/bin/bash
set -euxo pipefail

# Write date when exiting
trap 'date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Exiting $0"' EXIT

# Write date on onset
date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Entering $0"

# Process arguments
if [ $# -ne 1 ]
then
    echo "Usage: $0 SIM_NUM"
    exit 2
fi

sim_num="$1"

# Compute line number number
line_num="$((sim_num + 1))"

# Generate parameters and run simulation
cat pregenerated-parameters-mechanical-latin-targeted-random.txt | \
    sed -n "$line_num"p | tr ',' '\n' | \
    "$CHASTE_BUILD/projects/cell-competition/apps/CellCompetitionApp" \
    random-movement-parameter=0 \
    output-directory=mechanical-latin-targeted-random \
    simulation-time=250 \
    dt=0.005 \
    sampling-timestep-multiple=200 \
    num-cells-across=6 \
    num-cells-up=6
