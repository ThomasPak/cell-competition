#!/bin/bash
# Set error checking
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

# Generate parameters and run simulation
python3 generate-vertex-parameters-heterotypic-proliferation-regimes-cross-section-II.py $sim_num | \
            "$CHASTE_BUILD/projects/cell-competition/apps/DeathClockApp"
