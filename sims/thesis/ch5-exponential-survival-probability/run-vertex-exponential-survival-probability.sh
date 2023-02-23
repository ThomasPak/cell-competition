#!/bin/bash
set -euxo pipefail

trap 'date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Exiting $0"' EXIT

date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Entering $0"

# Process arguments
if [ $# -ne 1 ]
then
    echo "Usage: $0 SIM_NUM"
    exit 2
fi

sim_num="$1"

python3 generate-vertex-parameters-exponential-survival-probability.py $sim_num | \
            "$CHASTE_BUILD/projects/cell-competition/apps/DeathClockApp"
