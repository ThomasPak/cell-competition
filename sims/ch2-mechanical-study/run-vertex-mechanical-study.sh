#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=1:00:00
#SBATCH --partition=short
#SBATCH --job-name=mechanical-study
#SBATCH --output=mechanical-study_%A_%a.out
#SBATCH --error=mechanical-study_%A_%a.err

#SBATCH --mail-type=all
#SBATCH --mail-user=thomas.pak@linacre.ox.ac.uk

# Set error checking
set -euxo pipefail

# Write date when exiting
trap 'date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Exiting $0"' EXIT

# Write date on onset
date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Entering $0"

# Process arguments
if [ $# -ne 2 ]
then
    echo "Usage: $0 SIM_START SIM_END"
    exit 2
fi

sim_start="$1"
sim_end="$2"

# Assert that size of simulation range is equal to slurm array task count
if [ "$((sim_end - sim_start))" -ne "$SLURM_ARRAY_TASK_COUNT" ]
then
    echo "Error: size of simulation range [$sim_start, $sim_end) is not equal to SLURM_ARRAY_TASK_COUNT: $SLURM_ARRAY_TASK_COUNT"
    exit 2
fi

# Compute simulation number
sim_num="$((sim_start + (SLURM_ARRAY_TASK_ID - SLURM_ARRAY_TASK_MIN)))"
line_num="$((sim_num + 1))"

# Generate parameters and run simulation
python3 generate-vertex-parameters-mechanical-study.py | \
    sed -n "$line_num"p | tr ',' '\n' | \
    "$CHASTE_BUILD/projects/cell-competition/apps/CellCompetitionApp" \
    random-movement-parameter=0 \
    output-directory=mechanical-tests \
    simulation-time=250 \
    dt=0.005 \
    sampling-timestep-multiple=200 \
    num-cells-across=6 \
    num-cells-up=6 \
    cell-a-g1-duration=30.0 \
    cell-a-g2-duration=70.0
