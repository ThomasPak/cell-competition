#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --job-name=g1-proportion-tests
#SBATCH --output=g1-proportion-tests_%A_%a.out
#SBATCH --error=g1-proportion-tests_%A_%a.err
#SBATCH --array=0-699

#SBATCH --mail-type=all
#SBATCH --mail-user=thomas.pak@linacre.ox.ac.uk

set -euxo pipefail

trap 'date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Exiting $0"' EXIT

date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Entering $0"

python3 generate-parameters-g1-proportion-tests.py $SLURM_ARRAY_TASK_ID | \
            "$CHASTE_BUILD/projects/cell-competition/apps/DeathClockApp"
