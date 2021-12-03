#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --job-name=exponential-survival-probability
#SBATCH --output=exponential-survival-probability_%A_%a.out
#SBATCH --error=exponential-survival-probability_%A_%a.err
#SBATCH --array=0-299

#SBATCH --mail-type=all
#SBATCH --mail-user=thomas.pak@linacre.ox.ac.uk

set -euxo pipefail

trap 'date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Exiting $0"' EXIT

date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Entering $0"

python3 generate-vertex-parameters-exponential-survival-probability.py $SLURM_ARRAY_TASK_ID | \
            "$CHASTE_BUILD/projects/cell-competition/apps/DeathClockApp"
