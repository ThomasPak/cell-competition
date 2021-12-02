#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --job-name=proliferation-regime-1-tests
#SBATCH --output=proliferation-regime-1-tests_%A_%a.out
#SBATCH --error=proliferation-regime-1-tests_%A_%a.err
#SBATCH --array=0-899

#SBATCH --mail-type=all
#SBATCH --mail-user=thomas.pak@linacre.ox.ac.uk

set -euxo pipefail

trap 'date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Exiting $0"' EXIT

date "+%Y-%d-%m %T" | tr -d "\n" ; echo "::Entering $0"

python3 generate-parameters-proliferation-regime-1-tests.py $SLURM_ARRAY_TASK_ID | \
            "$CHASTE_BUILD/projects/cell-competition/apps/DeathClockApp"
