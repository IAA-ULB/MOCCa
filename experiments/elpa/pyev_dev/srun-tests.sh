#!/bin/bash
#SBATCH --job-name srun-tests.sh
#SBATCH --output=slurm-%x.%j.out
#SBATCH --time=0:30:00
#SBATCH --mem=0
#SBATCH --account=project_465000095
#SBATCH --nodes=1 --tasks-per-node=4 --cpus-per-task=1
#SBATCH --partition=small

. /pfs/lustrep4/projappl/project_465000095/entijske/tantalus_full/experiments/env/ml.sh

# use pytest -x for stopping at the first failure

export USE_PYELPA=1
echo "USE_PYELPA=${USE_PYELPA}"
srun -n 4 python -m pytest -s tests

# this one fails often because not all Elpa functionality is implemented in pyev (yet)
export USE_PYELPA=0
echo "USE_PYELPA=${USE_PYELPA}"
srun -n 4 python -m pytest -s tests
