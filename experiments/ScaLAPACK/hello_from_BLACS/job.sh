#!/bin/bash -l
#SBATCH --time=4:00:00      
#SBATCH --mem=0        
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=1 --tasks-per-node=64 --cpus-per-task=1 

# pass command line arguments to ml-vaughan
. ../../pasta-test/ml-vaughan.sh $@
ml

make

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun exe