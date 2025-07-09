#!/bin/bash -l
#SBATCH --time=0:5:00      
#SBATCH --mem=0        
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=1 --tasks-per-node=64 --cpus-per-task=1 
#SBATCH --partition=debug
#SBATCH --account=project_465000095

# pass command line arguments to ml-vaughan
. ../../env/lumi/ml.sh
ml

make

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun -n 6 python ../hello_from_PyScalapack.py