#!/bin/bash
#SBATCH --job-name job-14x14
#SBATCH --output=slurm-%x.%j.out
#SBATCH --time=4:00:00
#SBATCH --mem=0
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=4 --tasks-per-node=64 --cpus-per-task=1

. ../../../vaughan/ml.sh -p -v

mpirun -np 196 python scaling_small.py 4096 8 
