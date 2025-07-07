#!/bin/bash
#SBATCH --job-name job-16x16-131067
#SBATCH --output=slurm-%x.%j.out
#SBATCH --time=10:00:00
#SBATCH --mem=0
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=4 --tasks-per-node=64 --cpus-per-task=1

. ../../../vaughan/ml.sh -p -v

mpirun -np 256 python pdsyev.py 131067 8 
