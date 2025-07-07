#!/bin/bash
#SBATCH --job-name job-32x32-23170
#SBATCH --output=slurm-%x.%j.out
#SBATCH --time=10:00:00
#SBATCH --mem=0
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=16 --tasks-per-node=64 --cpus-per-task=1

. ../../../vaughan/ml.sh -p -v

mpirun -np 1024 python pdsyev.py 23170 8 
