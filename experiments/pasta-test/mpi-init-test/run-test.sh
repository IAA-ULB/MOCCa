#!/bin/bash -l
#SBATCH --job-name pasta-test
#SBATCH --time=2:00:00      
#SBATCH --mem=0        
#SBATCH --ntasks=128 --cpus-per-task=1 --nodes=2
#SBATCH --account=ap_calcua_epicure

. ./ml.sh

srun  ./test