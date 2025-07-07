#!/bin/bash
#SBATCH --job-name job-1
#SBATCH --time=4:00:00      
#SBATCH --mem=0        
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=1 --tasks-per-node=64 --cpus-per-task=1 

. ../../../vaughan/ml.sh -p -v

#  test how far we can go on a single node

timing_file='timings-1'

rm -f ${timing_file}.txt
touch ${timing_file}.txt

mpirun -np 49 python scaling_small.py    2048  8 
# mpirun -np 64 python scaling_small.py    4096  8 
# mpirun -np 64 python scaling_small.py   8192  8 >> ${timing_file}.txt
