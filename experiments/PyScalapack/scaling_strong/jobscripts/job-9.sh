#!/bin/bash
#SBATCH --job-name job-4
#SBATCH --time=4:00:00      
#SBATCH --mem=0        
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=9 --tasks-per-node=64 --cpus-per-task=1 

. ../../../vaughan/ml.sh -p -v

#  test how far we can go on a single node

# mpirun -np 256 python scaling_small.py    2048  8 
# mpirun -np 256 python scaling_small.py    4096  8 
# mpirun -np 256 python scaling_small.py    8192  8 
# mpirun -np 256 python scaling_small.py   16384  8 
mpirun -np 256 python scaling_small.py   32768  8 
