#!/bin/bash
#SBATCH --job-name pasta-test
#SBATCH --time=4:00:00      
#SBATCH --mem=0        
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=1 --tasks-per-node=64 --cpus-per-task=1 

. ../../../vaughan/ml.sh -p -v

#  test how far we can go on a single node

rm -f timings2.txt
touch timings2.txt

mpirun -np 64 python scaling_small.py    2048   8 >> timings.txt
mpirun -np 64 python scaling_small.py    4096  16 >> timings.txt
mpirun -np 64 python scaling_small.py    8192  32 >> timings.txt
mpirun -np 64 python scaling_small.py   16384  64 >> timings.txt
mpirun -np 64 python scaling_small.py   32768 128 >> timings.txt
mpirun -np 64 python scaling_small.py   65536 256 >> timings.txt
mpirun -np 64 python scaling_small.py  131072 512 >> timings.txt
