#!/bin/bash
#SBATCH --job-name pasta-test
#SBATCH --time=4:00:00      
#SBATCH --mem=0        
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=1 --tasks-per-node=64 --cpus-per-task=1 

. ../../../vaughan/ml.sh -p -v

rm -f timings.txt
touch timings.txt

mpirun -np  4 python scaling_small.py    64  8 >> timings.txt
mpirun -np 16 python scaling_small.py    64  8 >> timings.txt
mpirun -np 64 python scaling_small.py    64  8 >> timings.txt


mpirun -np  4 python scaling_small.py   128  8 >> timings.txt
mpirun -np 16 python scaling_small.py   128  8 >> timings.txt
mpirun -np 64 python scaling_small.py   128  8 >> timings.txt

mpirun -np  4 python scaling_small.py   128 16 >> timings.txt
mpirun -np 16 python scaling_small.py   128 16 >> timings.txt
mpirun -np 64 python scaling_small.py   128 16 >> timings.txt


mpirun -np  4 python scaling_small.py   256  8 >> timings.txt
mpirun -np 16 python scaling_small.py   256  8 >> timings.txt
mpirun -np 64 python scaling_small.py   256  8 >> timings.txt

mpirun -np  4 python scaling_small.py   256 16 >> timings.txt
mpirun -np 16 python scaling_small.py   256 16 >> timings.txt
mpirun -np 64 python scaling_small.py   256 16 >> timings.txt

mpirun -np  4 python scaling_small.py   256 32 >> timings.txt
mpirun -np 16 python scaling_small.py   256 32 >> timings.txt
mpirun -np 64 python scaling_small.py   256 32 >> timings.txt


mpirun -np  4 python scaling_small.py   512  8 >> timings.txt
mpirun -np 16 python scaling_small.py   512  8 >> timings.txt
mpirun -np 64 python scaling_small.py   512  8 >> timings.txt

mpirun -np  4 python scaling_small.py   512 16 >> timings.txt
mpirun -np 16 python scaling_small.py   512 16 >> timings.txt
mpirun -np 64 python scaling_small.py   512 16 >> timings.txt

mpirun -np  4 python scaling_small.py   512 32 >> timings.txt
mpirun -np 16 python scaling_small.py   512 32 >> timings.txt
mpirun -np 64 python scaling_small.py   512 32 >> timings.txt


mpirun -np  4 python scaling_small.py   1024  8 >> timings.txt
mpirun -np 16 python scaling_small.py   1024  8 >> timings.txt
mpirun -np 64 python scaling_small.py   1024  8 >> timings.txt

mpirun -np  4 python scaling_small.py   1024 16 >> timings.txt
mpirun -np 16 python scaling_small.py   1024 16 >> timings.txt
mpirun -np 64 python scaling_small.py   1024 16 >> timings.txt

mpirun -np  4 python scaling_small.py   1024 32 >> timings.txt
mpirun -np 16 python scaling_small.py   1024 32 >> timings.txt
mpirun -np 64 python scaling_small.py   1024 32 >> timings.txt

mpirun -np  4 python scaling_small.py   1024 32 >> timings.txt
mpirun -np 16 python scaling_small.py   1024 32 >> timings.txt
mpirun -np 64 python scaling_small.py   1024 32 >> timings.txt


mpirun -np  4 python scaling_small.py   2048  8 >> timings.txt
mpirun -np 16 python scaling_small.py   2048  8 >> timings.txt
mpirun -np 64 python scaling_small.py   2048  8 >> timings.txt

mpirun -np  4 python scaling_small.py   2048 16 >> timings.txt
mpirun -np 16 python scaling_small.py   2048 16 >> timings.txt
mpirun -np 64 python scaling_small.py   2048 16 >> timings.txt

mpirun -np  4 python scaling_small.py   2048 32 >> timings.txt
mpirun -np 16 python scaling_small.py   2048 32 >> timings.txt
mpirun -np 64 python scaling_small.py   2048 32 >> timings.txt

mpirun -np  4 python scaling_small.py   4096  8 >> timings.txt
mpirun -np 16 python scaling_small.py   4096  8 >> timings.txt
mpirun -np 64 python scaling_small.py   4096  8 >> timings.txt

mpirun -np  4 python scaling_small.py   4096 16 >> timings.txt
mpirun -np 16 python scaling_small.py   4096 16 >> timings.txt
mpirun -np 64 python scaling_small.py   4096 16 >> timings.txt

mpirun -np  4 python scaling_small.py   4096 32 >> timings.txt
mpirun -np 16 python scaling_small.py   4096 32 >> timings.txt
mpirun -np 64 python scaling_small.py   4096 32 >> timings.txt

mpirun -np  4 python scaling_small.py   4096 32 >> timings.txt
mpirun -np 16 python scaling_small.py   4096 32 >> timings.txt
mpirun -np 64 python scaling_small.py   4096 32 >> timings.txt

mpirun -np  4 python scaling_small.py   4096 64 >> timings.txt
mpirun -np 16 python scaling_small.py   4096 64 >> timings.txt
mpirun -np 64 python scaling_small.py   4096 64 >> timings.txt


mpirun -np  4 python scaling_small.py   8192  8 >> timings.txt
mpirun -np 16 python scaling_small.py   8192  8 >> timings.txt
mpirun -np 64 python scaling_small.py   8192  8 >> timings.txt

mpirun -np  4 python scaling_small.py   8192 16 >> timings.txt
mpirun -np 16 python scaling_small.py   8192 16 >> timings.txt
mpirun -np 64 python scaling_small.py   8192 16 >> timings.txt

mpirun -np  4 python scaling_small.py   8192 32 >> timings.txt
mpirun -np 16 python scaling_small.py   8192 32 >> timings.txt
mpirun -np 64 python scaling_small.py   8192 32 >> timings.txt

mpirun -np  4 python scaling_small.py   8192 32 >> timings.txt
mpirun -np 16 python scaling_small.py   8192 32 >> timings.txt
mpirun -np 64 python scaling_small.py   8192 32 >> timings.txt

mpirun -np  4 python scaling_small.py   8192 64 >> timings.txt
mpirun -np 16 python scaling_small.py   8192 64 >> timings.txt
mpirun -np 64 python scaling_small.py   8192 64 >> timings.txt

