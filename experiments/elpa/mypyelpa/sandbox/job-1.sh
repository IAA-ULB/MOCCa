#!/bin/bash
#SBATCH --job-name job-1
#SBATCH --time=0:10:00      
#SBATCH --mem=0        
#SBATCH --reservation=rocky9
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=1 --tasks-per-node=1 --cpus-per-task=1 

. ../../env/ml.sh -p -v

export LD_LIBRARY_PATH=./${LD_LIBRARY_PATH}
# make lib
python t.py