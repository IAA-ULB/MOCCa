#!/bin/bash -l
#SBATCH --job-name pasta-test
#SBATCH --time=4:00:00      
#SBATCH --mem=0        
#SBATCH --account=ap_calcua_epicure
#SBATCH --nodes=16 --ntasks-per-node=64 --cpus-per-task=1 

# pass command line arguments to ml-vaughan
. ./ml-vaughan.sh $@

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "#---------------------------------------------------------------------------------" 
echo "Date:       =  " `date`
echo "Host:       =  " `hostname`
echo "Directory:  =  " `pwd`
echo "JOB_ID:     =  " ${SLURM_JOB_ID}
echo "#---------------------------------------------------------------------------------" 

EXE="Tantalus.BXL.exe"
INP="data_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3"
echo 'Running ${EXE} < ${INP}'

srun ${EXE} < ${INP}

echo 'DONE!'
