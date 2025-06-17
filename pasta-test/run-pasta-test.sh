#!/bin/bash -l
#SBATCH --job-name pasta-test
#SBATCH --time=2:00:00      
#SBATCH --mem=0        
#SBATCH --ntasks=128 --cpus-per-task=1 --nodes=2
#SBATCH --account=ap_calcua_epicure

#--------------------------------------------------------------
# Template script for Tantalus runs on Vaughan
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ml calcua/2024a
ml ScaLAPACK
ml

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


echo "#---------------------------------------------------------------------------------" 
echo "Date:       =  " `date`
echo "Host:       =  " `hostname`
echo "Directory:  =  " `pwd`
echo "JOB_ID:     =  " ${SLURM_JOB_ID}
echo "#---------------------------------------------------------------------------------" 



echo 'Running Tantalus.BXL.exe < data_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3 > test-pasta.out'

# Executing the mean field-code
srun ./Tantalus.BXL.exe < data_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3 

echo 'DONE!'
