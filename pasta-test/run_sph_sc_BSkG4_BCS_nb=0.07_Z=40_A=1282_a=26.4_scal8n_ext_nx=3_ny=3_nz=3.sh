#!/bin/bash -l
#SBATCH --job-name sph_sc_0.07_333
#SBATCH --time=2:00:00      
#SBATCH --mem=0          
#SBATCH --ntasks=1024 --nodes=8  
#SBATCH --cpus-per-task=1        
#SBATCH --account=project_465001242 
#SBATCH --partition=standard

#--------------------------------------------------------------
# Template script for Tantalus runs on LUMI
# - cray compiler options everywhere
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Preliminary stuff
# a) Get the right Programming environment loaded, just in case
module load PrgEnv-cray
# b) Make sure Python can be used to run auxiliary scripts
module load cray-python
# c) Manage large allocations
module load craype-hugepages2M
# d) HDF5 library
module load cray-hdf5-parallel
# e) Make this explicit, just in case
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
#------------------------------------------------------------------------------------------
#Creating workspace
SCRATCHDIR=/scratch/project_465001242/
#------------------------------------------------------------------------------------------

#Checking if we are dealing with a SLURM job or not
if [ -z ${SLURM_JOB_ID} ]
then
  SLURM_JOB_ID='NOTSUBMITTED'
fi
# not an array job
TEMPDIR=$SCRATCHDIR/${SLURM_JOB_ID}


#making directory if does not exist
if [[ ! -d $TEMPDIR ]]
then
    mkdir $TEMPDIR
fi
#-----------------------------------------------------------------------------------------
#Going to a proper starting place
#usually the nucleus under consideration
cd /users/shchechi/pasta_mpi

echo "#---------------------------------------------------------------------------------" 
echo "Date:       =  " `date`
echo "Host:       =  " `hostname`
echo "Directory:  =  " `pwd`
echo "JOB_ID:     =  " ${SLURM_JOB_ID}
echo "#---------------------------------------------------------------------------------" 

#-----------------------------------------------------------------------------------------
#Copying needed files to temporary directory
cp $HOME/tantalus_full/exec/Tantalus.BXL.exe  $TEMPDIR/
cp data_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3  $TEMPDIR/
cp inp_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3.pot  $TEMPDIR/
cp BSkG4.param  $TEMPDIR/

#----------------------------------------------------------------------------------------
# Moving to the working directory
cd $TEMPDIR

echo 'Running Tantalus.BXL.exe < data_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3 > out_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3'
# Executing the mean field-code
srun ./Tantalus.BXL.exe < data_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3 > out_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3
echo 'DONE!'

#         /scratch/project_465001242/res_den_pot/fold_den_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3
#         /scratch/project_465001242/res_den_pot/fold_pot_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3
# -v iter*den   /scratch/project_465001242/res_den_pot/fold_den_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3/
# -v iter*pot   /scratch/project_465001242/res_den_pot/fold_pot_sph_sc_BSkG4_BCS_nb=0.07_Z=40_A=1282_a=26.4_scal8n_ext_nx=3_ny=3_nz=3/
cp -v out*    /users/shchechi/pasta_mpi/res_out/
cp -v den*    /scratch/project_465001242/res_den_pot/
cp -v pot*    /scratch/project_465001242/res_den_pot/
cp -v tant*   /scratch/project_465001242/res_wf/ 

cd ..
rm -r $TEMPDIR
echo '--------------------------------------------'
###########################################################################################