#--------------------------------------------------------------------------------------
# Perform a calculation starting from cylindrical shape with small number of particles.
# This checks that periodic bc works well and energy converges to the known value.
# Test starts from *.pot file and generates wf in hdf5 format.
# Calculations restart from hdf5 file and varify the converged status.
#
# Note: this calculation is not yet MPI capable, because on astropc19 I have not yet
#       managed to combine MPI and HDF5 support.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - total energy                        -489.575 MeV                 100 keV
#  - proton number                        20                          0.00001
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash spag.sh [EXESUFFIX]
#
# where EXESUFFIX specifies the  suffix of the executable to be used.
# Note: this test relies on the executable having been compiled with HDF5 support.
#            Example: "BXL" for "Tantalus.BXL.exe".
#
# Dependencies: inp_cyl.pot
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : nikolai.shchechilin@ulb.be
# Reference commit hash: e75e6f25c2e2e7223d59043ddf92a6756bb93793
#-------------------------------------------------------------------------------
# These are the hardcoded answers
refE=-489.575 # total energy of pasta in MeV
refZ=20.0     # number of protons

#set -e #this stops the calc if error is somewhere
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# A small function to write equivalent input data
write_data()
{
# Create runtime data
cat << EOF > tant.data
&nucleus
neutrons=80, protons=20
store_derivatives=.false.
/
&mesh
nx=5, ny=5, nz=5, dx=1.25992105
/
&func
name_param='BSkG4'
/
&pairing
Type='BCS'
/
&evolution
strategy='HBSANE'              ! This iterative strategy is set to mimic the
ortho_strategy = 'CHOLESKY'    !  defaults of MPI-enabled calculations.
subspace_rotation = .true.     !
freezeiter=$1
maxiter=$2
printiter=100
D2H_freeze=1000
/
&scfiteration
kerker_k0=0.3
/
&wfs
nwn = 100, nwp = 40
/
&IO
InputFilename=$3
Outputfilename='tant_cyl.hdf5'
/
&MomentParam
/
&Cranking
/
EOF
}
# Basic starting point of all testing scripts
source ../functions.sh
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Set up
setup_test_env_dep "spag" "$1" "BSkG4" "inp_cyl.pot"
write_data  40 400 "'inp_cyl.pot'"
# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
echo "Calculations done"
# .... and immediately check if Tantalus reported back some error codes
tantalus_check1=$?
# Starting the checking
# a) Get the total energy from the STDOUT file
E=$(get_total_energy_stdout $outfile)
echo "fast check energy" $E $refE
# b) Check that the number of protons is correct
protons=$(get_Z_stdout $outfile)
echo "fast check Z" $protons $refZ
# ... and compare with a tolerance of 100 keV to the expected answer
compare_floats $E $refE 0.1
check_energy1=$?
compare_floats $protons $refZ 0.00001
check_Z1=$?
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# restarting with hdf5 file
write_data  0 10 "'tant_cyl.hdf5'"
# Run the calculation
./$exe < tant.data > $outfile.b
# .... and immediately check if Tantalus reported back some error codes
tantalus_check2=$?
# Starting the checking
# a) Get the total energy from the STDOUT file
E=$(get_total_energy_stdout $outfile)
echo "fast check2 energy" $E $refE
# b) Check that the number of protons is correct
protons=$(get_Z_stdout $outfile)
echo "fast check2 Z" $protons $refZ
# ... and compare with a tolerance of 100 keV to the expected answer
compare_floats $E $refE 0.100
check_energy2=$?
compare_floats $protons $refZ 0.00001
check_Z2=$?
# remove working directory and traces of these calculations
teardown_test_env
# printing the checks
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
echo '------------------------------------------------'
echo ' Runtime checks                           '
echo '------------------------------------------------'
printf ' pasta start            ->  %1d \n' $tantalus_check1
printf ' pasta REstart          ->  %1d \n' $tantalus_check2
echo '------------------------------------------------'
echo '------------------------------------------------'
printf ' Total energy consistency            ->  %1d \n' $check_energy1
echo '------------------------------------------------'
printf ' Number of protons consistency       ->  %1d \n' $check_Z1
echo '------------------------------------------------'
echo '------------------------------------------------'
echo ' After restart'
echo '------------------------------------------------'
printf ' Total energy consistency            ->  %1d \n' $check_energy2
echo '------------------------------------------------'
printf ' Number of protons consistency       ->  %1d \n' $check_Z2
echo '------------------------------------------------'
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
exitcode=$(($tantalus_check1 || $check_energy1 || $check_Z1 || $tantalus_check2 || $check_energy2 || $check_Z2 ))
echo '------------------------------------------------'
printf ' Success?                     ->  %1d \n' $exitcode
echo '------------------------------------------'

# Return exit code 1 if any of the checks failed
exit $exitcode
