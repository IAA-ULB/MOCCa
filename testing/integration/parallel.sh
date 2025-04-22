#---------------------------------------------------------------------------------
# Testing the equivalence of MPI and serial versions of the code with a simple
# (constrained) Cr48 calculation.
#
# Note:
#      - this test relies implicitly on the mpirun command being compatible
#        with the MPI library you linked to in the compilation of the executable.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#
#  Quantity                              Answer                     Tolerance
#  --------                              ------                     ---------
#  Total energy                       < the first result obtained >   1 keV
#
#  TODO:
#  - add other observables
#  - enable this test for configurations with different symmetry options
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#
#  bash parallel.sh -p param -s exec1 -m exec2 -r np
#
# -p param: specify a parameterization name
# -s exec : specify the  suffix of the serial executable, i.e. the one without MPI
#           Example: "BXL" for "Tantalus.BXL.exe".
# -m exec : specify the  suffix of the parallel executable, i.e. the one with MPI
# -r  np  : the number of MPI ranks to use
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : wouter.ryssens@ulb.be
# Complexity           : medium
# Reference commit hash: 486028f3177530a238062fee11d5392677e83578
#--------------------------------------------------------------------------------
usage() { echo "Usage: $0 -p param -s exec -m exec -r np" 1>&2; exit 1; }

while getopts "p:s:m:r:" opt; do
    case "${opt}" in
        p)
            param=${OPTARG}
            ;;
        s)
            exec_serial=${OPTARG}
            ;;
        m)
            exec_mpi=${OPTARG}
            ;;
        r)
            MPI_RANKS=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# A small function to write equivalent input data
write_data()
{
cat << EOF > tant.data
&nucleus
neutrons=24, protons=24
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param="$1"
/
&pairing
Type='BCS'
/
&evolution
maxiter=200
freezeiter=0
/
&scfiteration
/
&wfs
nwn = 40, nwp = 40
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='trash'
/
&MomentParam
!moreconstraints=.true.
/
&MomentConstraint
!l=2
!m=0
!constraint=50
!moreconstraints=.true.
/
&MomentConstraint
!l=2
!m=2
!constraint=10
/
&Cranking
/
EOF
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

#----------------------------------------------------------------------------------
# Setting up the reference calculation
setup_test_env "serial" "$exec_serial" "$param"
write_data  $param
# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check_serial=$?
# Saving reference values
serialE=$(get_total_energy_stdout $outfile)
# ... and tear down this testing environment.
teardown_test_env
#----------------------------------------------------------------------------------
# Setting up the reference calculation
setup_test_env "parallel" "$exec_mpi" "$param"
write_data  $param
# Run the calculation
echo "Running $exe"
mpirun -n $MPI_RANKS ./$exe < tant.data > $outfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check_parallel=$?
# Saving reference values
parallelE=$(get_total_energy_stdout $outfile)
# ... and tear down this testing environment.
teardown_test_env


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
echo '------------------------------------------'
echo ' Runtime checks                           '
echo '------------------------------------------'
printf ' serial    ->  %1d \n' $tantalus_check_serial
printf ' parallel  ->  %1d \n' $tantalus_check_parallel
# a) Compare total energies with a tolerance of 1 keV
compare_floats $serialE              $parallelE 0.001
check_energy=$?

echo '------------------------------------------'
printf ' Energy consistency           ->  %1d \n' $check_energy
echo '------------------------------------------'

exitcode=$(( $tantalus_check_serial || tantalus_check_parallel || $check_energy ))

printf ' Success?                     ->  %1d \n' $exitcode
echo '------------------------------------------'

# Return exit code 1 if any of the checks failed
exit $exitcode
