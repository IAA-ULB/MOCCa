#---------------------------------------------------------------------------------
# Testing the equivalence of gradient and direct HFB solvers for a ground state
# calculation of an even-even nucleus.
#
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
#  bash hfb_direct_gradient.sh -p param -s exec1
#
# -p param: specify a parameterization name
# -s exec : specify the  suffix of the serial executable, i.e. the one without MPI
#           Example: "BXL" for "Tantalus.BXL.exe".
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : wouter.ryssens@ulb.be
# Reference commit hash: 7a785f352de2dc4d995dc190f66ced1afc6ff34
#--------------------------------------------------------------------------------
usage() { echo "Usage: $0 -p param -e exec" 1>&2; exit 1; }

while getopts "p:e:" opt; do
    case "${opt}" in
        p)
            param=${OPTARG}
            ;;
        e)
            exec=${OPTARG}
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
Type='HFB'
pairingscheme=$2
/
&evolution
maxiter=400
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
moreconstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=50
moreconstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=10
/
&Cranking
/
EOF
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

#----------------------------------------------------------------------------------
# Setting up the direct HFB calculation
setup_test_env "directHFB" "$exec" "$param"
write_data  $param 0
# Run the calculation
echo "Running $exe with the direct HFB solver"
./$exe < tant.data > $outfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check_direct=$?
# Saving reference values
directE=$(get_total_energy_stdout $outfile)
# ... and tear down this testing environment.
teardown_test_env
#----------------------------------------------------------------------------------
# Gradient calculation
setup_test_env "gradientHFB" "$exec" "$param"
write_data  $param 1
# Run the calculation
echo "Running $exe with the gradient HFB solver"
./$exe < tant.data > $outfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check_gradient=$?
# Saving reference values
gradientE=$(get_total_energy_stdout $outfile)
# ... and tear down this testing environment.
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
echo '------------------------------------------'
echo ' Runtime checks                           '
echo '------------------------------------------'
printf ' direct HFB solver  ->  %1d \n' $tantalus_check_direct
printf ' gradient HFB solver->  %1d \n' $tantalus_check_gradient
tantalus_check=$(( $tantalus_check_direct|| $tantalus_check_gradient ))

# a) Compare total energies with a tolerance of 1 keV
compare_floats $directE              $gradientE 0.001
check_energy=$?
echo '------------------------------------------'
printf ' Energy consistency           ->  %1d \n' $check_energy
echo '------------------------------------------'

exitcode=$(( $tantalus_check || $check_energy ))

printf ' Success?                     ->  %1d \n' $exitcode
echo '------------------------------------------'

# Return exit code 1 if any of the checks failed
exit $exitcode
