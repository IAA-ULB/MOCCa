#---------------------------------------------------------------------------------
# Testing the equivalence of results obtained with and without storing the
# derivatives of single-particle wavefunctions for 48Cr constrained to a
# triaxial deformation with HFB.
#
# Note that this test does not check the results w.r.t. to a hardcoded reference
# value, but rather checks that results are identical in both modes.
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
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#
#  bash store_derivatives.sh -p param -e exec [OTHER FLAGS]
#
# -p param: specify a parameterization name
# -e exec : specify the  suffix of the executable
#           Example: "BXL" for "MOCCa.BXL.exe".
#
#  Attention: the exe being called should be able to auto-initialise, i.e. to
#             start from scratch without reading a .wf file!
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: 22b6c4b464086ec74f57b72d95a15e11a37a06ec
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
store_derivatives=$1
/
&mesh
nx=12, ny=12, nz=20, dx=1.0
! A large enough number of points in the z-direction to ensure a broken parity executable will not run into trouble.
/
&func
name_param="$2"
/
&pairing
Type='HFB'
/
&evolution
maxiter=200
/
&scfiteration
/
&wfs
nwn = 60, nwp = 60
osc_freq = 0.2, 0.19, 0.18
! A sufficiently large nwn and nwp to ensure a broken time-reversal executable will not run into trouble.
/
&IO
InputFilename='init'
OutputFilename='trash'
allowtransform=.true.
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
# Setting up the reference calculation
setup_test_env "store_derivatives=true" "$exec" "$param"
write_data '.true.' $param
# Run the calculation
echo "Running $exe with store_derivatives=true"
./$exe < tant.data > $outfile
# .... and immediately check if MOCCa reported back some error codes
mocca_check=$?
# Saving reference values
refE=$(get_total_energy_stdout $outfile)
# ... and tear down this testing environment.
teardown_test_env
#----------------------------------------------------------------------------------
# Setting up the reference calculation
setup_test_env "store_derivatives=false" "$exec" "$param"
write_data '.false.' $param
# Run the calculation
echo "Running $exe with store_derivatives=false"
./$exe < tant.data > $outfile
# .... and immediately check if MOCCa reported back some error codes
mocca_check_false=$?
# Saving reference values
testE=$(get_total_energy_stdout $outfile)
# ... and tear down this testing environment.
teardown_test_env


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
echo '------------------------------------------'
echo ' Runtime checks                           '
echo '------------------------------------------'
printf ' store_derivatives = .true.   ->  %1d \n' $mocca_check
printf ' store_derivatives = .false.  ->  %1d \n' $mocca_check_false
# a) Compare total energies with a tolerance of 1 keV
compare_floats $testE              $refE 0.001
check_energy=$?

echo '------------------------------------------'
printf ' Energy consistency           ->  %1d \n' $check_energy
echo '------------------------------------------'

exitcode=$(( $mocca_check || mocca_check_false || $check_energy ))

printf ' Success?                     ->  %1d \n' $exitcode
echo '------------------------------------------'

# Return exit code 1 if any of the checks failed
exit $exitcode
