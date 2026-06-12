#------------------------------------------------------------------------------------
# Testing the equivalence of results obtained for O16 with periodic and antiperiodic
# boundary conditions, respectively. Note that this test does not check the results
# w.r.t. to a hardcoded reference value, but rather checks that results are identical
# in both modes.
#
# Notes
# 1. Because of the different implementation for terms with derivatives in the
#    functional, it is normal to have small differences between calculations
#    with periodic and antiperiodic boundary conditions for dx ~ 0.8 fm - 1.0 fm.
#    To get agreement at the level of a keV, you need to go to lower values of
#    mesh spacing; this test uses dx = 0.65 fm.
#
# 2. There are several parts of the total energy in a typical .param file that
#    would invalidate this test. Among them are
#    - corrections for spurious motion: rotational correction, vibrational correction,
#                                       center of mass corrections, ...
#    - Coulomb: the pasta calculations assume a homogeneous background of electrons
#               that the nuclear calculations do not.
#
#    Hence use extreme caution when selection a .param file to run this test for;
#    the BSkG3_naked.param is an example for which the test passes.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#
#  Quantity                              Answer                     Tolerance
#  --------                              ------                     ---------
#  Total energy                       < the first result obtained >   1 keV
#
#  TODO:
#  - add other observables
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#
#  bash boundary_conditions.sh -p param -e exec [OTHER FLAGS]
#
# -p param: specify a parameterization name
# -e exec : specify the  suffix of the executable with antiperiodic boundary conditions
#            Example: "BXL.NUCLEI" for "MOCCa.BXL.NUCLEI.exe".
# -f exec  : specify the  suffix of the executable with periodic boundary conditions
#            Example: "BXL.PASTA" for "MOCCa.BXL.PASTA.exe".
#
#  Attention: the exe being called should be able to auto-initialise, i.e. to
#             start from scratch without reading a .wf file!
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : wouter.ryssens@ulb.be
# Reference commit hash: bfaa9cbbad6f1aa28e2d70613ee05a4da2a04af5
#--------------------------------------------------------------------------------
usage() { echo "Usage: $0 -p param -e exec" 1>&2; exit 1; }

while getopts "p:e:f:" opt; do
    case "${opt}" in
        p)
            param=${OPTARG}
            ;;
        e)
            exec_anti=${OPTARG}
            ;;
        f)
            exec_period=${OPTARG}
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
neutrons=8, protons=8
/
&mesh
nx=15, ny=15, nz=15, dx=0.65
/
&func
name_param="$1"
/
&pairing
Type='HF'
/
&evolution
maxiter=200
/
&scfiteration
/
&wfs
ini_strategy='NILSSON'
nwn = 60, nwp = 60
osc_freq = 0.2, 0.2, 0.2
! A sufficiently large nwn and nwp to ensure a broken time-reversal executable will not run into trouble.
/
&IO
InputFilename='init'
OutputFilename='trash'
allowtransform=.true.
/
&MomentParam
/
&MomentConstraint
/
&MomentConstraint
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
setup_test_env "boundary_conditions" "$exec_anti" "$param"
write_data  $param
# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
# .... and immediately check if MOCCa reported back some error codes
mocca_check=$?
# Saving reference values
refE=$(get_total_energy_stdout $outfile)
echo $refE $outfile
# ... and tear down this testing environment.
teardown_test_env
#----------------------------------------------------------------------------------
# Setting up the reference calculation
setup_test_env "boundary_conditions" "$exec_period" "$param"
write_data  $param
# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
# .... and immediately check if MOCCa reported back some error codes
mocca_check_false=$?
# Saving reference values
testE=$(get_total_energy_stdout $outfile)
# ... and tear down this testing environment.
teardown_test_env

echo $testE $outfile

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
echo '------------------------------------------------'
echo ' Runtime checks                           '
echo '------------------------------------------------'
printf ' antiperiodic boundary conditions   ->  %1d \n' $mocca_check
printf ' periodic     boundary conditions   ->  %1d \n' $mocca_check_false
# a) Compare total energies with a tolerance of 1 keV
compare_floats $testE              $refE 0.001
check_energy=$?

echo '------------------------------------------------'
printf ' Energy consistency                 ->  %1d \n' $check_energy
echo '------------------------------------------------'

exitcode=$(( $mocca_check || mocca_check_false || $check_energy ))
echo '------------------------------------------------'
printf ' Success?                     ->  %1d \n' $exitcode
echo '------------------------------------------'

# Return exit code 1 if any of the checks failed
exit $exitcode
