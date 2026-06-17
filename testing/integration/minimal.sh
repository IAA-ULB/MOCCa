#-------------------------------------------------------------------------------
# Perform a spherical calculation of O16 with SLy4 in a minimal box to compare
# to known values for the energy and quadrupole deformation.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - total energy                        -128.513125 MeV            1 keV
#  - quadrupole deformation \beta_20        0.0                     0.0001
#  - quadrupole deformation \beta_22        0.0                     0.0001
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash minimal.sh [EXESUFFIX]
#
# where EXESUFFIX specifies the  suffix of the executable to be used
#            Example: "BXL" for "MOCCa.BXL.exe".
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: commit 24470533c34a47af11223514bbecd500b212151a
#-------------------------------------------------------------------------------
# These are the hardcoded answers
refE=-128.513125 # Total energy of O16 in MeV
refB20=0.0       # this nucleus should be REALLY spherical
refB22=0.0

set -e
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env "minimal" "$1" "SLy4"

# Create runtime data
cat << EOF > mocca.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy4'
/
&pairing
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='trash'
/
&MomentParam
/
&Cranking
/
EOF

# Run the calculation
./$exe < mocca.data > $outfile
# .... and immediately check if MOCCa reported back some error codes
mocca_check=$?

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
# a) Get the total energy from the STDOUT file
E=$(get_total_energy_stdout $outfile)
# ... and compare with a tolerance of 1 keV to the expected answer
compare_floats $E $refE 0.001
check_energy=$?
# b) Check that the quadrupole deformation is zero
B20=$(get_B20_stdout $outfile)
compare_floats $B20 $refB20 0.00001
check_B20=$?
B22=$(get_B22_stdout $outfile)
compare_floats $B22 $refB22 0.00001
check_B22=$?
# remove working directory and traces of these calculations
teardown_test_env
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
exit $(($mocca_check || $check_energy || $check_B20 || $check_B22 ))
