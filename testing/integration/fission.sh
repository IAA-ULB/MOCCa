#-------------------------------------------------------------------------------
# Perform a deformed calculation of Pu240 with BSkG3 to compare to known values 
# for the energy, axial quadrupole and collective inertia tensor.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                                Target                    Tolerance
#  ------------------------------------------------------------------------------
#  - total energy  (MeV)                 -1805.647188               0.050
#  - \beta_20                                 0.68421               0.001
#  - I_2020 ( MeV^{-1} b^{-l} [hbar^2])       0.20311               0.001
#  - I_2030 ( MeV^{-1} b^{-l} [hbar^2])       0.00000               0.000
#  - I_3030 ( MeV^{-1} b^{-l} [hbar^2])       0.04728               0.001
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Usage
# ------
#   bash fission.sh [EXESUFFIX]
# 
# where EXESUFFIX identifies the suffix to be used in selecting the executable;
# should be a BSkG-capable one. 
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: commit 1c649912ed5f0e25cea9384e0bb580f0a8eba3f5
#-------------------------------------------------------------------------------
# These are the hardcoded answers
refE=-1805.647188
refB20=0.68421
refI_2020=0.197829
refI_2030=0.000000
refI_3030=0.045706

#set -e
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env "fission" "$1" "BSkG3"

# Create runtime data
cat << EOF > tant.data
&nucleus
neutrons=146, protons=94
energy_prec=1e-06
moment_prec=0.0001
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='BSkG3'
/
&pairing
Type="HFB"
pairingscheme=1
bogofromfile=.false.
/
&evolution
maxiter=500
gradient_safety=1
/
&scfiteration
/
&wfs
nwp=130
nwn=220
/
&IO
InputFilename='init'
OutputFilename='trash'
N_inertia=7
/
&inertia
Inertia_l=2,2,3,3,4,4,4
Inertia_m=0,2,0,2,0,2,4
/
&MomentParam
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=0
Constraint=2180.129895613483
multfromfile=.true.
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=2
Constraint=0.0
multfromfile=.true.
/
&Cranking
/
EOF

# Run the calculation
./$exe < tant.data > $outfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check=$?
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
# a) Get the total energy from the STDOUT file
E=$(get_total_energy_stdout $outfile)
# ... and compare with a tolerance of 50 keV to the expected answer
compare_floats $E $refE 0.050
check_energy=$?
# b) Check that the quadrupole deformation is zero
B20=$(get_B20_stdout $outfile)
compare_floats $B20 $refB20 0.001
check_B20=$?
# c) And finally a few components of the inertia tensor
read I_2020 I_2030 I_3030 <<< $(get_inertia_components $outfile)

# Compare each component with its reference
compare_floats $I_2020 $refI_2020 0.001
check_I_2020=$?

compare_floats $I_2030 $refI_2030 0.001
check_I_2030=$?

compare_floats $I_3030 $refI_3030 0.001
check_I_3030=$?

# remove working directory and traces of these calculations
teardown_test_env
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
overall_check=$(($tantalus_check || $check_energy || $check_B20 || $check_I_2020 || $check_I_2030 || $check_I_3030))

# Print final status
if [ $overall_check -eq 0 ]; then
    echo "All checks passed: SUCCESS"
else
    echo "Some checks failed: FAIL"

    # Specify which ones failed
    [ $tantalus_check -ne 0 ] && echo "  - Tantalus check failed"
    [ $check_energy -ne 0 ] && echo "  - Energy check failed"
    [ $check_B20 -ne 0 ] && echo "  - B20 check failed"
    [ $check_I_2020 -ne 0 ] && echo "  - Inertia component I_2020 check failed"
    [ $check_I_2030 -ne 0 ] && echo "  - Inertia component I_2030 check failed"
    [ $check_I_3030 -ne 0 ] && echo "  - Inertia component I_3030 check failed"
fi

# Exit with the combined status
exit $overall_check
