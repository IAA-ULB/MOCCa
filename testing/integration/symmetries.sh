#---------------------------------------------------------------------------------
# Testing the equivalence of different symmetry modes with multiple calculations
# of Mg24 constrained to a triaxial shape.
#
# More precisely, this testing script runs the exact same calculation with four
# different Tantalus configuration files:
#    - X   : maximally symmetric mode
#    - X-P : broken parity
#    - X-T : broken time-reversal
#    - X-TP: broken time-reversal AND broken parity
#
# The script does not check the results w.r.t. a hard-coded answer but rather
# verifies that the outcomes of all calculations are identical.
#
# The precise configurations can be chosen by the user as input, i.e. "X" in
# the list above can be many different things.
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
#  - improve input to accept longer options?
#     -> this will mean moving away from getopts
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#
# 1.  the "default"
#     Invoking without arguments, i.e.
#        bash symmetries.sh [no arguments]
#     will perform the test for a BXL-style CONFIG with BSkG3 for the flags
#     requested.
#
# 2.  the "manual"
#     Full control through (example)
#        bash symmetries.sh param exec exec_P exec_T exec_TP
#     where
#        - param  : the name of a parameterisation, i.e. BSkG3
#        - exec   : the suffix of the maximally symmetric executable on
#                   your system. Example: "BXL" for "Tantalus.BXL.exe".
#        - exec_P : the suffix of the parity-broken executable on your
#                   system. Example: "BXL-P" for "Tantalus.BXL-P.exe".
#        - exec_T : the suffix of the time-reversal-broken executable
#                   on your system.
#                   Example: "BXL-T" for "Tantalus.BXL-T.exe".
#        - exec_TP: the suffix of the parity+time-reversal-broken executable
#                   on your system.
#                   Example: "BXL-TP" for "Tantalus.BXL-TP.exe".
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : wouter.ryssens@ulb.be
# Complexity           : medium
# Reference commit hash: ec5d48ee9e546d795cc74d6160bf7812cd1d3b93
#--------------------------------------------------------------------------------
logfiletag='symmetries'
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Input checking: this script requires either 5 or 0 positional arguments!
if [[ "$#" -ne 0 && "$#" -ne 5 ]]; then
  echo " Illegal number of parameters!"
  echo " Correct useage is either"
  echo " > bash symmetries.sh"
  echo " OR "
  echo " > bash symmetries.sh param exec exec_P exec_T exec_TP"
  echo " See the comments in the script for more info."
  exit 1
fi

if [[ "$#" -eq 5 ]]; then
# The parameterisation to test is the first argument
param=$1
# The relevant executable suffixes are the next ones
exec=$2     # maximally symmetric configuration
exec_P=$3   # parity-broken configuration
exec_T=$4   # time-reversal broken configuration
exec_TP=$5  # parity AND time-reversal broken configuration
else
param='BSkG3'
exec="BXL"
exec_P="BXL-P"
exec_T="BXL-T"
exec_TP="BXL-TP"
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up the reference calculation: BXL
setup_test_env "$logfiletag" "$exec" "$param"

cat << EOF > tant.data
&nucleus
neutrons=12, protons=12
store_derivatives=$store_derivatives
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type='HFB'
/
&evolution
maxiter=$iter
/
&scfiteration
/
&wfs
nwn = 24, nwp = 24
osc_freq = 0.2, 0.19, 0.18
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
constraint=20
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

# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check=$?
# Saving reference values
refE=$(get_total_energy_stdout $outfile)
# ... and tear down this testing environment.
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Performing a parity-broken calculation
setup_test_env "$logfiletag" "$exec_P" "$param"

cat << EOF > tant.data
&nucleus
neutrons=12, protons=12
store_derivatives=$store_derivatives
/
&mesh
nx=12, ny=12, nz=24, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type='HFB'
/
&evolution
maxiter=$iter
/
&scfiteration
/
&wfs
nwn = 24, nwp = 24
osc_freq = 0.2, 0.19, 0.18
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
constraint=20
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

# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
# ... immediately check if Tantalus reported back some error codes
tantalus_check_P=$?
# Saving reference values
E_parity=$(get_total_energy_stdout $outfile)
# ... remove all trace of these calculations
teardown_test_env
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Performing a time-reversal-broken calculation
setup_test_env "$logfiletag" "$exec_T" "$param"

cat << EOF > tant.data
&nucleus
neutrons=12, protons=12
store_derivatives=$store_derivatives
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type='HFB'
/
&evolution
maxiter=$iter
/
&scfiteration
/
&wfs
nwn = 48, nwp = 48
osc_freq = 0.2, 0.19, 0.18
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
constraint=20
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
# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
# ... immediately check if Tantalus reported back some error codes
tantalus_check_T=$?
# Saving reference values
E_timereversal=$(get_total_energy_stdout $outfile)
# ... remove all trace of these calculations
teardown_test_env
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Performing a time-reversal- and parity-broken calculation
# We need to set-up an initial run
setup_test_env "$logfiletag" "$exec_P" "$param"

cat << EOF > tant.init.data
&nucleus
neutrons=12, protons=12
store_derivatives=$store_derivatives
/
&mesh
nx=12, ny=12, nz=24, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type='HFB'
/
&evolution
maxiter=0
/
&scfiteration
/
&wfs
nwn = 24, nwp = 24
osc_freq = 0.2, 0.19, 0.18
/
&IO
InputFilename='init'
OutputFilename='setting_up.wf'
allowtransform=.true.
/
&MomentParam
moreconstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=20
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
# Run the initial run calculation (we don't keep the output)
./$exe < tant.init.data > /dev/null
# ... but we do keep the wf file!
mv setting_up.wf ../
teardown_test_env

# Change the environment
setup_test_env "$logfiletag" "$exec_TP" "$param"
mv ../setting_up.wf .

# ... and now put the actual data!
cat << EOF > tant.data
&nucleus
neutrons=12, protons=12
store_derivatives=$store_derivatives
/
&mesh
nx=12, ny=12, nz=24, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type='HFB'
/
&evolution
maxiter=$iter
/
&scfiteration
/
&wfs
nwn = 48, nwp = 48
osc_freq = 0.2, 0.19, 0.18
/
&IO
InputFilename='setting_up.wf'
OutputFilename='trash'
allowtransform=.true.
/
&MomentParam
moreconstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=20
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
# Run the calculation
echo "Running $exe"
./$exe < tant.data > $outfile
# ... immediately check if Tantalus reported back some error codes
tantalus_check_TP=$?
# Saving reference values
E_timereversal_parity=$(get_total_energy_stdout $outfile)
# ... remove all trace of these calculations
teardown_test_env
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
echo '---------------------------------'
echo ' Runtime checks                  '
echo '---------------------------------'
printf ' Did %-10s run?       %1d \n' $exec    $tantalus_check
printf ' Did %-10s run?       %1d \n' $exec_P  $tantalus_check_P
printf ' Did %-10s run?       %1d \n' $exec_T  $tantalus_check_T
printf ' Did %-10s run?       %1d \n' $exec_TP $tantalus_check_TP
# a) Compare total energies with a tolerance of 1 keV
# a.1) parity-broken calculation
compare_floats $E_parity              $refE 0.001
check_energy_P=$?
# a.2) time-reversal broken calculation
compare_floats $E_timereversal        $refE 0.001
check_energy_T=$?
# a.3) time-reversal broken calculation
compare_floats $E_timereversal_parity $refE 0.001
check_energy_TP=$?

echo '---------------------------------'
echo ' Energy consistency              '
echo '---------------------------------'
printf " %-10s = %-10s?  %1d \n"  $exec $exec_P  $check_energy_P
printf " %-10s = %-10s?  %1d \n"  $exec $exec_T  $check_energy_P
printf " %-10s = %-10s?  %1d \n"  $exec $exec_TP $check_energy_TP
echo '---------------------------------'

# t_check = Global exit code for correct endings of executables
t_check=$(( $tantalus_check || $tantalus_check_P || $tantalus_check_T || $tantalus_check_TP ))
# e_check = Global exit code for energy comparisons
e_check=$(( $check_energy_P || $check_energy_T || $check_energy_TP ))
# Global exit code: everything needs to pass!
exitcode=$(( $t_check || $e_check ))

echo ' SUCCESS?      ' $exitcode
echo '---------------------------------'

# Return exit code 1 if any of the checks failed
exit $exitcode
