#!/usr/bin/env sh

set -e # exit immediately if command gives non-zero exit status
#-------------------------------------------------------------------------------
# Perform a FAM calculation for the isoscalar dipole motion and verify that
# the spurious mode subtraction removes essentially all of the strength.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  Strength real part after subtraction   0.0                       1e-14
#  Strength imag part after subtraction   0.0                       1e-14
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_spurious_mode_subtraction.sh [MF_EXESUFFIX] [FAM_EXESUFFIX] [PARAM] [PAIRING] [-v/--verbose]
#
# where
#   MF_EXESUFFIX   : executable suffix for mean-field calculation
#   FAM_EXESUFFIX   : executable suffix for FAM calculation
#   PARAM          : parameterization name (e.g., t0t3)
#   PAIRING        : 'HF' or 'HFB'
#   -v or --verbose : print the values which are compared.
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
#-------------------------------------------------------------------------------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Default values for options
verbose=false

# Temporary array to hold arguments
args=()

# Parse all arguments for -v or --verbose
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--verbose)
      verbose=true
      shift
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done

# Check we have enough arguments
if [ ${#args[@]} -lt 4 ]; then
  echo "Error: Not enough arguments provided."
  echo "Usage: bash fam_spurious_mode_subtraction.sh [MF_EXESUFFIX] [FAM_EXESUFFIX] [PARAM] [PAIRING] [-v/--verbose]"
  exit 1
fi

# Extract arguments
mf_exesuffix="${args[0]}"
fam_exesuffix="${args[1]}"
param="${args[2]}"
pairing="${args[3]}"

# Validate pairing
if [[ "$pairing" != "HF" && "$pairing" != "HFB" ]]; then
  echo "Error: PAIRING must be either 'HF' or 'HFB'"
  exit 1
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh
LC_NUMERIC="en_US.UTF-8"

# Set up with meaningful naming that includes all parameters
setup_test_env_fam "fam_spurious_mode_subtraction.${pairing}.${param}" "$mf_exesuffix" "$fam_exesuffix" "$param"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1) Run the mean-field calculation

cat << EOF > mf.data
&nucleus
neutrons=8, protons=8
energy_prec=1e-16
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='$param'
/
&pairing
type='$pairing'
/
&evolution
maxiter=300
/
&scfiteration
/
&wfs
nwn = 14, nwp = 14
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='mf.wf'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
EOF

# Run the calculation
./$exe < mf.data > $mfoutfile
# .... and immediately check if MOCCa reported back some error codes
mocca_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2) Run the FAM calculation

cat << EOF > fam.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=8, ny=8, nz=16, dx=0.8
/
&func
name_param='$param'
/
&pairing
type='$pairing'
/
&evolution
/
&scfiteration
/
&wfs
nwn = 14, nwp = 14
/
&IO
InputFilename='mf.wf'
OutputFilename='trash'
allowtransform=.true.
famfile='S_10.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega_min=1.0
omega_max=1.0
smear=1.0
l=1
m=0
maxiter=30
fam_precision=1e-8
!eff_charge_n = +1
!eff_charge_p = -1
/
EOF
# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if MOCCa reported back some error codes
fam_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3) Check that real and imaginary parts of strength after spurious mode subtraction
#     are both vanishingly small (< 1e-14)

# Extract real and imaginary parts of the strength from the S_10.fam file
# Column 4 = S_complex_re, Column 5 = S_complex_im
# Use LC_NUMERIC to ensure proper handling of scientific notation
strength_real=$(LC_NUMERIC="en_US.UTF-8" awk 'NR>0 {if($1 ~ /^[+-]?[0-9]/) {print $4; exit}}' "S_10.fam")
strength_imag=$(LC_NUMERIC="en_US.UTF-8" awk 'NR>0 {if($1 ~ /^[+-]?[0-9]/) {print $5; exit}}' "S_10.fam")

#echo "Strength real part from S_10.fam: $strength_real"
#echo "Strength imag part from S_10.fam: $strength_imag"

# Check if both real and imaginary parts are < 1e-14 in absolute value
# Use awk for comparison since it handles scientific notation natively
ifail=0

# Check real part
if ! LC_NUMERIC="en_US.UTF-8" awk -v val="$strength_real" 'BEGIN {if (val < 0) val = -val; exit (val < 1e-14) ? 0 : 1}'; then
  echo "Result: FAIL - Real part ($strength_real) is >= 1e-14"
  ifail=1
fi

# Check imaginary part
if ! LC_NUMERIC="en_US.UTF-8" awk -v val="$strength_imag" 'BEGIN {if (val < 0) val = -val; exit (val < 1e-14) ? 0 : 1}'; then
  echo "Result: FAIL - Imaginary part ($strength_imag) is >= 1e-14"
  ifail=1
fi

if (( $ifail == 0 )); then
  echo "Result: PASS - Both real ($strength_real) and imaginary ($strength_imag) parts are < 1e-14"
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
fail=$(($mocca_check || $fam_check || $ifail))

if (($fail == 0)) ; then
  echo -e "test FAM spurious mode subtraction :\033[1;32m success \033[0m"
else
  echo -e "test FAM spurious mode subtraction :\033[1;31m failed ! exit status : mf = $mocca_check, fam = $fam_check, strength = $ifail \033[0m"
fi

exit $fail
