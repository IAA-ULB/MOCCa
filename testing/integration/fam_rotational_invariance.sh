#!/usr/bin/env sh
#--------------------------------------------------------------------------------
# Verify that FAMQRPA calculations are rotationally invariant for a spherically
# symmetric configuration of Ti42 in a small modelspace by comparing Q20 and
# Q2+2 components.
#--------------------------------------------------------------------------------
# This script tests:
#
#   |(complex) Strength Q20 - Strength Q22| < | 0.15 |
#
# Constraints are imposed on the MF calculation to get a solution closer to
# spherical symmetry.
#
# TODO: check if this precision can be improved; perhaps the modelspace
#       is too small/not perfectly spherical? The latter seems the case,
#       as the accuracy for rotational invariance seems to deteriorate
#       with increasing \omega.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
# bash fam_rotational_invariance.sh [EXESUFFIX] [PARAM] [PAIRING] [OMEGA] [ETA] [-v/--verbose] [--free]
#
# where
#  - EXESUFFIX specifies the executable to be used
#  - PARAM the parameterisation.
#  - PAIRING = 'HF' or 'HFB'
#  - OMEGA   = the frequency [MeV]
#  - ETA     = the smearing  [MeV]
#
# Calling the script with flags
#    -v or --verbose will print the values which are compared and details on symmetry
#    --free will execute the test for the free response, i.e. without induced mean-fields
# Dependencies:
#   None
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner   : W. Ryssens [wouter.ryssens@ulb.be]
#--------------------------------------------------------------------------------
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh
LC_NUMERIC="en_US.UTF-8" # Important to deal with floats and particularly with scientific notation

# Initialize verbose mode as false by default
verbose=false
free=false

# Temporary array to hold arguments
args=()

# Parse all arguments for -v or --verbose
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--verbose)
      verbose=true
      shift
      ;;
    --free)
      free=true
      shift
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done


# Set up FAM testing environment
setup_test_env_fam "fam_rotational_invariance.${args[2]}.om=${args[3]}.eta=${args[4]}" "${args[0]}" "${args[0]}" "${args[1]}"

mfoutfile=../logs/fam_rotational_invariance.${args[2]}.om=${args[3]}.eta=${args[4]}.${args[0]}.${args[1]}.out

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1a) Run a standard mean-field calculation
# Create runtime data
cat << EOF > mf.data
&nucleus
neutrons=20, protons=22
energy_prec=1e-16
/
&mesh
nx=8, ny=8, nz=8, dx=1.0
/
&func
name_param="${args[1]}"
/
&pairing
type="${args[2]}"
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 40, nwp = 40
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='mf_hfb.wf'
allowtransform=.true.
/
&MomentParam
moreconstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=0
moreconstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=0
/
&Cranking
/
EOF

# Run the calculation
./$exe < mf.data > $mfoutfile
# .... and immediately check if MOCCa reported back some error codes
mocca_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2b) Run the LO QFAM calculations

for m in 0 2
do

# Create runtime data
cat << EOF > qfam.data
&nucleus
neutrons=20, protons=22
/
&mesh
nx=8, ny=8, nz=8, dx=1.0
/
&func
name_param="${args[1]}"
/
&pairing
type="${args[2]}"
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 40, nwp = 40
/
&IO
InputFilename='mf_hfb.wf'
OutputFilename='trash'
famfile="S_2$m.fam"
/
&MomentParam
/
&Cranking
/
&fam
omega=$omega${args[3]}
smear=$smear${args[4]}
l=2
m=$m
maxiter=$(if [ "$free" = true ]; then echo "0"; else echo "30"; fi)
fam_precision=1e-8
/
EOF

#Redefine famoutfile to keep all logs
famoutfile=../logs/fam_rotational_invariance.${args[2]}.om=${args[3]}.eta=${args[4]}.${args[0]}.${args[1]}.fam.Q2$m.out

# Run the calculation
./$exefam < qfam.data > $famoutfile
done

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3) Parse and compare the FAM output files
# Extract strength components from both files
# Each returns: real_part imaginary_part
components_Q20=($(extract_strength_components "S_20.fam"))
components_Q2p2=($(extract_strength_components "S_22.fam"))

# Assign to readable variable names
Q20_re=${components_Q20[0]}
Q20_im=${components_Q20[1]}
Q2p2_re=${components_Q2p2[0]}
Q2p2_im=${components_Q2p2[1]}

# Verbose output
if [ "$verbose" = true ]; then
  echo ""
  echo "Strength function components (real, imaginary):"
  echo "  Q20 : $Q20_re  $Q20_im"
  echo "  Q2+2: $Q2p2_re  $Q2p2_im"
  echo ""
fi

tolerance=0.15

# Check rotational invariance: Q20 == Q2+2 within tolerance
all_pass=true

# Compare Q20 vs Q2+2
compare_floats $Q20_re $Q2p2_re $tolerance
real_check=$?

compare_floats $Q20_im $Q2p2_im $tolerance
imaginary_check=$?


if [[ $real_check -ne 0 ]] ; then
  echo "FAIL: Real parts of Q20 ($Q20_re) and Q2+2 ($Q2p2_re) differ by more than $tolerance"
  all_pass=false
fi
if [[ $imaginary_check -ne 0 ]] ; then
  echo "FAIL: Imaginary parts of Q20 ($Q20_im) and Q2+2 ($Q2p2_im) differ by more than $tolerance"
  all_pass=false
fi

if [ "$all_pass" = true ]; then
  echo "PASS: Both components (Q20, Q2+2) are equal within tolerance of $tolerance"
fi

# Exit with error if any comparison failed
if [ "$all_pass" = false ]; then
  exit 1
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Clean up the test environment
teardown_test_env
