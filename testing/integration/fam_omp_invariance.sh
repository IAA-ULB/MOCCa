#!/usr/bin/env sh
#--------------------------------------------------------------------------------
# Verify that FAMQRPA calculations are invariant under serial and OpenMP execution
# by comparing the strength functions from serial, 1-thread and 4-thread executions.
#--------------------------------------------------------------------------------
# This script tests:
#
#   |(complex) Strength_serial - Strength_omp_1thread| < | TOLERANCE |
#   |(complex) Strength_serial - Strength_omp_4threads| < | TOLERANCE |
#
# for a spherically symmetric configuration of Ti42 in a small modelspace.
# The same mean-field wavefunction is used as input for all FAM calculations
# to isolate the linear response part.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
# bash fam_omp_invariance.sh [EXESUFFIX] [PARAM] [PAIRING] [OMEGA] [ETA] [OMPSUFFIX] [-v/--verbose] [--free]
#
# where
#  - EXESUFFIX specifies the serial executable suffix to be used
#  - PARAM the parameterisation.
#  - PAIRING = 'HF' or 'HFB'
#  - OMEGA   = the frequency [MeV]
#  - ETA     = the smearing  [MeV]
#  - OMPSUFFIX = the OpenMP executable suffix (optional, defaults to EXESUFFIX-OMP)
#
# Calling the script with flags
#    -v or --verbose will print the values which are compared
#    --free will execute the test for the free response, i.e. without induced mean-fields
# Dependencies:
#   None
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner   : W. Ryssens [wouter.ryssens@ulb.be]
#--------------------------------------------------------------------------------
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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

# Set default OMP suffix if not provided
# Args: EXESUFFIX, PARAM, PAIRING, OMEGA, ETA, [OMPSUFFIX]
if [ ${#args[@]} -ge 6 ]; then
  omp_suffix="${args[5]}"
else
  omp_suffix="${args[0]}-OMP"
fi

# Tolerance for comparing results
tolerance=0.001

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1) Run a standard mean-field calculation
# We use the serial executable for the mean-field calculation

# Set up mean-field test environment
setup_test_env "mf" "${args[0]}" "${args[1]}"

mfoutfile=../logs/fam_omp_invariance.${args[2]}.om=${args[3]}.eta=${args[4]}.${args[0]}.${args[1]}.mf.out

# Create runtime data for mean-field
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

# Run the mean-field calculation
./$exe < mf.data > $mfoutfile
mocca_check_mf=$?

# Save the mean-field wavefunction for FAM calculations
cp mf_hfb.wf ../mf_hfb.wf

# Clean up mean-field environment
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2) Run the LO QFAM calculations with serial executable

# Set up serial FAM environment
setup_test_env_fam "fam_omp_serial" "${args[0]}" "${args[0]}" "${args[1]}"

omega=${args[3]}
smear=${args[4]}

# Create runtime data for FAM
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
InputFilename='../mf_hfb.wf'
OutputFilename='trash'
famfile="S_20_serial.fam"
/
&MomentParam
/
&Cranking
/
&fam
omega=$omega
smear=$smear
l=2
m=0
maxiter=$(if [ "$free" = true ]; then echo "0"; else echo "30"; fi)
fam_precision=1e-8
/
EOF

famoutfile_serial=../logs/fam_omp_invariance.${args[2]}.om=${args[3]}.eta=${args[4]}.${args[0]}.${args[1]}.serial.fam.out
./$exefam < qfam.data > $famoutfile_serial
mocca_check_fam_serial=$?

# Extract strength components from all three files
# Each returns: real_part imaginary_part
components_serial=($(extract_strength_components "S_20_serial.fam"))

# Clean up serial FAM environment
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3) Run the LO QFAM calculations with OpenMP executable using 1 thread

# Set up OpenMP FAM environment with 1 thread
setup_test_env_fam "fam_omp_1thread" "${args[0]}" "$omp_suffix" "${args[1]}"

# Create runtime data for FAM with 1 thread
export OMP_NUM_THREADS=1

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
InputFilename='../mf_hfb.wf'
OutputFilename='trash'
famfile="S_20_omp1.fam"
/
&MomentParam
/
&Cranking
/
&fam
omega=$omega
smear=$smear
l=2
m=0
maxiter=$(if [ "$free" = true ]; then echo "0"; else echo "30"; fi)
fam_precision=1e-8
/
EOF

famoutfile_omp1=../logs/fam_omp_invariance.${args[2]}.om=${args[3]}.eta=${args[4]}.${omp_suffix}.${args[1]}.omp1.fam.out
./$exefam < qfam.data > $famoutfile_omp1
mocca_check_fam_omp1=$?

components_omp1=($(extract_strength_components "S_20_omp1.fam"))

# Clean up OpenMP 1-thread FAM environment
unset OMP_NUM_THREADS
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (4) Run the LO QFAM calculations with OpenMP executable using 4 threads

# Set up OpenMP FAM environment with 4 threads
setup_test_env_fam "fam_omp_4threads" "${args[0]}" "$omp_suffix" "${args[1]}"

# Create runtime data for FAM with 4 threads
export OMP_NUM_THREADS=4

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
InputFilename='../mf_hfb.wf'
OutputFilename='trash'
famfile="S_20_omp4.fam"
/
&MomentParam
/
&Cranking
/
&fam
omega=$omega
smear=$smear
l=2
m=0
maxiter=$(if [ "$free" = true ]; then echo "0"; else echo "30"; fi)
fam_precision=1e-8
/
EOF

famoutfile_omp4=../logs/fam_omp_invariance.${args[2]}.om=${args[3]}.eta=${args[4]}.${omp_suffix}.${args[1]}.omp4.fam.out
./$exefam < qfam.data > $famoutfile_omp4
mocca_check_fam_omp4=$?

components_omp4=($(extract_strength_components "S_20_omp4.fam"))

# Clean up OpenMP 4-thread FAM environment
unset OMP_NUM_THREADS
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (5) Parse and compare the FAM output files

# Assign to readable variable names
serial_re=${components_serial[0]}
serial_im=${components_serial[1]}
omp1_re=${components_omp1[0]}
omp1_im=${components_omp1[1]}
omp4_re=${components_omp4[0]}
omp4_im=${components_omp4[1]}

# Verbose output
if [ "$verbose" = true ]; then
  echo ""
  echo "Strength function components (real, imaginary):"
  echo "  Serial     : $serial_re  $serial_im"
  echo "  OpenMP 1T  : $omp1_re  $omp1_im"
  echo "  OpenMP 4T  : $omp4_re  $omp4_im"
  echo ""
fi

# Check OpenMP invariance: serial == OpenMP (1 thread and 4 threads) within tolerance
all_pass=true

# Compare serial vs OpenMP 1 thread
compare_floats $serial_re $omp1_re $tolerance
real_check_omp1=$?
compare_floats $serial_im $omp1_im $tolerance
imag_check_omp1=$?

# Compare serial vs OpenMP 4 threads
compare_floats $serial_re $omp4_re $tolerance
real_check_omp4=$?
compare_floats $serial_im $omp4_im $tolerance
imag_check_omp4=$?

# Check runtime errors
runtime_check=$(( $mocca_check_mf || $mocca_check_fam_serial || $mocca_check_fam_omp1 || $mocca_check_fam_omp4 ))

if [[ $real_check_omp1 -ne 0 ]] ; then
  echo "FAIL: Real parts of serial ($serial_re) and OpenMP 1T ($omp1_re) differ by more than $tolerance"
  all_pass=false
fi
if [[ $imag_check_omp1 -ne 0 ]] ; then
  echo "FAIL: Imaginary parts of serial ($serial_im) and OpenMP 1T ($omp1_im) differ by more than $tolerance"
  all_pass=false
fi
if [[ $real_check_omp4 -ne 0 ]] ; then
  echo "FAIL: Real parts of serial ($serial_re) and OpenMP 4T ($omp4_re) differ by more than $tolerance"
  all_pass=false
fi
if [[ $imag_check_omp4 -ne 0 ]] ; then
  echo "FAIL: Imaginary parts of serial ($serial_im) and OpenMP 4T ($omp4_im) differ by more than $tolerance"
  all_pass=false
fi
if [[ $runtime_check -ne 0 ]] ; then
  echo "FAIL: One or more calculations reported runtime errors"
  all_pass=false
fi

if [ "$all_pass" = true ]; then
  echo "PASS: Serial and OpenMP results (1T and 4T) are all equal within tolerance of $tolerance"
fi

# Clean up saved files
rm -f  ../mf_hfb.wf

# Exit with error if any comparison failed
exitcode=$(( $runtime_check || $real_check_omp1 || $imag_check_omp1 || $real_check_omp4 || $imag_check_omp4 ))
exit $exitcode
