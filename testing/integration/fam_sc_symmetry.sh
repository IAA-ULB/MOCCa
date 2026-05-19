#!/usr/bin/env sh
#-------------------------------------------------------------------------------
# Perform mean-field + linear response calculations of Ti22 with t0t3 in a
# minimal box with different self-consistent symmetry options and check that
# the resulting monopole strengths are identical.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  Monopole strength                     [result of EXE-T]          1e-6
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_symmetry.sh [EXESUFFIX] [EXESUFFIX_T] [EXESUFFIX_P] [--pairing HF|HFB] [--parameterisation PARAM] [-v/--verbose]
#
# where
# - EXESUFFIX     : maximally symmetric executable
# - EXESUFFIX_T   : time-reversal breaking executable
# - EXESUFFIX_P   : parity breaking executable
# -- pairing      : specifies the pairing type (HF or HFB, default: HFB)
# --parameterisation : specifies the parameterization (default: t0t3)
# -v or --verbose : print the values which are compared.
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: c44b9e6a4b6173a87afcce8c73a175bb05e12182
#-------------------------------------------------------------------------------
# Initialize verbose mode as false by default
verbose=false
pairing="HFB"  # Default value
parameterisation="t0t3"  # Default value

# Temporary array to hold arguments
args=()

# Parse all arguments for -v or --verbose or --pairing
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--verbose)
      verbose=true
      shift
      ;;
    --pairing)
      shift
      if [[ "$1" == "HF" || "$1" == "HFB" ]]; then
        pairing="$1"
      else
        echo "Error: --pairing must be either 'HF' or 'HFB'"
        exit 1
      fi
      shift
      ;;
    --parameterisation)
      shift
      parameterisation="$1"
      shift
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh


# Set up
setup_test_env_fam "fam_sc_symmetry.$pairing" "${args[0]}" "${args[0]}" "$parameterisation"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1) Run the mean-field calculation

# Create runtime data
cat << EOF > mf.data
&nucleus
neutrons=20, protons=22
energy_prec=1e-16
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='$parameterisation'
/
&pairing
type="$pairing"
/
&evolution
maxiter=100
dt=0.0209, momentum=0.5746
Estimateparams=.false.
/
&scfiteration
/
&wfs
nwn = 20, nwp = 20
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
# .... and immediately check if Tantalus reported back some error codes
tantalus_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2) Run a calculation that freezes the potentials just to get a
#     robust set of virtual states
cat << EOF > mf.data
&nucleus
neutrons=20, protons=22
energy_prec=1e-20
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='$parameterisation'
/
&pairing
type="$pairing"
/
&evolution
maxiter=1000
dt=0.0209, momentum=0.5746
Estimateparams=.false.
freezeiter=1000
/
&scfiteration
/
&wfs
nwn = 20, nwp = 20
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='mf.wf'
OutputFilename='mf.wf'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
EOF
./$exe < mf.data > $mfoutfile.bis

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3) Run the LO FAM calculation

# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=20, protons=22
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='$parameterisation'
/
&pairing
type="$pairing"
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 20, nwp = 20
/
&IO
InputFilename='mf.wf'
OutputFilename='trash'
allowtransform=.true.
famfile='S_20.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega=25.0
smear=1.0
l=0
m=0
maxiter=30
fam_precision=1e-8
mixingscheme=0
/
EOF

# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
fam_check=$?

# We keep the wf file!
mv mf.wf ../
# Extract strength value before cleanup
omega=25.0
strength_fam=$(get_strength "S_20.fam" "$omega")
# ... but otherwise clean-up
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (4) Run the LO-T FAM calculation
setup_test_env_fam "fam_sc_symmetry.$pairing" "${args[0]}" "${args[1]}" "$parameterisation"
cp ../mf.wf .

# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=20, protons=22
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='$parameterisation'
/
&pairing
type="$pairing"
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
InputFilename='mf.wf'
OutputFilename='trash'
allowtransform=.true.
famfile='S_20.T.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega=25.0
smear=1.0
l=0
m=0
maxiter=30
fam_precision=1e-8
mixingscheme=0
/
EOF
# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
fam_T_check=$?

# Extract strength value before cleanup
strength_fam_t=$(get_strength "S_20.T.fam" "$omega")
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (4) Run the LO-P FAM calculation
setup_test_env_fam "fam_sc_symmetry.$pairing" "${args[0]}" "${args[2]}" "$parameterisation"
cp ../mf.wf .

# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=20, protons=22
/
&mesh
nx=8, ny=8, nz=16, dx=0.8
/
&func
name_param='$parameterisation'
/
&pairing
type="$pairing"
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 20, nwp = 20
/
&IO
InputFilename='mf.wf'
OutputFilename='trash'
allowtransform=.true.
famfile='S_20.P.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega=25.0
smear=1.0
l=0
m=0
maxiter=30
fam_precision=1e-10
mixingscheme=0
/
EOF
# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
fam_P_check=$?

# Extract strength value before cleanup
strength_fam_p=$(get_strength "S_20.P.fam" "$omega")
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (5) Parse and compare the FAM output files
# Define tolerance for comparison
tolerance=0.0000001

# Print values if verbose mode is enabled
if [ "$verbose" = true ]; then
    echo "Strength function values at omega = $omega:"
    LC_NUMERIC="en_US.UTF-8" printf "QFAM   strength: %12.8f\n" $strength_fam
    LC_NUMERIC="en_US.UTF-8" printf "QFAM-T strength: %12.8f (reference)\n" $strength_fam_t
    LC_NUMERIC="en_US.UTF-8" printf "QFAM-P strength: %12.8f\n" $strength_fam_p
    echo ""
fi

# Compare FAM vs FAM-T (treat FAM-T as correct)
compare_floats "$strength_fam" "$strength_fam_t" "$tolerance"
fam_vs_famt=$?

# Compare FAM-P vs FAM-T (treat FAM-T as correct)
compare_floats "$strength_fam_p" "$strength_fam_t" "$tolerance"
fam_p_vs_famt=$?

# Overall test results
if [ $fam_vs_famt -eq 0 ] && [ $fam_p_vs_famt -eq 0 ]; then
    echo "✓ ALL COMPARISONS PASSED"
    echo "  FAM strength matches FAM-T reference"
    echo "  FAM-P strength matches FAM-T reference"
    exit 0
else
    echo "✗ SOME COMPARISONS FAILED"
    if [ $fam_vs_famt -ne 0 ]; then
        echo "  FAM strength does NOT match FAM-T reference"
        LC_NUMERIC="en_US.UTF-8" printf "  FAM: %.8f vs FAM-T: %.8f\n" $strength_fam $strength_fam_t
    fi
    if [ $fam_p_vs_famt -ne 0 ]; then
        echo "  FAM-P strength does NOT match FAM-T reference"
        LC_NUMERIC="en_US.UTF-8" printf "  FAM-P: %.8f vs FAM-T: %.8f\n" $strength_fam_p $strength_fam_t
    fi
    exit 1
fi
