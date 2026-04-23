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
#   bash fam_unit_test.sh [EXESUFFIX] [--pairing HF|HFB] [--parameterisation PARAM] [--nw NW] [-v/--verbose]
#
# where
# - EXESUFFIX     : maximally symmetric executable
# -- pairing      : specifies the pairing type (HF or HFB, default: HFB)
# --parameterisation : specifies the parameterization (default: t0t3)
# --nw           : specifies the number of wavefunctions (default: 20)
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
nw=20  # Default value for number of wavefunctions

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
    --nw)
      shift
      if [[ "$1" =~ ^[0-9]+$ ]]; then
        nw="$1"
      else
        echo "Error: --nw must be a positive integer"
        exit 1
      fi
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
setup_test_env_fam "fam_units_tests" "${args[0]}" "${args[0]}" "$parameterisation"
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1) Run the mean-field calculation

# Create runtime data
cat << EOF > mf.data
&nucleus
neutrons=10, protons=8
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
nwn = $nw, nwp = $nw
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
neutrons=10, protons=8
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
nwn = $nw, nwp = $nw
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
# (3) Run the unit tests

# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=10, protons=8
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
nwn = $nw, nwp = $nw
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
unit_test=.true.
/
EOF

# Run the calculation
./$exefam < fam.data | tee $famoutfile
# .... and immediately check if Tantalus reported back some error codes
fam_check=$?

teardown_test_env
