#!/usr/bin/env sh
#-------------------------------------------------------------------------------
# Run the unit tests of the FAM solver for a 18O nucleus in a minimal box.
#  These tests (and their success conditions) are defined in the Fortran code.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_unit_test.sh [EXESUFFIX] [--pairing HF|HFB] [--parameterisation PARAM] [--nw NW]
#
# where
# EXESUFFIX           : executable to use
# -- pairing          : specifies the pairing type (HF or HFB, default: HFB)
# -- parameterisation : specifies the parameterization (default: t0t3)
# -- nw               : specifies the number of wavefunctions (default: 15)
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: c44b9e6a4b6173a87afcce8c73a175bb05e12182
#-------------------------------------------------------------------------------
# Initialize verbose mode as false by default
pairing="HFB"  # Default value
parameterisation="t0t3"  # Default value
nw=15  # Default value for number of wavefunctions

# Temporary array to hold arguments
args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
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
setup_test_env_fam "fam_unit_tests" "${args[0]}" "${args[0]}" "$parameterisation"
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
if [ $tantalus_check -ne 0 ]; then
    echo "ERROR: Mean-field calculation failed with exit status $tantalus_check"
fi

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
if [ $fam_check -ne 0 ]; then
    echo "ERROR: FAM calculation failed with exit status $fam_check"
fi

teardown_test_env

# Return overall exit status - fail if any calculation failed
exit_status=0
if [ $tantalus_check -ne 0 ] || [ $fam_check -ne 0 ]; then
    exit_status=1
fi

echo ""
echo "========= Test Summary ========"
echo "Mean-field exit status : $tantalus_check"
echo "FAM unit tests status  : $fam_check"
echo "Overall exit status    : $exit_status"
echo "==============================="

exit $exit_status
