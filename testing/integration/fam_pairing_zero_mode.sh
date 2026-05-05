#!/usr/bin/env sh
#-------------------------------------------------------------------------------
# Perform HFB + linear response calculations of Ti42 with t0t3 and check
# that the proton number particle operator is the momentum of a zero-mode, i.e.
# that the
#     S(N_p, \omega) = 0      if \omega != 0
#                    = - M_n  if \omega  = 0
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_pairing_zero_mode.sh [EXESUFFIX] [--pairing HF|HFB] [--parameterisation PARAM] [-v/--verbose]
#
# where
# - EXESUFFIX     : maximally symmetric executable
# -- pairing      : specifies the pairing type (HF or HFB, default: HFB)
# --parameterisation : specifies the parameterization (default: t0t3)
# -v or --verbose : print the values which are compared.
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash:
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
setup_test_env_fam "fam_pairing_zero_mode" "${args[0]}" "${args[0]}" "$parameterisation"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1) Run the mean-field calculation

# Create runtime data
cat << EOF > mf.data
&nucleus
neutrons=20, protons=22
energy_prec=1e-16
/
&mesh
nx=12, ny=12, nz=12, dx=1.2
/
&func
name_param='$parameterisation'
/
&pairing
type="$pairing"
/
&evolution
maxiter=1000
!dt=0.0180, momentum=0.6
!Estimateparams=.false.
/
&scfiteration
/
&wfs
nwn = 30, nwp = 30
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
nx=12, ny=12, nz=12, dx=1.2
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
nwn = 30, nwp = 30
osc_freq = 0.2, 0.2, 0.18
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
nx=12, ny=12, nz=12, dx=1.2
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
nwn = 30, nwp = 30
/
&IO
InputFilename='mf.wf'
OutputFilename='trash'
allowtransform=.true.
famfile='N.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega_min=-1.0
omega_max=1.0
omega_step=0.1
smear=0.00
!l=0
!m=0
operator_type='particle number'
maxiter=30
fam_precision=1e-8
eff_charge_n=0.0
/
EOF

# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
fam_check=$?
