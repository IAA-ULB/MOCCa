#!/usr/bin/env sh
#-------------------------------------------------------------------------------
# Perform HFB + linear response calculations of Ti42 with t0t3 and check
# that the proton number particle operator is the momentum of a zero-mode.
# Fits a zero-mode model and verifies the fitted frequency components are small.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  Zero-mode omega_r (popt[1])            0 MeV                      0.1 MeV
#  Zero-mode omega_i (popt[2])            0 MeV                      0.1 MeV
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_pairing_zero_mode.sh [EXESUFFIX] [--pairing HF|HFB] [--parameterisation PARAM] [-v/--verbose] [--nx NX] [--nw NW]
#
# where
#   EXESUFFIX        : executable to test
# --pairing          : specifies the pairing type (HF or HFB, default: HFB)
# --parameterisation : specifies the parameterization (default: t0t3)
# --nx               : number of mesh points (in all 3 dimensions)
# --nw               : number of single-particle wavefunctions (of one species)
# -v or --verbose    : print the values which are compared.
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash:
#-------------------------------------------------------------------------------

# - - - - - - - - - - - - - -
# Default values for options
verbose=false
pairing="HFB"
parameterisation="t0t3"
nx=8
nw=30
# - - - - - - - - - - - - - -

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
    --nx)
      shift
      nx="$1"
      shift
      ;;
    --nw)
      shift
      nw="$1"
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
nx=$nx, ny=$nx, nz=$nx, dx=1.0
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

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2) Run the LO FAM calculation

# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=20, protons=22
/
&mesh
nx=$nx,ny=$nx,nz=$nx, dx=1.0
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
famfile='N.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega_min=0.0
omega_max=+0.1
omega_step=0.01
smear=0.0
!l=0
!m=0
operator_type='particle number'
maxiter=30
maxhist=31
eff_charge_n=0.0
/
EOF

# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if MOCCa reported back some error codes
fam_check=$?


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (4) Check that strength at non-zero frequencies is quite a lot smaller
#     than the one at zero frequency
#
# Quick inline python script
cat << EOF > analyse.py
import numpy as np
from scipy.optimize import curve_fit

def zero_mode(omega, M, omega_ng_r, omega_ng_i):
    omega_ng_sq = (omega_ng_r**2 - omega_ng_i**2) + 2j * omega_ng_r * omega_ng_i
    denominator = omega**2 - omega_ng_sq
    # Avoid division by zero near pole
    with np.errstate(divide="ignore", invalid="ignore"):
        f = M * omega_ng_sq / denominator
    return f.real

dat = np.loadtxt("N.fam")
popt, __ = curve_fit(zero_mode, dat[:, 0], -dat[:, 3], p0=[dat[0,3], 0.001, 0.01])

if np.abs(popt[1]) > 0.1 or np.abs(popt[2]) > 0.1:
    ifail = 1
else:
    ifail = 0

# Print results for bash to capture: omega_r omega_i ifail
print(f"{popt[1]} {popt[2]} {ifail}")
EOF

# Run python script and capture all three outputs
read omega_r omega_i ifail <<< $(python analyse.py)

# Report results
if $verbose; then
  echo "----------------------------------------"
  echo "Zero-mode strength check:"
  echo "Fitted omega_r: $omega_r MeV"
  echo "Fitted omega_i: $omega_i MeV"
  echo "Failure code (ifail): $ifail"
  if (($ifail == 0)) ; then
    echo "Result: PASS - Zero-mode frequency components are small (< 0.1 MeV)"
  else
    echo "Result: FAIL - Zero-mode frequency components are too large (>= 0.1 MeV)"
  fi
  echo "----------------------------------------"
fi

# Clean up
#teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
fail=$(($mocca_check || $fam_check || $ifail))

if (($fail == 0)) ; then
  echo -e "test FAM pairing zero mode :\033[1;32m success \033[0m"
else
  echo -e "test FAM pairing zero mode :\033[1;31m failed ! exit status : tant = $mocca_check, fam = $fam_check, zero_mode = $ifail \033[0m"
fi

exit $fail
