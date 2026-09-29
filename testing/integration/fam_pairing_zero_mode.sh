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
export OMP_NUM_THREADS=4
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

# ------------------------------------------------------------------
# (4) Fit a zero-mode model and verify the fitted frequency components
cat <<'PYEOF' > analyse.py
import numpy as np
from scipy.optimize import curve_fit

# Optional plotting for manual debugging
#import matplotlib
#matplotlib.use('Agg')
#import matplotlib.pyplot as plt

def zero_mode(omega, M, omega_ng_r, omega_ng_i):
    omega_ng_sq = (omega_ng_r**2 - omega_ng_i**2) + 2j * omega_ng_r * omega_ng_i
    denominator = omega**2 - omega_ng_sq
    with np.errstate(divide='ignore', invalid='ignore'):
        f = M * omega_ng_sq / denominator
    return f.real

# Load data: omega, ..., strength (4th column = index 3)
dat = np.loadtxt('N.fam')
omega = dat[:, 0]
strength = -dat[:, 3]          # negative because we fit -dat[:,3]

# --- Robust fitting strategy -------------------------------------------------
#  Multiple initial guesses, including (ω_r, ω_i) = (0, 0)
initial_guesses = [
    [strength[0], 0.0,      0.0],       # actual zero mode
    [strength[0], 0.001,    0.0],       # small real
    [strength[0], 0.0,      0.01],      # small imag
    [strength[0], -0.001,   0.0],       # negative real
    [strength[0], 0.0,      -0.01],     # negative imag
    [strength[0], 0.001,    +0.01],     # small both
]

best_popt = None
best_residual = float('inf')

for p0 in initial_guesses:
    try:
        popt, pcov = curve_fit(
            zero_mode, omega, strength,
            p0=p0,  maxfev=10000
        )
        residual = np.sum((zero_mode(omega, *popt) - strength) ** 2)
        if residual < best_residual:
            best_residual = residual
            best_popt = popt
    except (RuntimeError, TypeError, ValueError):
        continue  # Try next guess

# --- Fallback if all guesses fail ------------------------------------------
#if best_popt is None:
#    # Use a safe default that will fail the tolerance check
#    best_popt = [strength[0], 1.0, 1.0]

# --- Draw plot --------------------------------------------------------------
# Generate fitted curve
#omega_fit = np.linspace(min(omega), max(omega), 500)
#strength_fit = zero_mode(omega_fit, *best_popt)

#plt.figure(figsize=(8, 6))
#plt.scatter(omega, strength, color='blue', label='FAM strength', s=20, zorder=5)
#plt.plot(omega_fit, strength_fit, color='red', linewidth=2,
#         label=f'Fit: M={best_popt[0]:.2f}, $\omega_{{r}}$={best_popt[1]:.4f}, $\omega_{{i}}$={best_popt[2]:.4f}')
#plt.axvline(x=0.0, color='black', linestyle=':', linewidth=1)
#plt.axhline(y=0.0, color='black', linestyle=':', linewidth=1)
#plt.xlabel('$\omega$ (MeV)')
#plt.ylabel('Strength')
#plt.title('Zero-mode fit for particle number operator')
#plt.grid(True, alpha=0.3)
#plt.legend()
#plt.tight_layout()
#plt.savefig('zero_mode_fit.png', dpi=150)
#plt.close()

# --- Decide pass/fail --------------------------------------------------------
ifail = 0 if (abs(best_popt[1]) <= 0.1 and abs(best_popt[2]) <= 0.1) else 1

# Output: one value per line (ω_r, ω_i, ifail)
print(best_popt[1])
print(best_popt[2])
print(ifail)
PYEOF

# Run analysis and *only* capture stdout; warnings go to the terminal.
if ! python3 analyse.py > analyse.out; then
    echo "Error: Python analysis script failed"
    ifail=1
    omega_r=""
    omega_i=""
else
    # Read the three lines explicitly
    omega_r=$(head -n1 analyse.out)
    omega_i=$(sed -n '2p' analyse.out)
    ifail=$(sed -n '3p' analyse.out)
    # Sanity: if any line is missing we treat as failure
    if [ -z "$omega_r" ] || [ -z "$omega_i" ] || [ -z "$ifail" ]; then
        echo "Error: Python script produced incomplete output"
        echo "Content of analyse.out:"
        cat analyse.out
        ifail=1
    fi
fi

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

# Clean up - keeping N.fam and plot!
cp N.fam ../logs/fam_pairing_zero_mode.${args[0]}.N.fam
cp zero_mode_fit.png ../logs/fam_pairing_zero_mode.${args[0]}.zero_mode_fit.png
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
