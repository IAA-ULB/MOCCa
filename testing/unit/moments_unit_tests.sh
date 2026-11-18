#!/usr/bin/env sh
#-------------------------------------------------------------------------------
# Run the moments unit tests for an O16 nucleus in a minimal box.
#  These tests (and their success conditions) are defined in the Fortran code
#  of the moments_testing module: a dense dump of all multipole matrices with
#  l <= 4 and the trace check of those matrices against the many-body moment
#  values calculated in moments.f90.
#
#  The unit tests are armed with unit_test=.true. in the &mesh/ namelist: the
#  code runs a short Hartree-Fock evolution for O16, executes the tests on the
#  final mean-field state and stops with the exit status of the test suite.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash moments_unit_tests.sh [EXESUFFIX] [--pairing HF|BCS|HFB]
#
# where
# EXESUFFIX       : executable to use, e.g. "NLO" for MOCCa.NLO.exe
#                   (default: NLO)
# --pairing       : pairing type of the mean-field state (default: HF).
#                   Note that the trace check is skipped for HFB, whose
#                   density matrix is diagonal in the canonical basis.
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: this script lands with the moments_testing module
#-------------------------------------------------------------------------------
# Initialize default values
suffix="NLO"   # Default executable suffix
pairing="HF"   # Default pairing type

# Temporary array to hold arguments
args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pairing)
      shift
      if [[ "$1" == "HF" || "$1" == "BCS" || "$1" == "HFB" ]]; then
        pairing="$1"
      else
        echo "Error: --pairing must be HF, BCS or HFB"
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

# Restore the positional arguments (the executable suffix)
set -- "${args[@]}"
if [[ "$1" != "" ]]; then
  suffix="$1"
fi

#- - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env "moments_unit_tests" "$suffix" "SLy4"

# Create runtime data
cat << EOF > mocca.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
unit_test=.true.
/
&func
name_param='SLy4'
/
&pairing
type="$pairing"
/
&evolution
maxiter=25
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='trash'
/
&MomentParam
/
&Cranking
/
EOF

# Run the calculation: the code exits with the status of the test suite
./$exe < mocca.data > $outfile
mocca_check=$?

teardown_test_env

#- - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - -
echo ""
echo "========= Test Summary ========"
echo "Moments unit tests status : $mocca_check"
echo "Overall exit status        : $mocca_check"
echo "Full output at testing/unit/logs/moments_unit_tests.$suffix.out"
echo "==============================="

exit $mocca_check
