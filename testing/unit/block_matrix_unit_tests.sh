#!/usr/bin/env bash
#-------------------------------------------------------------------------------
# Run the unit tests of the symmetry-block matrix module, phase 1:
#  construction of the block list from MOCCa's HFBlocks bookkeeping.
#
# The tests are a standalone program that is compiled on the fly together
#  with the module itself (which has no $-markers and no dependencies on
#  other MOCCa modules), then executed and its exit status reported.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash block_matrix_unit_tests.sh
#
# Dependencies: gfortran; no MOCCa executable needs to be built.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
#-------------------------------------------------------------------------------

set -e

# The script is meant to be run from inside the testing/unit directory, but
#  do not rely on it: resolve the directories from the script's location.
udir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
sdir="$(cd "${udir}/../.." && pwd)"

work="$(mktemp -d)"
trap 'rm -rf "${work}"' EXIT

gfortran -O2 -J "${work}" -o "${work}/block_matrix_unit_tests.exe" \
  "${sdir}/src_orig/block_matrix_module.f90" \
  "${udir}/block_matrix_unit_tests.f90"

"${work}/block_matrix_unit_tests.exe"
