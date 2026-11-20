#!/usr/bin/env bash
#-------------------------------------------------------------------------------
# Run the unit tests of the symmetry-block matrix module, phases 1 and 2
#  plus the basis rotations of design section 6a: construction of the
#  block list from MOCCa's HFBlocks bookkeeping, population of the stored
#  blocks through the fill callbacks, the partner-reconstructing dense
#  expansion of the design's section 13, and the unitary (sector-diagonal
#  real orthogonal) transformations of stored matrices - plus the
#  honest-trace contract of D-7 (design 15.8 amendment, 2026-11-20).
#
# The tests are standalone programs that are compiled on the fly together
#  with the module itself (which has no $-markers and no dependencies on
#  other MOCCa modules), then executed and the exit statuses reported:
#   *) block_matrix_unit_tests    : the in-process check suite
#   *) block_matrix_refusal_tests : the rotation structure guards, which
#      stop the process - every scenario is driven separately and must
#      exit with status 1
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

gfortran -O2 -J "${work}" -o "${work}/block_matrix_refusal_tests.exe" \
  "${sdir}/src_orig/block_matrix_module.f90" \
  "${udir}/block_matrix_refusal_tests.f90"

# Every refusal scenario must be stopped by the module's guards: exit
#  status 1 = refused as expected, anything else is a failure.
refusal_failed=0
scenario=1
while [ "${scenario}" -le 9 ]; do
  rc=0
  "${work}/block_matrix_refusal_tests.exe" "${scenario}" || rc=$?
  if [ "${rc}" -ne 1 ]; then
    echo "refusal scenario ${scenario}: expected exit 1, got ${rc}: FAIL"
    refusal_failed=1
  else
    echo "refusal scenario ${scenario}: refused as expected"
  fi
  scenario=$((scenario + 1))
done

if [ "${refusal_failed}" -ne 0 ]; then
  echo "block_matrix refusal tests FAILED"
  exit 1
fi
echo "All block_matrix refusal tests passed."
