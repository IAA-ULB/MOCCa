#!/usr/bin/env bash
#-------------------------------------------------------------------------------
# Unit tests of spectrum_count_module: count the eigenvalues of a small mesh
# hamiltonian by Chebyshev trace estimation and compare against the exact
# count from a dense diagonalization. Two tiers:
#   - 8^3 smoke test (fast, suitable for CI)
#   - 16^3 reference tier (dense level spectrum near the cutoff; run locally)
# Usage:
#   bash chebyshev_count_unit_tests.sh [8|16]
#-------------------------------------------------------------------------------
MESH=${1:-8}
if [ "$MESH" != "8" ] && [ "$MESH" != "16" ]; then
  echo "Usage: bash chebyshev_count_unit_tests.sh [8|16]"
  exit 1
fi

# Run from inside this directory, per the testing conventions.
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

# Transient work directory, destroyed at teardown
WORK=$PWD/work/chebyshev_count.$MESH
rm -rf $WORK ; mkdir -p $WORK/src

# Push the marker-bearing sources through Hephaestos; the counter module is
# marker-free but goes through the same pass for uniformity. Hephaestos
# resolves configs/ relative to the repository root.
cd $ROOT
for f in compilation geninfo spectrum_count_module; do
  python3 Hephaestos.py default mf 1 $WORK/src src_orig/$f.f90 > /dev/null || exit 1
done

cd $WORK
gfortran -O2 -cpp -DUSE_MPI=0 -J . -c src/compilation.f90
gfortran -O2 -cpp -DUSE_MPI=0 -J . -c src/geninfo.f90
gfortran -O2 -g -fcheck=all -fbacktrace -J . -c src/spectrum_count_module.f90
gfortran -O2 -J . -o chebyshev_count_unit_tests \
         ../../chebyshev_count_unit_tests.f90 compilation.o geninfo.o spectrum_count_module.o -llapack -lblas

./chebyshev_count_unit_tests $MESH

cd ..
rm -rf $WORK
