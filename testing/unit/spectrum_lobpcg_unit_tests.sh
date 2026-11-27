#!/usr/bin/env bash
#-------------------------------------------------------------------------------
# Unit tests of spectrum_lobpcg_module: solve for the lowest eigenpairs of a
# small mesh hamiltonian by fixed-block LOBPCG and compare against the exact
# lowest-n eigenvalues from a dense diagonalization. Tier:
#   - 8^3 smoke test (fast, suitable for CI)
#   - 16^3 reference tier: level spacing near the cutoff like a real run;
#     run locally before committing.
# Usage:
#   bash spectrum_lobpcg_unit_tests.sh [8|16]
#-------------------------------------------------------------------------------
MESH=${1:-8}
if [ "$MESH" != "8" ] && [ "$MESH" != "16" ]; then
  echo "Usage: bash spectrum_lobpcg_unit_tests.sh [8|16]"
  exit 1
fi

# Run from inside this directory, per the testing conventions.
cd "$(dirname "$0")"
ROOT=$(cd ../.. && pwd)

# Transient work directory, destroyed at teardown
WORK=$PWD/work/spectrum_lobpcg.$MESH
rm -rf $WORK ; mkdir -p $WORK/src

# Push the marker-bearing sources through Hephaestos; the LOBPCG module is
# marker-free but goes through the same pass for uniformity. Hephaestos
# resolves configs/ relative to the repository root.
cd $ROOT
for f in compilation geninfo spectrum_lobpcg_module; do
  python3 Hephaestos.py default mf 1 $WORK/src src_orig/$f.f90 > /dev/null || exit 1
done

cd $WORK
gfortran -O2 -cpp -DUSE_MPI=0 -J . -c src/compilation.f90
gfortran -O2 -cpp -DUSE_MPI=0 -J . -c src/geninfo.f90
gfortran -O2 -cpp -DUSE_MPI=0 -J . -c src/spectrum_lobpcg_module.f90
gfortran -O2 -J . -o spectrum_lobpcg_unit_tests \
         ../../spectrum_lobpcg_unit_tests.f90 compilation.o geninfo.o spectrum_lobpcg_module.o -llapack -lblas

./spectrum_lobpcg_unit_tests $MESH

cd ..
rm -rf $WORK
