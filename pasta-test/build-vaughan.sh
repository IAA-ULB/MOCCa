#! /bin/bash
# Use as 
# > . ./pasta-test/build-vaughan.sh
#   - uses the same command line arguments as ml-vaughan.sh
#   - must be run in directory tantalus_full.

# Verify present working directory
pwd=${PWD##*/}
if [[ "$pwd" != "tantalus_full" ]]; then
    >&2 echo "error: Script "build-vaughan" must be run in directory "tantalus_full" (top level)."
    >&2 echo "error:   current directory is '$pwd'."
    return 1
fi

# Load modules needed for building and running tantalus
. ./pasta-test/ml-vaughan.sh $@

# modules needed by Hephaestos
ml Python
ml SciPy-bundle
ml make
ml 
>&2 echo

# directories needed by hephaestos
mkdir -p src
mkdir -p mod
mkdir -p obj

case "$TOOLCHAIN" in
    intel)  include="intel-parallel-hdf5--vaughan";;
    # hdf5 must be added first!
    # gnu  )  include="gnu-parallel--vaughan";;
    *    )  >&2 echo "Toolchain unknown: ${TOOLCHAIN}."; return 1;;
esac

make INCLUDE=$include CALCTYPE=PASTA CONFIG=BXL

if [ $? = 0 ] ; then echo "Built: ${include}." >&2 ;fi