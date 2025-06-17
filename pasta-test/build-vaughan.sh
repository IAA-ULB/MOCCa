#! /bin/bash
# Use as 
# > . ./pasta-test/build-vaughan.sh
# (in directory tantalus_full.)

pwd=${PWD##*/}
if [[ "$pwd" != "tantalus_full" ]]; then
    >&2 echo "error: Script "build-vaughan" must be run in directory "TANTALUS_FULL" (top level)."
    >&2 echo "error:   current directory is '$pwd'."
    return 1
fi
. ./pasta-test/ml-vaughan.sh

# modules needed by Hephaestos
ml Python
ml SciPy-bundle
ml make

ml 

mkdir -p src
mkdir -p mod
mkdir -p obj

if [[ "$1" = "intel" ]]; then 
    include="intel-parallel--vaughan"
else
    include="gnu-parallel--vaughan"
fi

make INCLUDE=$include CALCTYPE=PASTA CONFIG=BXL
