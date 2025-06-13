#! /bin/bash

# First try to build the pasta code

ml calcua/2024a
ml ScaLAPACK
ml Python
ml SciPy-bundle
ml make

ml 

mkdir -p src
mkdir -p mod
mkdir -p obj

make INCLUDE=gnu-parallel-4Vaughan CALCTYPE=PASTA CONFIG=BXL
