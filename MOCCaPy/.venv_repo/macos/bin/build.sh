#!/bin/bash
# Build Tantalus.BXL.exe on my mac
# gfortran installed with `brew install gcc`
# hdf5 installed with `brew install hdf5`
mkdir -p src

cp make_include/make.inc.macos-gnu-debug-hdf5 make.inc
make mf CONFIG=BXL CALCTYPE=NUCLEI
mv ./exec/Tantalus.BXL.exe ./MOCCaPy/.venv/bin/mocca.dbg

cp make_include/make.inc.macos-gnu-hdf5 make.inc
make mf CONFIG=BXL CALCTYPE=NUCLEI
mv ./exec/Tantalus.BXL.exe ./MOCCaPy/.venv/bin/mocca.exe
