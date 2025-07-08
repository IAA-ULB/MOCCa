#!/bin/bash

rm -rf _build

# python -m numpy.f2py --f90exec="/opt/cray/pe/craype/2.7.31.11/bin/ftn" --f77exec="/opt/cray/pe/craype/2.7.31.11/bin/ftn" -c hello.F90 -m hello --build-dir _build
python -m numpy.f2py  -c hello.F90 -m hello --build-dir _build


# python test_hello.py