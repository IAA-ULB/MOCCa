#!/bin/bash

python -m numpy.f2py --f90exec=/opt/cray/pe/craype/2.7.31.11/bin/ftn --f77exec=/opt/cray/pe/craype/2.7.31.11/bin/ftn -c fib.F90 -m fib 

python -c "import fib; print(fib.fib(8))"