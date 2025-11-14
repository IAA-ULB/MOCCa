#!/usr/bin/env bash

python -m numpy.f2py -c test_f2py.f90 -m fill_f90
# add -llapack