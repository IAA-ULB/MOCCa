#!/usr/bin/env bash

#rm -rf _build_nil8_
python -m numpy.f2py -c nil8.f90 -m nil8_f90 -llapack #--build-dir _build_nil8_
if [ $? -eq 0 ]; then
    echo "SUCCESS"
    cp nil8_f90.cpython-*.so ../mocca/mean_field
else
    echo FAIL
fi