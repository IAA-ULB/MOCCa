#!/usr/bin/env bash

module=nil8
python -m numpy.f2py -c ${module}.f90 -m ${module}_f90 -llapack #--build-dir _build_nil8_
if [ $? -eq 0 ]; then
    echo "SUCCESS : ${module}_f90"
    cp ${module}_f90.cpython-*.so ../mocca/mean_field
else
    echo "FAILED : ${module}_f90"
fi
