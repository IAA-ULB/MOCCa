#!/bin/bash
export USE_PYELPA=0
echo "USE_PYELPA=${USE_PYELPA}"
mpirun -n 4 python -m pytest -s -x tests

# export USE_PYELPA=1
# echo "USE_PYELPA=${USE_PYELPA}"
# mpirun -n 4 python -m pytest -s -x tests