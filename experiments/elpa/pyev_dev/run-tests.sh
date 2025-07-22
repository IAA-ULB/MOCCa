#!/bin/bash
# use pytest -x for stopping at the first failure

export USE_PYELPA=1
echo "USE_PYELPA=${USE_PYELPA}"
mpirun -n 4 python -m pytest -s tests

# this one fails often because not all Elpa functionality is implemented in pyev
export USE_PYELPA=0
echo "USE_PYELPA=${USE_PYELPA}"
mpirun -n 4 python -m pytest -s tests
