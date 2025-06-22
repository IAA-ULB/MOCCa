#!/bin/bash

# . ../../../vaughan/ml.sh -v -p

command='mpirun -n 4 python prace-tut.py'
echo
echo ${command}
${command}

