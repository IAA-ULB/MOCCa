#!/bin/bash
# Forward commandline arguments to ./exe
make && mpirun -np 4 ./exe $@