#!/bin/bash

make && mpirun -np 4 ./exe $@