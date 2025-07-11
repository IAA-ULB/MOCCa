#!/usr/bin/env python 

# Simple python script that takes the same command line arguments as test_real_example_nonrandom.f90
# to check its solution.
# this script is sequential

import numpy as np
import sys

if __name__ == "__main__":

    # Default parameters
    na = 8
    nblk = 2
    
    # Command line arguments
    nargs = len(sys.argv)
    if nargs > 1:
        na = int(sys.argv[1])
    if nargs > 2:
        nblk = int(sys.argv[2])
    
    # Fix nev
    nev = na
    
    a = np.zeros((na,na), dtype=float, order='F')

    for i in range(na):
        for j in range(na):
            if i==j:
                a[i,i] = 1.5
            else:
                a[i,j] = 1/np.sqrt(np.abs(i-j))
        print(a[i,:])
    eigenvalues, _ = np.linalg.eigh(a)

    print("Eigenvalues:")
    for i in range(na):
        print(f"{i} {eigenvalues[i]}")