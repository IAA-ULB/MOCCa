from timeit import default_timer as timer
from time import sleep
import sys

import numpy as np
import math

from mpi4py import MPI

# we ar trying to pass the raw mpi communicator to PyScalapack
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
nranks = comm.Get_size()

if rank==0:
    # print(f"({rank}/{nranks}: {comm.__dir__()})")
    print(f"{comm}")
    print(f"MPI.COMM_WORLD.py2f() = {MPI.COMM_WORLD.py2f()}")

import PyScalapack
scalapack = PyScalapack(
    "/apps/antwerpen/zen2/rocky8/impi/2021.13.0-intel-compilers-2024.2.0/mpi/2021.13/lib/libmpi.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_core.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_sequential.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_intel_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_blacs_intelmpi_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_scalapack_lp64.so",
    )

with scalapack(b'C', 1, 1) as context0:
    # This must be executed on all ranks! 
    # Verify that nranks is a perfect square.

    if context0:
        print(f"context0 ({context0.rank.value}/{context0.size.value}) -> ({context0.myrow.value},{context0.mycol.value})")
    
with scalapack(b'C', 2, 2) as context:
    # This must be executed on all ranks! 
    # Verify that nranks is a perfect square.

    print(f"context* ({context.rank.value}/{context.size.value}) -> ({context.myrow.value},{context.mycol.value})")

print("end of Python script")
sleep(5)
