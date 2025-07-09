# the equivalent of hello_from_BLACS.F90
import sys
if len(sys.argv) == 1:
    ___MPI__ = False
else:
    ___MPI__ = bool(sys.argv[1])
print(f"Using mpi4py: {___MPI__}")    

import numpy as np

if ___MPI__:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()
    if rank==0:
        print(f"mpi4py: ({rank}/{size}: {comm.handle=})")
else:
    print(f"mpi4py: is not used.")

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
    # Verify that size is a perfect square.

    if context0:
        print(f"context0 ({context0.rank.value}/{context0.size.value}) -> ({context0.myrow.value},{context0.mycol.value})")
    
with scalapack(b'C', 2, 3) as context:
    # This must be executed on all ranks! 
    # Verify that size is a perfect square.

    print(f"context* ({context.rank.value}/{context.size.value}) -> ({context.myrow.value},{context.mycol.value})")

print("end of Python script")
