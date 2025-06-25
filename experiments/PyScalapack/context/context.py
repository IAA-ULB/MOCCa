import numpy as np
import PyScalapack
import sys
nprow = int(sys.argv[1])
npcol = int(sys.argv[2])
import os
size = os.environ['MPI_LOCALNRANKS']
if size == 1:
    nprow = npcol = 1
elif size == 4:
    nprow = npcol = 4

# print(f"{size=}")
# for k,v in os.environ.items():
#     if 'MPI' in k:
#         print(f"{k}: {v}")

scalapack = PyScalapack(
    "/apps/antwerpen/zen2/rocky8/impi/2021.13.0-intel-compilers-2024.2.0/mpi/2021.13/lib/libmpi.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_core.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_sequential.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_intel_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_blacs_intelmpi_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_scalapack_lp64.so",
    ) 

# Fortran uses Column-major order 
# Python and C++ use Row-major order
FORTRAN_ORDER = b'C' # numpy uses order='F'
CPP_ORDER     = b'R' # numpy uses order='C'

context_order = FORTRAN_ORDER

# matrix size: n x n 
n = 64

bf = 8
with scalapack(context_order, nprow, npcol) as context:
    print(f"({context.rank.value},{context.size.value})")
