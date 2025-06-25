import numpy as np
import PyScalapack

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

n  = 8
nb = 8

context_order = FORTRAN_ORDER

with scalapack(context_order, 1, 1) as context:
    print(f"({context.rank.value},{context.size.value})")    
        # Create the full matrix

    A_ful = context.array(n, n, n, n, dtype=float)
    # fill the matrix symmetrically
    for i in range(n):
        for j in range(n):
            # print(f"[{i},{j}] -> {10.0*(i+1) + (j+1)}")
            if i <= j:
                A_ful.data[i,j] = 10.0*(i+1) + (j+1)
            else:
                A_ful.data[i,j] = 10.0*(j+1) + (i+1)

    for i in range(n):
        print(A_ful.data[i,:])

    
    eigenvalues  = np.zeros(n,dtype=float)
    eigenvectors = context.array(n, n, n, n, dtype=float)
    work         = np.zeros(1,dtype=float)

    info = -1
    info = scalapack.pdsyev(
        b'V', b'L', n,
        *A_ful.scalapack_params(),
        eigenvalues,
        *eigenvectors.scalapack_params(),
        work, -1,
        info
    )
    lwork = int(work[0])
    work = np.zeros(lwork,dtype=float)
    info = scalapack.pdsyev(
        b'V', b'L', n,
        *A_ful.scalapack_params(),
        eigenvalues,
        *eigenvectors.scalapack_params(),
        work, lwork,
        info
    )
    print(f"{eigenvalues=}")
