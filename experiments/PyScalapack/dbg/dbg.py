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

m = 5
n = 7
mb = nb = 1

context_order = FORTRAN_ORDER

with (
        scalapack(context_order, 1, 1) as context,
        scalapack(context_order, 1, 1) as context0,
):
    # Create the full matrix
    A_full = context0.array(m, n, 1, 1, dtype=float)
    # print(f"{context0.array.__doc__=}")
    if context0:
        print(f"context0: rank {context0.rank.value}/{context0.size.value}")
        A_full.data[...] = np.zeros(shape=(m,n),dtype=float)
        print(A_full.data.shape)
        print(A_full.data.flags)
        for i in range(m):
            for j in range(n):
                # print(f"[{i},{j}] -> {10.0*(i+1) + (j+1)}")
                A_full.data[i,j] = 10.0*(i+1) + (j+1)

        if context_order == FORTRAN_ORDER and not np.isfortran(A_full.data):
            raise RuntimeError('context_order (column-major) is not respected!')
        elif context_order == CPP_ORDER and np.isfortran(A_full.data):
            raise RuntimeError('context_order (row-major) is not respected!')

        print(A_full.data)

    # Create the local matrix 
    A_sub = context.array(m, n, mb, nb, dtype=float)
    print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")
    scalapack.pdgemr2d(
        *(m,n),
        *A_full.scalapack_params(),
        *A_sub.scalapack_params(),
        context.ictxt,
    )
    print("before")
    a = np.array([1.,0,2.,0,3.,0],dtype=float)
    print(f"{a=}")
    sum = scalapack.pcmax1(3, a, 1);
    print("after")
    print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")