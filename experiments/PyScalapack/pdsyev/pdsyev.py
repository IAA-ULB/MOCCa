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

m  = n  = 8
mb = nb = 2

context_order = FORTRAN_ORDER

with scalapack(context_order, 2, 2) as context:
    
    # Create the local matrix 
    A_sub = context.array(m, n, mb, nb, dtype=float)
    A_sub.data[...] = np.random.randn(*A_sub.data.shape) # we fill the entire submatrix because that is easier. 
    Z_sub = context.array(m, n, mb, nb, dtype=float)
    # 
    print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")
    # print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")

    # work array
    # size of workarray
    lel1 = nb**2*(n/(nb*context.nprow.value)+n/(nb*context.npcol.value)+3)
    lel2 = nb*((n-1)/(nb*context.nprow.value*context.npcol.value)+1)
    lwork = int(n*5 + max(2*n,lel1) + n*lel2 + 1)
    print(f"{lwork=}")
    work = context.array(lwork,1,mb,nb,dtype=float)

    w = context.array(n,1,mb,nb,dtype=float)

    info = 0
    scalapack.pdsyev(
        *A_sub.data.shape,
        b'N', # N or V
        b'U', # U or L
        n,
        A_sub, 1, 1, *A_sub.scalapack_params(),
        w,
        Z_sub, 1, 1, *Z_sub.scalapack_params(),
        work, lwork,
        info
    )
    