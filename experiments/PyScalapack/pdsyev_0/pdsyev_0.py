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



context_order = FORTRAN_ORDER

m = n = 8

mb = nb = -1
with scalapack(context_order, 1, 1) as context0:
    if context0.size.value == 1:
        mb = nb = 1
    elif context0.size.value == 4:
        mb = nb = 2


with (
        scalapack(context_order, mb, nb) as context,
        scalapack(context_order, 1, 1) as context0,
):
    print(f"{nb=}")
    # Create the full matrix
    A_full = context0.array(m, n, 1, 1, dtype=float)
    # print(f"{context0.array.__doc__=}")
    if context0:
        print(f"context0: rank {context0.rank.value}/{context0.size.value}")
        A_full.data[...] = np.zeros(shape=(m,n),dtype=float)
        print(A_full.data.shape)
        print(A_full.data.flags)
        # fill the matrix symmetrically
        for i in range(m):
            for j in range(n):
                # print(f"[{i},{j}] -> {10.0*(i+1) + (j+1)}")
                if i <= j:
                    A_full.data[i,j] = 10.0*(i+1) + (j+1)
                else:
                    A_full.data[i,j] = 10.0*(j+1) + (i+1)

        if context_order == FORTRAN_ORDER and not np.isfortran(A_full.data):
            raise RuntimeError('context_order (column-major) is not respected!')
        elif context_order == CPP_ORDER and np.isfortran(A_full.data):
            raise RuntimeError('context_order (row-major) is not respected!')

        print(A_full.data)

        # Compute sequential solution
        eigenvalues, eigenvectors = np.linalg.eigh(A_full.data)
        print(f"sequential solution: eigenvalues:\n{eigenvalues}")
        print(f"sequential solution: eigenvectors:\n{eigenvectors}")

    # Create the local matrix 
    A_sub = context.array(m, n, mb, nb, dtype=float)
    print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")
    scalapack.pdgemr2d(
        *(m,n),
        *A_full.scalapack_params(),
        *A_sub.scalapack_params(),
        context.ictxt,
    )
    print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")

    eigenvectors_sub = context.array(m, n, mb, nb, dtype=float)
    # 
    print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")
    # print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")

    # Compute parallel solution 
    lel1 = nb**2*(n/(nb*context.nprow.value)+n/(nb*context.npcol.value)+3)
    lel2 = nb*((n-1)/(nb*context.nprow.value*context.npcol.value)+1)
    lwork = int(n*5 + max(2*n,lel1) + n*lel2 + 1)
    # print(f"{lwork=}")
    # c_double_p = scalapack.ctypes.POINTER(scalapack.ctypes.c_double)
    eigenvalues = np.zeros(n,dtype=float)
    work = np.zeros(1,dtype=float,order='F')
    lwork = -1
    info = -1
    info = scalapack.pdsyev(
        b'V', b'L', n,
        *A_sub.scalapack_params(),
        eigenvalues,
        *eigenvectors_sub.scalapack_params(),
        work, -1,
        info
    )
    lwork = int(work[0])
    print(f"{lwork=}")
    work = np.zeros(lwork,dtype=float,order='F')
    info = scalapack.pdsyev(
        b'V', b'L', n,
        *A_sub.scalapack_params(),
        eigenvalues,
        *eigenvectors_sub.scalapack_params(),
        work, lwork,
        info
    )
    print(f"{eigenvalues=}")
