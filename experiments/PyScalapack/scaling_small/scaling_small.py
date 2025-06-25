from timeit import default_timer as timer
import sys

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

# matrix size: n x n 
nprow = npcol = -1

n  = int(sys.argv[1])
bf = int(sys.argv[2])

with scalapack(context_order, 1, 1) as context0:
    # This must be executed on all ranks! 
    nranks = context0.size.value
    if nranks in (1,4,9,16,25,36,49,64):
        nprow = npcol = int(np.sqrt(nranks))
    else:
        raise ValueError("bad number of ranks")
    if context0:
        print(f"Initializing {nprow}x{npcol} processor grid", file=sys.stderr)

    with( scalapack(context_order, nprow, npcol) as context  # distributed system on all ranks
        , scalapack(context_order,     1,     1) as context0 # full system on rank 0
        ): 
        # print(f"({context.rank.value},{context.size.value})")

        # Create the full matri
        A_ful = context0.array(n, n, 1, 1, dtype=float)
        # print(f"{context0.array.__doc__=}")
        if context0:

            print(f"context0: rank {context0.rank.value}/{context0.size.value} initializing full matrix", file=sys.stderr)
            # fill the matrix symmetrically
            for i in range(n):
                for j in range(n):
                    # print(f"[{i},{j}] -> {10.0*(i+1) + (j+1)}")
                    if i == j:
                        A_ful.data[i,j] = 1.5
                    # we fill only the lower triangle
                    elif i<j:
                        continue
                    #     A_ful.data[i,j] = np.sqrt(1.0/(j-i))
                    else:
                        A_ful.data[i,j] = np.sqrt(1.0/(i-j))

            # print(A_ful.data)

            # Compute sequential solution
            eigenvalues0, _ = np.linalg.eigh(A_ful.data, UPLO='L')
            # UPLO='L' is default 
            # print(f"sequential solution: eigenvalues:\n{eigenvalues0}")
            # print(f"sequential solution: eigenvectors:\n{eigenvectors}")

        # Create the distributed matrix 
        A_sub = context.array(n, n, bf, bf, dtype=float)
        # print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")
        scalapack.pdgemr2d(
            *(n,n),
            *A_ful.scalapack_params(),
            *A_sub.scalapack_params(),
            context.ictxt,
        )
        # print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")

        eigenvectors_sub = context.array(n, n, bf, bf, dtype=float)
        # 
        # print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")
        # print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")

        # Compute parallel solution 
        start = timer()
        eigenvalues = np.zeros(n,dtype=float)
        work = np.zeros(1,dtype=float,order='F')
        lwork = -1
        info = -1
        info = scalapack.pdsyev(
            b'V', b'L', n,
            *A_sub.scalapack_params(),
            eigenvalues,
            *eigenvectors_sub.scalapack_params(),
            work, lwork,
            info
        )
        lwork = int(work[0])
        # print(f"{lwork=}")
        work = np.zeros(lwork,dtype=float,order='F')
        if False:
            info = scalapack.pdsyev(
                b'V', b'L', n,
                *A_sub.scalapack_params(),
                eigenvalues,
                *eigenvectors_sub.scalapack_params(),
                work, lwork,
                info
            )
        stop = timer()
        dt = stop - start
        print(f"{context.rank.value}, {context.size.value}, {n}, {bf}, {dt}")
        # context.size.value, context.rank.value, n, bf, dt
        # with open("timings.txt", mode='a') as f:
        #     f.write(f"{context.size.value}, {context.rank.value}, {n}, {bf}, {dt}\n")
        if context0:
            n_doubles  = n*n                                # full matrix on rank 0
            n_doubles += n*n * 2 / context0.size.value      # distributed matrix and eigenvectors
            n_doubles += n + lwork                          # eigenvalues and work
            
            print(f"{n}x{n} {bf=}: difference norm = {np.linalg.norm(eigenvalues0-eigenvalues)} {dt}s\n"
                  f"# {n}x{n} : {n_doubles=}", file=sys.stderr)
