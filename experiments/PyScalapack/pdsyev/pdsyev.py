from timeit import default_timer as timer
import sys

import numpy as np
import math
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

n  = int(sys.argv[1])   # The dimension of the square matrix (n -> nxn matrix)
bf = int(sys.argv[2])   # The blocking factor
TEST = False            # Allow verifying that the sequential and the distributed system are the same and give the same result.
if len(sys.argv) > 3:
    TEST = sys.argv[3] == 'TEST'
# For testing the correctness of the script, TEST sets up the full matrix on rank 0 and solves the problem sequentially.

with scalapack(context_order, 1, 1) as context0:
    # This must be executed on all ranks! 
    # Verify that nranks is a perfect square.
    nranks = context0.size.value
    nprow = int(math.sqrt(nranks))
    if nranks == nprow*nprow:
        npcol = nprow
        
    else:
        raise ValueError(f"bad number of {ranks=}")
    
    if context0:
        print(f"\nInitializing {nprow}x{npcol}={nranks} processor grid for {n}x{n} system with {bf=}.", file=sys.stderr)

    with( scalapack(context_order, nprow, npcol) as context  # distributed system on all ranks
        , scalapack(context_order,     1,     1) as context0 # full system on rank 0
        ): 
        # print(f"({context.rank.value},{context.size.value})")

        if TEST:
            # Create and solve the full matrix on rank 0 for verifying the correctness of teh script.
            # Create the full matriX
            A_ful = context0.array(n, n, 1, 1, dtype=float)
            
            # print(f"{context0.array.__doc__=}")
            if context0:
                # print(f"context0: rank {context0.rank.value}/{context0.size.value} initializing full matrix", file=sys.stderr)
                # fill the matrix symmetrically
                for i in range(n):
                    for j in range(n):
                        # print(f"[{i},{j}] -> {10.0*(i+1) + (j+1)}")
                        if i == j:
                            A_ful.data[i,j] = 1.5
                        else:
                            A_ful.data[i,j] = np.sqrt(1.0/np.abs(i-j))
                            
                # print(A_ful.data)

                # Compute sequential solution
                try:
                    eigenvalues0, _ = np.linalg.eigh(A_ful.data, UPLO='L')
                    # UPLO='L' is default 
                    # print(f"sequential solution: eigenvalues:\n{eigenvalues0}")
                except np.LinAlgError:
                    print("Sequential solution did not converge", file=sys.stderr)

            # Create a distributed matrix to distribute the full matrix to, just to verify that directly setting A_sub does the same thing
            A_sub_dgemr2d = context.array(n, n, bf, bf, dtype=float)
            # print(f"context: rank {context.rank.value}/{context.size.value} = ({context.myrow.value},{context.mycol.value}): A_sub\n{A_sub.data}")
            scalapack.pdgemr2d(
                *(n,n),
                *A_ful.scalapack_params(),
                *A_sub_dgemr2d.scalapack_params(),
                context.ictxt,
            )

        # Create the distributed matrix and initialize it directly.
        A_sub = context.array(n, n, bf, bf, dtype=float)

        # fortran code from https://info.gwdg.de/wiki/doku.php?id=wiki:hpc:scalapack
        # ! initialize in parallel the local parts of the distributed matrix
        #       mll = NUMROC( m, mb, myrow, 0, nprow )
        #       nll = NUMROC( n, nb, mycol, 0, npcol )
        #       do jl = 1 , nll
        #          j = INDXL2G(jl,nb,mycol,0,npcol)
        #          do il = 1 , mll
        #             i = INDXL2G(il,mb,myrow,0,nprow)
        #             al(il+mll*(jl-1))= matval(i,j)
        #          end do
        #       end do        

        mll = scalapack.numroc( n, bf, context.myrow.value, 0, context.nprow.value )
        nll = scalapack.numroc( n, bf, context.mycol.value, 0, context.npcol.value )        
        # print(f"([{context.rank.value}/{context.size.value}]=({context.myrow.value},{context.mycol.value}) : {mll=} {nll=}")
        # print(f"([{context.rank.value}/{context.size.value}]=({context.myrow.value},{context.mycol.value}) : {A_sub.data.shape=}")
        for jl in range(1,1+nll):
            # This is how it works: 
            # . the indices fed into indxl2g (il and jl) are Fortran indices, i.e. they start from 1, not 0
            # . the indices that are returned from indxl2g (i and j) are also Fortran indices, i.e. they start from 1, not 0
            # BUT (!) the distributed sub-matrices are numpy arrays, and they are indexed in the Python way, i.e. starting from 0.
            j = scalapack.indxl2g(jl,bf,context.mycol.value,0,npcol)
            for il in range(1,1+mll):
                i = scalapack.indxl2g(il,bf,context.myrow.value,0,nprow)
                if i == j:
                    A_sub.data[il-1,jl-1] = 1.5
                else:
                    A_sub.data[il-1,jl-1] = np.sqrt(1.0/np.abs(j-i))
                # print(f"([{context.rank.value}/{context.size.value}]({context.myrow.value},{context.mycol.value}) : L({il},{jl}) -> G({i},{j}) = {A_sub.data[il-1,jl-1]}")

                if TEST:
                    # verify the equality of A_sub and A_sub_dgemr2d
                    if not A_sub.data[il-1,jl-1] == A_sub_dgemr2d.data[il-1,jl-1]: # Python indexing (0-based), hence il-1 and jl-1.
                        raise RuntimeError(f"{i-1},{j-1}")
                
        # for il0 in range(mll): # il0 is 0-based
        #     print(f"([{context.rank.value}/{context.size.value}]=({context.myrow.value},{context.mycol.value}) A_sub[{il0},:] = {A_sub.data[il0,:]}")

        # distributed vector for the eigenvectors
        eigenvectors_sub = context.array(n, n, bf, bf, dtype=float)

        # Compute parallel solution 
        start = timer()
        eigenvalues = np.zeros(n,dtype=float)
        work = np.zeros(1,dtype=float,order='F')
        lwork = -1
        info = np.array([n*n]) # an ordinary Python variable cannot be used as an output argument.
        scalapack.pdsyev(
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
        scalapack.pdsyev(
            b'V', b'L', n,
            *A_sub.scalapack_params(),
            eigenvalues,
            *eigenvectors_sub.scalapack_params(),
            work, lwork,
            info
        )
        if info[0] != 0:
            print(f"[{context.rank.value}/{context.size.value}]({context.myrow.value},{context.mycol.value}): distributed solution did not converge {info=}.", file=sys.stderr)

        stop = timer()
        wt = stop - start
        # Compute the number of doubles allocated
        # The full matrix is not accounted for because it is only used if TEST==True
        nDP = mll*nll*2       # distributed matrix and eigenvectors
        nDP += n + lwork      # eigenvalues and work
        print(f"> {context.rank.value}, {context.size.value}, {context.nprow.value}, {context.npcol.value}, {n}, {bf}, {wt}, {nDP}")

        wt_nDP = np.array([wt,nDP])
        # Compute the element-wise sum of wt_nDP over all processes in context:
        scalapack.dgsum2d( context.ictxt, b'A', b'1', 1, 2, wt_nDP, 1, 0, 0 )
        if context0:
            print(f">>{context.rank.value}, {context.size.value}, {context.nprow.value}, {context.npcol.value}, {n}, {bf}, {wt_nDP[0]}, {int(wt_nDP[1])}")
        
        if TEST and context0:
            print(f"Test: {eigenvalues0=}")
            print(f"Test: {eigenvalues=}")
            diffnorm = np.linalg.norm(eigenvalues0-eigenvalues)
            print(f"Test: {n}x{n} {bf=}: difference norm = {diffnorm}", file=sys.stderr)
            if not diffnorm < 1e-12:
                raise RuntimeError(f"The sequential and the distributed solution differ.")