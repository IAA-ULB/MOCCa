# Simple python script that takes the same command line arguments as test_real_example_nonrandom.f90
# to check its solution.
# this script is sequential

import numpy as np
import sys
import pyscalapack

scalapack = pyscalapack(
    "/apps/antwerpen/zen2/rocky8/impi/2021.13.0-intel-compilers-2024.2.0/mpi/2021.13/lib/libmpi.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_core.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_sequential.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_intel_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_blacs_intelmpi_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_scalapack_lp64.so",
)  

__symmetric__ = True

def aij(i,j):
    if not i>0:
        raise ValueError
    if not j>0:
        raise ValueError
    if __symmetric__:
        return 1.5 if i==j else 1/np.sqrt(np.abs(i-j))
    else:
        # print(f"--aij({i},{j})")
        return 10*i+j

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
    
    with ( scalapack(b'C', 2, 2) as context
         , scalapack(b'C', 1, 1) as context0
         ):
            
        nranks = context0.size.value
        nprow = int(np.sqrt(nranks))
        if nranks == nprow*nprow:
            npcol = nprow            
        else:
            raise ValueError(f"bad number of {ranks=}")
        
        # Global matrix (numpy array, on all ranks)
        a0n = np.zeros((na,na), dtype=float, order='F')
        # Global matrix (pyscalapack array only on context0)
        a0s = context0.array(na, na, nblk, nblk, dtype=float)
        if context0:
            print(f"a0s {a0s.local_m}x{a0s.local_n}")
            if not a0s.local_m>0:
                raise ValueError
        # initialize a0n on all ranks and a0s on context0 only
        for i in range(na):     # 0-based index
            for j in range(na): # 0-based index
                a0n[i,j] = aij(i+1,j+1) # on the left hand we need 0-based indices to index the numpy array
                                        # on the left hand we need 1-based indices for the underlying Fortran representation
                if context0:
                    a0s.data[i,j] = a0n[i,j] # indexing numpy arrays requires 0-based indices
            if context0:
                print(f"a0n[i,:] {a0n[i,:]}")
                print(f"a0s[i,:] {a0s.data[i,:]}")
            
        if context0:
            print(f"Solving sequentially on rank={context0.rank.value}.")
            for i in range(na):
                print(a0n[i,:])
            eigenvalues0, _ = np.linalg.eigh(a0n)

            print("*** np.linalg.eigh(a) solution ***")
            print("Eigenvalues:")
            for i in range(na):
                print(f"{i} {eigenvalues0[i]}")

        if a0s.context:
            print(f"[{context.rank.value}/{context.size.value}] a0s\n{a0s.data}")

        # distribute a0s into a0d
        a0d = context.array(na, na, nblk, nblk, dtype=float)
        # scalapack.pdgemr2d(
        #     *(na,na),
        #     *a0s.scalapack_params(),
        #     *a0d.scalapack_params(),
        #     context.ictxt,
        # )
        a0s.pdgemr2d(a0d)
        print(f"[{context.rank.value}/{context.size.value}] a0d\n{a0d.data}")

        # local matrix
        a = context.array(na, na, nblk, nblk, dtype=float)
        # print(f"[{context.rank.value}/{context.size.value}] a = {a.local_m}x{a.local_n}")

        for il in range(a.local_m):     # il is 0-based
            i = a.row_indxl2g(il+1)     # il+1 and i are 1-based
            for jl in range(a.local_n): # jl = 0-based
                j = a.col_indxl2g(jl+1) # jl+1 and j are 1-based
                # print(f"--a({i},{j})")
                a.data[il,jl] = aij(i,j)

        print(f"[{context.rank.value}/{context.size.value}] a\n{a.data}")

        for il in range(a.local_m):
            for jl in range(a.local_n):
                print(f"a[{il},{jl}] : {a.data[il,jl]} =? {a0d.data[il,jl]} : {a.data[il,jl] == a0d.data[il,jl]}")

        #--------------------------------------------------------------------------------------------------------------
        # Compute the parallel solution

        preallocate_eigenvalues_eigenvectors = False
        if preallocate_eigenvalues_eigenvectors:
            eigenvalues = np.zeros(na,dtype=float)
            # distributed array for the eigenvectors
            eigenvectors = context.array(na, na, nblk, nblk, dtype=float)
            a.pdsyev(eigenvalues=eigenvalues, eigenvectors=eigenvectors)        
        else:
            eigenvalues,eigenvectors,info = a.pdsyev()        

        if context0:
            print("*** pyscalapack.pdsyev solution ***")
            print("Eigenvalues:")
            for i in range(na):
                print(f"{i} {eigenvalues[i]}")
            norm_diff = np.linalg.norm(eigenvalues-eigenvalues0)
            print(f"{norm_diff=}")
            if norm_diff>1e-12:
                print("Parallel and sequential solutions might differ.")
