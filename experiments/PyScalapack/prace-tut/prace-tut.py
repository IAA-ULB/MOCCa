# The PyScalapack version of `tantalus_full/experiments/ScaLAPACK/prace-tut/src.f90`
# `
# Run this script as: 
#   > mpirun -n 4 python prace-tut.py 
# Output:
#   Matrix dimension is (23, 47)
#   Matrix local dimension at process (0, 0) is (13, 25)
#   Matrix local dimension at process (1, 0) is (10, 25)
#   Matrix local dimension at process (0, 1) is (13, 22)
#   Matrix local dimension at process (1, 1) is (10, 22)
# The output matches the documentation

import PyScalapack
# Note: the order of the libraries below is important.
scalapack = PyScalapack(
    "/apps/antwerpen/zen2/rocky8/impi/2021.13.0-intel-compilers-2024.2.0/mpi/2021.13/lib/libmpi.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_core.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_sequential.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_intel_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_blacs_intelmpi_lp64.so",
    "/apps/antwerpen/zen2/rocky8/imkl/2024.2.0/mkl/2024.2/lib/intel64/libmkl_scalapack_lp64.so",
    ) 
import numpy as np

"""
Matrix dimension is (5, 5)
rank=0/4 : processor grid coordinate=(0,0), local matrix dimension=3x3
(i=0,j=0) -> (il=0,0) on grid coordinate (0,0)
(i=0,j=1) -> (il=0,1) on grid coordinate (0,0)
(i=0,j=2) -> (il=0,2) on grid coordinate (0,0)
                                                (i=0,j=3) -> (il=0,1) on grid coordinate (0,1)
                                                (i=0,j=4) -> (il=0,2) on grid coordinate (0,1)
(i=1,j=0) -> (il=1,0) on grid coordinate (0,0)
(i=1,j=1) -> (il=1,1) on grid coordinate (0,0)
(i=1,j=2) -> (il=1,2) on grid coordinate (0,0)
                                                (i=1,j=3) -> (il=1,1) on grid coordinate (0,1)
                                                (i=1,j=4) -> (il=1,2) on grid coordinate (0,1)
(i=2,j=0) -> (il=2,0) on grid coordinate (0,0)
(i=2,j=1) -> (il=2,1) on grid coordinate (0,0)
(i=2,j=2) -> (il=2,2) on grid coordinate (0,0)
                                                (i=2,j=3) -> (il=2,1) on grid coordinate (0,1)
                                                (i=2,j=4) -> (il=2,2) on grid coordinate (0,1)
                                                (i=3,j=0) -> (il=1,0) on grid coordinate (1,0)
                                                (i=3,j=1) -> (il=1,1) on grid coordinate (1,0)
                                                (i=3,j=2) -> (il=1,2) on grid coordinate (1,0)
                                                (i=3,j=3) -> (il=1,1) on grid coordinate (1,1)
                                                (i=3,j=4) -> (il=1,2) on grid coordinate (1,1)
                                                (i=4,j=0) -> (il=2,0) on grid coordinate (1,0)
                                                (i=4,j=1) -> (il=2,1) on grid coordinate (1,0)
                                                (i=4,j=2) -> (il=2,2) on grid coordinate (1,0)
                                                (i=4,j=3) -> (il=2,1) on grid coordinate (1,1)
                                                (i=4,j=4) -> (il=2,2) on grid coordinate (1,1)

(i=0,j=0) -> (il=0,0) on grid coordinate (0,0)
(i=0,j=1) -> (il=0,1) on grid coordinate (0,0)
(i=0,j=2) -> (il=0,2) on grid coordinate (0,0)
                                                (i=0,j=3) -> (il=0,5) on grid coordinate (0,1)
                                                (i=0,j=4) -> (il=0,6) on grid coordinate (0,1)
(i=1,j=0) -> (il=1,0) on grid coordinate (0,0)
(i=1,j=1) -> (il=1,1) on grid coordinate (0,0)
(i=1,j=2) -> (il=1,2) on grid coordinate (0,0)
                                                (i=1,j=3) -> (il=1,5) on grid coordinate (0,1)
                                                (i=1,j=4) -> (il=1,6) on grid coordinate (0,1)
(i=2,j=0) -> (il=2,0) on grid coordinate (0,0)
(i=2,j=1) -> (il=2,1) on grid coordinate (0,0)
(i=2,j=2) -> (il=2,2) on grid coordinate (0,0)
(i=2,j=3) -> (il=2,5) on grid coordinate (0,1)
(i=2,j=4) -> (il=2,6) on grid coordinate (0,1)
(i=3,j=0) -> (il=5,0) on grid coordinate (1,0)
(i=3,j=1) -> (il=5,1) on grid coordinate (1,0)
(i=3,j=2) -> (il=5,2) on grid coordinate (1,0)
(i=3,j=3) -> (il=5,5) on grid coordinate (1,1)
(i=3,j=4) -> (il=5,6) on grid coordinate (1,1)
(i=4,j=0) -> (il=6,0) on grid coordinate (1,0)
(i=4,j=1) -> (il=6,1) on grid coordinate (1,0)
(i=4,j=2) -> (il=6,2) on grid coordinate (1,0)
(i=4,j=3) -> (il=6,5) on grid coordinate (1,1)
(i=4,j=4) -> (il=6,6) on grid coordinate (1,1)
                                                
"""

with scalapack(b'R', 2, 2) as context:
    Al = context.array(m=5, n=5, mb=2, nb=2, dtype=np.float64)

    if context.rank.value == 0:
        print(f"Matrix dimension is ({Al.m}, {Al.n})")

    # print '("[",i0,"/",i0,"] myrow="i0" mycol="i0" nprow="i0" npcol="i0"")', myid, nproc, myrow, mycol, nprow, npcol
        print(f"rank={context.rank.value}/{context.size.value} : processor grid coordinate=({context.myrow.value},{context.mycol.value}), local matrix dimension={Al.local_m}x{Al.local_n}")

    # if context.rank.value == 0:
    #     for i in range(5):
    #         il  = scalapack.indxg2l(i,2,0,0,2)
    #         ipr = scalapack.indxg2p(i,2,0,0,2)
    #         for j in range(5):
    #             jl  = scalapack.indxl2g(j,2,0,0,2)
    #             ipc = scalapack.indxg2p(j,2,0,0,2)
    #             print(f"({i=},{j=}) -> ({il=},{jl}) on grid coordinate ({ipr},{ipc})")
    for i in range(5):
        for j in range(5):
            print(f"rank={context.rank.value}/{context.size.value} : {i=} {j=}")
            scalapack.pdelset( Al, i, j, 10*(i+1) + j+1 )    

    # for rank in range(4):
    #     if context.rank.value == rank:
    #         print(f"\nMatrix local dimension at process " +  
    #               f"({context.myrow.value}, {context.mycol.value})" +  
    #               f" is ({Al.local_m}, {Al.local_n})"
    #               )
    #         for il in range(Al.local_m):
    #             i = scalapack.indxl2g(il,2,context.mycol.value,0,4)
    #             for jl in range(Al.local_n):
    #                j = scalapack.indxl2g(jl,2,context.myrow.value,0,4)
    #                print(f"({il=}'{jl=}) -> ({i=},{j=})")

            # for i in range(5):
            #     il  = scalapack.indxg2l(i,2,0,0,2)
            #     ipr = scalapack.indxg2p(i,2,0,0,2)
            #     for j in range(5):
            #         jl  = scalapack.indxg2l(j,2,0,0,2)
            #         ipc = scalapack.indxg2p(j,2,0,0,2)
            #         if ipr==context.myrow.value and ipc==context.mycol.value:
            #             print(f"r{context.rank.value}: ({ipr=},{ipc=}) ({il=},{jl=}) ({i=},{j=})")

    #         Aij = 10.0*i + j
    #         if ipr==context.myrow.value and ipc==context.mycol.value:
    #             Al.data[il,jl] = Aij
