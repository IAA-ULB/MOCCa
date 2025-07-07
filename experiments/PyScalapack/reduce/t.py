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

m = 4
n = 2
mb = 1
nb = 1

context_order = FORTRAN_ORDER

with scalapack(context_order, 2, 2) as context:
    # print("hello")
    nranks = context.size.value
    rank = context.rank.value
    # if context0:xx
    #     R0 = context0.array(2*nranks, 2, 1, 1, dtype=float)
    #     print(f"[{context.rank.value}/{context.size.value}] R0=\n{R0.data}")
    # Ra = context0.array(2*nranks, 2, 1, 1, dtype=float)
    # print(f"[{context.rank.value}/{context.size.value}] Ra=\n{Ra.data}")
    a = np.zeros(2,dtype=float)
    a[0] = rank+1
    a[1] = (rank+1)*10
    
    print(f"[{context.rank.value}/{context.size.value}] 0 {a=}")

    scalapack.dgsum2d( context.ictxt, b'A', b'1', 1, 2, a, 1, 0, 0 )
    
    print(f"[{context.rank.value}/{context.size.value}] 1 {a=}")

