# This script copies the code in https://pypi.org/project/PyScalapack/, section 'Documents/Create an array`.
# `
# Run this script as: 
#   > mpirun -n 4 python Create_an_array.py 
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

with scalapack(b'C', 2, 2) as context:
    array = context.array(m=23, n=47, mb=5, nb=5, dtype=np.float64)

    if context.rank.value == 0:
        print(f"Matrix dimension is ({array.m}, {array.n})")
    print(f"Matrix local dimension at process " +  #
    f"({context.myrow.value}, {context.mycol.value})" +  #
    f" is ({array.local_m}, {array.local_n})")
