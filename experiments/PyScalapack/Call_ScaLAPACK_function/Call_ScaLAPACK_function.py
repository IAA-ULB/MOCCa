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

L1 = 128
L2 = 512
with (
        scalapack(b'C', 2, 2) as context,
        scalapack(b'C', 1, 1) as context0,
):
    array0 = context0.array(L1, L2, 1, 1, dtype=np.float64)
    if context0:
        array0.data[...] = np.random.randn(*array0.data.shape)

    array = context.array(L1, L2, 1, 1, dtype=np.float64)
    scalapack.pgemr2d["D"](
        *(L1, L2),
        *array0.scalapack_params(),
        *array.scalapack_params(),
        context.ictxt,
    )

    result = context.array(L1, L1, 1, 1, dtype=np.float64)
    scalapack.pdgemm(
        b'N',
        b'T',
        *(L1, L1, L2),
        scalapack.d_one,
        *array.scalapack_params(),
        *array.scalapack_params(),
        scalapack.d_zero,
        *result.scalapack_params(),
    )

    result0 = context0.array(L1, L1, 1, 1, dtype=np.float64)
    scalapack.pgemr2d["D"](
        *(L1, L1),
        *result.scalapack_params(),
        *result0.scalapack_params(),
        context.ictxt,
    )

    if context0:
        error = result0.data - array0.data @ array0.data.T
        print(np.linalg.norm(error))