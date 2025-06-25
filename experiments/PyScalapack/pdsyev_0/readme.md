Building on the experiments/PyScalapack/gwdg_pdgemr2d_ok case to create an 8x8 eigenvector problem that we solve in two ways:
- sequentially using `numpy.linalg.eigh`  ()
- in parallel on a 4x4 grid using `PyScalapack`
