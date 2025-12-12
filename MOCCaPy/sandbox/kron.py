import scipy
import numpy as np

def D1D(stencil, n):
    d0 = stencil[0] * np.ones(n, np.int8)
    d1 = stencil[1] * np.ones(n - 1, np.int8)
    d = scipy.sparse.diags_array((d1, d0, d1), offsets=(-1, 0, 1))
    return d


if __name__ == "__main__":
    Nx = 5
    Ny = 4
    Nz = 3
    stencil = (-2, 1)
    L_3D = scipy.sparse.kron(scipy.sparse.eye(Nz), scipy.sparse.kron(scipy.sparse.eye(Ny),     D1D(stencil, Nx))) \
         + scipy.sparse.kron(scipy.sparse.eye(Nz), scipy.sparse.kron(    D1D(stencil, Ny), scipy.sparse.eye(Nx))) \
         + scipy.sparse.kron(    D1D(stencil, Nz), scipy.sparse.kron(scipy.sparse.eye(Ny), scipy.sparse.eye(Nx)))
    print(L_3D)
    L_3D_full = L_3D.toarray()
    print("   0---+----+----+----+----+----+----+----+----+----+----+----+")
    for ixyz in range(Nx * Ny * Nz):
        s = f"{ixyz:2} "
        for jxyz in range(Nx * Ny * Nz):
            s += str(int(abs(L_3D_full[ixyz, jxyz])))
        print(s)