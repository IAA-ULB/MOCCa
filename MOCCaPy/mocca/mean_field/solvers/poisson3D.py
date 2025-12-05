r"""
This module is meant to solve the 3D Poisson equation arising in MOCCaPy.
It is built after `shenfun/demos/poisson3D.py`
"""
import os
import sympy as sp
import numpy as np
from shenfun import inner, div, grad, TestFunction, TrialFunction, \
    Array, Function, FunctionSpace, TensorProductSpace, la, \
    dx, comm, cleanup

from mocca.util.timer import Timer


def poisson3D(rhs, N=None, family='Chebyshev', out=None, verbose=False):
    """Solver Poisson equation with non-periodic boundary conditions.
    Args:
        rhs: a MeshQuantity representing the rhs of the Poisson equation function.
            This information must be transferred to the quadrature points of
            TensorProductSpace.mesh().
        N: number of spectral basis functions
        family: family of spectral basis functions, `shenfun` supports: `Chebyshev`|`C`,
            `Chebyshevu`|`U`, `Legendre`|`L`, `Fourier`|`F`, `Laguerre`|`La,`, `Hermite`|`H`,
            `Jacobi`|`J`, `Ultraspherical`|`Q`
        out: output array, ndarray of length `rhs.mesh.linear_size`. If `out='rhs'` the rhs (`rhs.data`) is
            overwritten with the solution.
    """
    if N is None:
        N = [len(g1d) for g1d in rhs.mesh.g1D]
    elif isinstance(N, int):
        N = 3*[N]
    else:
        assert len(N) == 3

    # rhs_symmetry: a tuple (±1,±1,±1) implying that rhs is symmetric (+1) or skew-symmetric on the
    #     corresponding coordinate axis. Ignored if the axis is not reduced. This determines the BC
    #     on the left end of reduced axes. A symmetric rhs implies a Neumann BC (derivative = 0), a
    #     skew-symmetric rhs implies a Dirichlet BC (value = 0).
    mesh = rhs.mesh
    assert mesh.dim == 3
    # Domain and BC
    domain = []
    bc = []
    for id in range(rhs.mesh.dim):
        w = mesh.d[id] * (mesh.M[id] // 2)
        if mesh.reduced[id]:
            left, right = 0, w
            bc_left  = {'N': 0.} if (rhs.symmetry[id,0] == 1) else {'D': 0.}
            # TODO: use eq 23 from PhysRevC92
            bc_right = {'D': 0.}
        else:
            left, right = -w, w
            bc_left  = {'D': 0.}
            bc_right = {'D': 0.}
        domain.append((left, right))
        bc.append({
            'left': bc_left,
            'right': bc_right
        })

    SD1 = FunctionSpace(N[0], family=family, bc=bc[0], domain=domain[0], alpha=1, beta=1) # alpha and beta are neglected for all but Jacobi
    SD2 = FunctionSpace(N[1], family=family, bc=bc[1], domain=domain[1], alpha=1, beta=1) # alpha and beta are neglected for all but Jacobi
    SD3 = FunctionSpace(N[2], family=family, bc=bc[2], domain=domain[2], alpha=1, beta=1) # alpha and beta are neglected for all but Jacobi

    # Try the uncommon approach of squeezing SD between the two Fourier spaces
    T = TensorProductSpace(comm, (SD1, SD2, SD3))
    B = T.get_testspace(kind='G')

    u = TrialFunction(T)
    v = TestFunction(B)

    # Get f on quad points
    with Timer("Poisson3D.1 Transfer problem to shenfun (interpolate quadrature points)") as t:
        xyz_qp = T.mesh() # the mesh of quadrature points: a list of three ndarrays of shape [(N,1,1), (1,N,1), (1,1,N)]
        # TODO: interpolate T.mesh()
        rhs_qp = np.empty((*N,1), np.float64)
        x = xyz_qp[0]
        y = xyz_qp[1]
        z = xyz_qp[2]
        r = np.empty((N[2],3), np.float64)
        r[:,2] = z
        for i in range (N[0]):
            r[:,0] = x[i]
            for j in range (N[1]):
                r[:,1] = y[0,j,0]
                rhs_qp[i, j, :] = mesh.interpolate(rhs, r, out=rhs_qp[i, j, :, :])

    if verbose: print(t)

    fj = Array(B, buffer=rhs_qp)
    f_hat = Function(B)
    f_hat = inner(v, fj, output_array=f_hat)

    # Get left hand side of Poisson equation
    matrices = inner(v, div(grad(u)))

    # Create linear algebra solver
    Solver = la.Solver3D
    H = Solver(matrices)

    # Solve and transform to real space
    with Timer("Poisson3D.2 Compute shenfun solution") as t:
        u_hat = Function(T)
        u_hat = H(f_hat, u_hat)
    # uq = u_hat.backward()
    if verbose: print(t)

    gr = np.transpose(mesh.grid)
    # Verify that transpose does not make a copy
    assert np.shares_memory(gr, mesh.grid)
    # u_gr = np.empty(mesh.linear_size, dtype=np.float64)
    # u_hat(gr, output_array=u_gr)

    with Timer("Poisson3D.3 Transfer shenfun solution to LagrangeMesh" ) as t:
        if out is None:
            result = u_hat(gr)
        elif out == 'rhs':
            result = u_hat(gr, output_array=rhs.data[:,0])
        else:
            result = u_hat(gr, output_array=out)

    print(t)
    cleanup((T, B))
    return result


