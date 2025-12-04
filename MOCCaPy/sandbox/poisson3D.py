r"""
This sandbox module is meant to solve the 3D Poisson equation arising in MOCCaPy.
It is built after `shenfun/demos/poisson3D.py`

We use a domain as for the LagrangeMesh
```
mesh = LagrangeMesh(M=30, d=.8, reduced=True)
```
    $$
    W = 15*0.8 = 12
    $$
    $$
    \nabla^2 u(x, y, z) = f(x, y, z), \quad (x, y, z) \in [0, W] \times [0,W] \times [0, W]
    $$

This module is used to test all sorts of boundary conditions
for bases using orthogonal polynomials, like::

    - Chebyshev first and second kind
    - Legendre
    - Jacobi

"""
import os
import sympy as sp
import numpy as np
from shenfun import inner, div, grad, TestFunction, TrialFunction, \
    Array, Function, FunctionSpace, TensorProductSpace, la, \
    dx, comm, cleanup
from mpi4py_fft.pencil import Subcomm
from mocca.mesh import LagrangeMesh

# Use sympy to compute a rhs, given an analytical solution
x, y, z = sp.symbols("x,y,z", real=True)
sigma = 2.0
ue = sp.exp(-0.5*(x**2 + y**2 + z**2)/sigma**2) # a gaussian curve
fe = ue.diff(x, 2) + ue.diff(y, 2) + ue.diff(z, 2)

d = 0.8
M = 30
mesh = LagrangeMesh(M=30, d=d, reduced=True)
w = 0.8*(M//2)

domain = ((0, w),
          (0, w),
          (0, w))
bc = {
    0: { 'left' : {'N': 0.}
       , 'right': {'D': 0.}
       },
    1: { 'left' : {'N': 0.}
       , 'right': {'D': 0.}
       },
    2: { 'left' : {'N': 0.}
       , 'right': {'D': 0.}
       },
}

def main(N, family, bc):

    SD1 = FunctionSpace(N, family=family, bc=bc[0], domain=domain[0], alpha=1, beta=1) # alpha and beta are neglected for all but Jacobi
    SD2 = FunctionSpace(N, family=family, bc=bc[1], domain=domain[1], alpha=1, beta=1) # alpha and beta are neglected for all but Jacobi
    SD3 = FunctionSpace(N, family=family, bc=bc[2], domain=domain[2], alpha=1, beta=1) # alpha and beta are neglected for all but Jacobi

    # Try the uncommon approach of squeezing SD between the two Fourier spaces
    T = TensorProductSpace(comm, (SD1,SD2,SD3))
    B = T.get_testspace(kind='G')
    X = T.mesh()

    u = TrialFunction(T)
    v = TestFunction(B)

    # Get f on quad points
    fj = Array(B, buffer=fe)

    # Compute right hand side of Poisson equation
    f_hat = Function(B)
    f_hat = inner(v, fj, output_array=f_hat)

    # Get left hand side of Poisson equation
    matrices = inner(v, div(grad(u)))

    # Create linear algebra solver
    Solver = la.Solver3D
    H = Solver(matrices)

    # Solve and transform to real space
    u_hat = Function(T)
    u_hat = H(f_hat, u_hat)
    uq = u_hat.backward()

    # Compare with analytical solution
    uj = Array(T, buffer=ue)
    error = np.sqrt(inner(1, (uj-uq)**2))
    if comm.Get_rank() == 0:
        print(f'poisson3D {family:s} L2 error = {error:2.6e}')
    if 'pytest 'in os.environ:
        assert error < 1e-6
    u_hat_000 = u_hat(0.4*np.ones((3,1), dtype=np.float64))
    print(u_hat_000)
    gr = np.transpose(mesh.grid)
    # Verify that transpose does not make a copy, the line below should print True.
    print(f"{np.shares_memory(gr, mesh.grid)=}")
    u_gr = np.empty(mesh.linear_size, dtype=np.float64)
    u_hat(gr, output_array=u_gr)
    # we expect this to print ~0
    print(u_hat_000-u_gr[0])
    cleanup((T, B))
    
if __name__ == '__main__':
    for family in ('legendre', 'chebyshev', 'chebyshevu', 'jacobi'):
        main(24, family, bc)
