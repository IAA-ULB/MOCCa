"""
Solve the generalized Poisson equation

$$
( a \Delta + b ) u(r) = f(r)
$$

using the Finite Difference Method.

The current approach uses a 1D stencil in each dimension, and is NOT ready for stencils using diagonal points.
"""

import numpy as np
import scipy.sparse as sparse
from tabulate import tabulate
import pytest
from mocca.mesh import LagrangeMesh, MeshQuantity

# TODO: implement generalized poisson
# TODO: time the methods
# TODO: use out= for matvec ?

iterative_linear_solvers = {
    'cg'    : sparse.linalg.cg,
    'cgs'   : sparse.linalg.cgs,
    'gmres' : sparse.linalg.gmres,
    'lgmres': sparse.linalg.lgmres,
    'minres': sparse.linalg.minres,
}

class FDStencil:
    """Base class for system of linear equations resulting from a particular 1D Finite Difference stencil.
    """
    def __init__(self, mesh, boundary_width):
        """
        Args:
            mesh: LagrangeMesh instance on which the solution of the Poisson equation is sought.
            boundary_width: width of the boundary of the stencil expression. Provided by the
                derived class
        """
        self.mesh = mesh # this is the embedded mesh, not the embedding mesh
        self.boundary_width = boundary_width

        self.h_uniform = mesh.d[0]
        for idim in range(1, mesh.dim):
            if self.h_uniform != mesh.d[idim]:
                self.h_uniform = 0
                break

        # Create an embedding grid for mesh
        bmesh = LagrangeMesh(
            dim=mesh.dim,
            M=mesh.M, d=mesh.d, reduced=mesh.reduced,
            highest_derivative_order=0,
            boundary_width=self.boundary_width,
        )

        if mesh.dim > 1:
            # Collect boundary points and symmetry points
            def remove_corners(condition, idim):
                """Remove corners from symmetry points.
                Remark:
                    This should not be called for ND stencils with diagonal points.
                """
                for jdim in range(self.mesh.dim):
                    if jdim != idim:
                        L, R = bmesh.g1D[jdim][ self.boundary_width-1], \
                               bmesh.g1D[jdim][-self.boundary_width  ]
                        np.logical_and(condition, bmesh.grid[:,jdim] > L, out=condition)
                        np.logical_and(condition, bmesh.grid[:,jdim] < R, out=condition)
                return condition

            bp = []
            self.sp_l = [[None]*bmesh.dim for ib in range(self.boundary_width)]
            for idim in range(bmesh.dim):
                for ib in range(self.boundary_width):
                    l, r = bmesh.g1D[idim][0+ib], bmesh.g1D[idim][-1-ib]
                    if not bmesh.reduced[idim]:
                        # There are only boundary points on the left side if the axis is not reduced.
                        bp.append(np.nonzero(bmesh.grid[:, idim] == l)[0])
                    else:
                        # Otherwise they are symmetry points,
                        # - which we must keep separate for each boundary layer and each axis.
                        # - and we must avoid corners because they have no counterpart in the original mesh
                        # we remove the corners afterwards
                        self.sp_l[ib][idim] = np.nonzero(remove_corners(bmesh.grid[:, idim] == l, idim))[0]

                    bp.append(np.nonzero(bmesh.grid[:, idim] == r)[0])

            self.bp_l = np.sort(np.unique(np.concat(bp))) # ndarray with linear indices of the boundary points
            self.bp_xyz = bmesh.grid[self.bp_l]           # ndarray with coordinates    of the boundary points

            # determine the linear indices of the points in u corresponding to the symmetry points.
            self.sp_lu = [[None]*bmesh.dim for ib in range(self.boundary_width)]
            for idim in range(bmesh.dim):
                if bmesh.reduced[idim]:
                    for ib in range(self.boundary_width):
                        iu = self.boundary_width - ib -1
                        l = self.mesh.g1D[idim][iu]
                        self.sp_lu[ib][idim] = np.nonzero(mesh.grid[:, idim] == l)[0]

        else:
            bp = []
            self.sp_l = [[None] * bmesh.dim for ib in range(self.boundary_width)]
            for ib in range(self.boundary_width):
                l, r = bmesh.g1D[0][0 + ib], bmesh.g1D[0][-1 - ib]
                if not bmesh.reduced[0]:
                    # There are only boundary points on the left side if the axis is not reduced.
                    bp.append(np.nonzero(bmesh.grid == l)[0])
                else:
                    # Otherwise they are symmetry points,
                    # - which we mut keep separate for each boundary layer and each axis.
                    # - and we must avoid corners because they have no counterpart in the original mesh
                    # we remove the corners afterwards
                    self.sp_l[ib][0] = np.nonzero(bmesh.grid == l)[0]

                bp.append(np.nonzero(bmesh.grid == r)[0])

            self.bp_l = np.sort(np.unique(np.concat(bp)))  # ndarray with linear indices of the boundary points
            self.bp_xyz = bmesh.grid[self.bp_l]  # ndarray with coordinates    of the boundary points

            # determine the linear indices of the points in u corresponding to the symmetry points.
            self.sp_lu = [[None]*bmesh.dim for ib in range(self.boundary_width)]
            if bmesh.reduced[0]:
                for ib in range(self.boundary_width):
                    iu = self.boundary_width - ib -1
                    l = self.mesh.g1D[0][iu]
                    self.sp_lu[ib][0] = np.nonzero(mesh.grid == l)[0]

        # MeshQuantity for embedding the solution u
        self.v = MeshQuantity(bmesh, n_components=1, symmetry=1)
        self.v_ = self.v.data[:,0]

        # Create LinearOperator instance.
        self.LO = sparse.linalg.LinearOperator((mesh.linear_size, mesh.linear_size), matvec=self.matvec, dtype=np.float64)

    def set_symmetry(self, symmetry):
        """set the symmetry that the solution must obey."""
        self.v.set_symmetry(symmetry)

    def denominator(self, idim=0):
        """to be overridden by derived class
        Returns:
            the denominator of the stencil expression for the Laplacian
        """
        raise NotImplementedError()

    def Astar_v(self):
        """To be overridden by derived class (FD only)
        Used by matvec and assemble_rhs.

        Returns:
            the matrix-vector product of $A^*$ and `self.v`, the linear representation of
            the embedding grid.
        """
        raise NotImplementedError()

    def matvec(self, u):
        """to be overridden by derived class
        Embed the current iteration u in the embedding grid, take symmetries into account, and return
        the matrix-vector product of the current iteration u and the embedded solution.
        """
        raise NotImplementedError()

    def assemble_rhs(self, f, boundary_value):
        """Assemble the right hand side of the linear system:
          - apply scaling to rhs in case of uniform grid spacing
          - applying the Dirichlet boundary conditions to the rhs.

        Args:
            f (ndarray): right hand side of the generalized Poisson equation. Modified during
                `assemble` and `solve` methods
            boundary_value (Callable): function that yields the Dirichlet boundary value when
                applied to the boundary points of the mesh.
        """
        assert f.shape == (self.mesh.linear_size,)
        self.f = f

        # scale rhs f
        if self.h_uniform != 0:
            self.f *= self.denominator()

        # set all interior nodes to 0 (but it is more efficient to set all nodes to 0)
        self.v_[:] = 0.
        # set boundary values on boundary nodes
        if self.mesh.dim == 3:
            self.v_[self.bp_l] = boundary_value(self.bp_xyz[:,0], self.bp_xyz[:,1], self.bp_xyz[:,2])
        elif self.mesh.dim == 2:
            self.v_[self.bp_l] = boundary_value(self.bp_xyz[:,0], self.bp_xyz[:,1])
        else:
            self.v_[self.bp_l] = boundary_value(self.bp_xyz)

        # this is probably not the most efficient way to do this, but it is only called once,
        # and it has the tremendous advantage of being generic and avoiding corner cases.
        b = self.Astar_v()
        self.f -= b

        self.v_[self.bp_l] = 0.  # required by the iterative solvers when calling matvec

    def solve(self, method='cg', **kwargs):
        """Solve the generalized Poisson equation an iterative solver.
        Args:
            method: one of the methods in `iterative_linear_solvers`
            **kwargs: keyword arguments to be passed to the iterative solver.
        """
        iterative_linear_solver = iterative_linear_solvers[method]
        u = iterative_linear_solver(self.LO, self.f, **kwargs)
        return u

class L3p(FDStencil):
    """Classical 3 point 1D stencil for the Laplacian with coefficients [1,-2,1]/h**2."""
    def __init__(self, mesh):
        """
        Args:
            mesh: LagrangeMesh instance on which the solution to the Poisson equation is sougth.
        """
        super().__init__(mesh, boundary_width=1)

    def denominator(self, idim=0):
        """The denominator of the stencil expression for the Laplacian."""
        h = self.h_uniform if (self.h_uniform != 0) else self.v.mesh.d[idim]
        h = h**2
        return h

    def Astar_v(self):
        """
        Returns:
            the matrix-vector product of $A^*$ and `self.v`, the linear representation of
            the embedding grid.

        Used by matvec and assemble_rhs.
        """
        if self.mesh.dim == 3:
            if self.h_uniform:
                asv  = -6 * self.v.dataG[1:-1, 1:-1, 1:-1] \
                          +(self.v.dataG[ :-2, 1:-1, 1:-1] \
                          + self.v.dataG[2:  , 1:-1, 1:-1])\
                          +(self.v.dataG[1:-1,  :-2, 1:-1] \
                           +self.v.dataG[1:-1, 2:  , 1:-1])\
                          +(self.v.dataG[1:-1, 1:-1,  :-2] \
                           +self.v.dataG[1:-1, 1:-1, 2:  ])

            else:
                invh2_0, invh2_1, invh2_2 = [1 / self.denominator(i) for i in range(3)]

                asv = -2 * (invh2_0 + invh2_1 + invh2_2) *  self.v.dataG[1:-1, 1:-1, 1:-1] \
                       +    invh2_0                      * (self.v.dataG[ :-2, 1:-1, 1:-1] \
                                                           +self.v.dataG[2:  , 1:-1, 1:-1])\
                       +              invh2_1            * (self.v.dataG[1:-1,  :-2, 1:-1] \
                                                           +self.v.dataG[1:-1, 2:  , 1:-1])\
                       +                        invh2_2  * (self.v.dataG[1:-1, 1:-1,  :-2] \
                                                           +self.v.dataG[1:-1, 1:-1, 2:  ])

            return asv.reshape(self.mesh.linear_size, order='F')

        elif self.mesh.dim == 2:
            if self.h_uniform:
                # np.multiply(-4.0, self.v.dataG[bw  :-bw  , bw  :-bw  ], out=self.uG)
                asv  = -4 * self.v.dataG[1:-1, 1:-1] \
                          +(self.v.dataG[ :-2, 1:-1] \
                          + self.v.dataG[2:  , 1:-1])\
                          +(self.v.dataG[1:-1,  :-2] \
                           +self.v.dataG[1:-1, 2:  ])

            else:
                invh2_0, invh2_1 = [1 / self.denominator(i) for i in range(2)]

                asv = -2 * (invh2_0 + invh2_1) *  self.v.dataG[1:-1, 1:-1] \
                       +    invh2_0            * (self.v.dataG[ :-2, 1:-1] \
                                                 +self.v.dataG[2:  , 1:-1])\
                       +              invh2_1  * (self.v.dataG[1:-1,  :-2] \
                                                 +self.v.dataG[1:-1, 2:  ])

            return asv.reshape(self.mesh.linear_size, order='F')

        elif self.mesh.dim == 1:
            return -2 * self.v_[1:-1] \
                    +   self.v_[0:-2] \
                    +   self.v_[2:  ]

    def matvec(self, u):
        """Embed the current iteration u in the embedding grid, take symmetries into account, and return
        the matrix-vector product of the current iteration u and the embedded solution.
        """
        if self.mesh.dim >= 2:
            if self.mesh.dim == 3:
                self.v.dataG[1:-1, 1:-1, 1:-1, 0] = u.reshape(self.mesh.N, order='F')
            else: # self.mesh.dim == 2
                self.v.dataG[1:-1, 1:-1, 0] = u.reshape(self.mesh.N, order='F')

            for idim in range(self.mesh.dim):
                if self.mesh.reduced[idim]:
                    self.v_[self.sp_l[0][idim]] = self.v.symmetry[0,idim] * u[self.sp_lu[0][idim]]

        elif self.mesh.dim == 1:
            self.v_[1:-1] = u
            if self.mesh.reduced[0]:
                self.v_[0] = self.v.symmetry[0]*self.v_[1]

        return self.Astar_v()


class L5p(FDStencil):
    """Classical 5 point 1D stencil for the Laplacian with coefficients [-1, 16, -30, 16, -1]/12*h**2."""
    def __init__(self, mesh):
        """
        Args:
            v MeshQuantity on embedding mesh
        """
        super().__init__(mesh, boundary_width=2)

    def denominator(self, idim=0):
        """The denominator of the stencil expression for the Laplacian."""
        h = self.h_uniform if (self.h_uniform != 0) else self.v.mesh.d[idim]
        h = 12 * h**2
        return h

    def Astar_v(self):
        """
        Returns:
            the matrix-vector product of $A^*$ and `self.v`, the linear representation of
            the embedding grid.

        Used by matvec and assemble_rhs.
        """
        if self.mesh.dim == 3:
            if self.h_uniform:
                # np.multiply(-4.0, self.v.dataG[bw  :-bw  , bw  :-bw  ], out=self.uG)
                asv  =  -90 *  self.v.dataG[2:-2, 2:-2, 2:-2] \
                       + 16 * (self.v.dataG[3:-1, 2:-2, 2:-2] \
                              +self.v.dataG[1:-3, 2:-2, 2:-2])\
                            - (self.v.dataG[4:  , 2:-2, 2:-2] \
                              +self.v.dataG[ :-4, 2:-2, 2:-2])\
                       + 16 * (self.v.dataG[2:-2, 3:-1, 2:-2] \
                              +self.v.dataG[2:-2, 1:-3, 2:-2])\
                            - (self.v.dataG[2:-2, 4:  , 2:-2] \
                              +self.v.dataG[2:-2,  :-4, 2:-2])\
                       + 16 * (self.v.dataG[2:-2, 2:-2, 3:-1] \
                              +self.v.dataG[2:-2, 2:-2, 1:-3])\
                            - (self.v.dataG[2:-2, 2:-2, 4:  ] \
                              +self.v.dataG[2:-2, 2:-2,  :-4])

            else:
                invh2_0, invh2_1, invh2_2 = [1 / self.denominator(i) for i in range(3)]

                asv  =  -30 * (invh2_0 + invh2_1 + invh2_2) *  self.v.dataG[2:-2, 2:-2, 2:-2] \
                       + 16 *  invh2_0                      * (self.v.dataG[3:-1, 2:-2, 2:-2] \
                                                              +self.v.dataG[1:-3, 2:-2, 2:-2])\
                       -       invh2_0                      * (self.v.dataG[4:  , 2:-2, 2:-2] \
                                                              +self.v.dataG[ :-4, 2:-2, 2:-2])\
                       + 16 *             invh2_1           * (self.v.dataG[2:-2, 3:-1, 2:-2] \
                                                              +self.v.dataG[2:-2, 1:-3, 2:-2])\
                       -                  invh2_1           * (self.v.dataG[2:-2, 4:  , 2:-2] \
                                                              +self.v.dataG[2:-2,  :-4, 2:-2])\
                       + 16 *                      invh2_2  * (self.v.dataG[2:-2, 2:-2, 3:-1] \
                                                              +self.v.dataG[2:-2, 2:-2, 1:-3])\
                       -                           invh2_2  * (self.v.dataG[2:-2, 2:-2, 4:  ] \
                                                              +self.v.dataG[2:-2, 2:-2,  :-4])

            return asv.reshape(self.mesh.linear_size, order='F')

        elif self.mesh.dim == 2:
            if self.h_uniform:
                asv  =  -60 *  self.v.dataG[2:-2, 2:-2] \
                       + 16 * (self.v.dataG[3:-1, 2:-2] \
                              +self.v.dataG[1:-3, 2:-2])\
                            - (self.v.dataG[4:  , 2:-2] \
                              +self.v.dataG[ :-4, 2:-2])\
                       + 16 * (self.v.dataG[2:-2, 3:-1] \
                              +self.v.dataG[2:-2, 1:-3])\
                            - (self.v.dataG[2:-2, 4:  ] \
                              +self.v.dataG[2:-2,  :-4])

            else:
                invh2_0, invh2_1 = [1 / self.denominator(i) for i in range(2)]

                asv  =  -30 * (invh2_0 + invh2_1) *  self.v.dataG[2:-2, 2:-2] \
                       + 16 *  invh2_0            * (self.v.dataG[3:-1, 2:-2] \
                                                    +self.v.dataG[1:-3, 2:-2])\
                       -       invh2_0            * (self.v.dataG[4:  , 2:-2] \
                                                    +self.v.dataG[ :-4, 2:-2])\
                       + 16 *             invh2_1 * (self.v.dataG[2:-2, 3:-1] \
                                                    +self.v.dataG[2:-2, 1:-3])\
                       -                  invh2_1 * (self.v.dataG[2:-2, 4:  ] \
                                                    +self.v.dataG[2:-2,  :-4])

            return asv.reshape(self.mesh.linear_size, order='F')

        elif self.mesh.dim == 1:
            return  -30 *  self.v_[2:-2] \
                   + 16 * (self.v_[3:-1] \
                          +self.v_[1:-3])\
                        - (self.v_[4:  ] \
                          +self.v_[ :-4])

    def matvec(self, u):
        """Embed the current iteration u in the embedding grid, take symmetries into account, and return
        the matrix-vector product of the current iteration u and the embedded solution.
        """
        if self.mesh.dim == 3:
            self.v.dataG[2:-2, 2:-2, 2:-2, 0] = u.reshape(self.mesh.N, order='F')
            for idim in range(self.mesh.dim):
                if self.mesh.reduced[idim]:
                    self.v_[self.sp_l[0][idim]] = self.v.symmetry[0,idim] * u[self.sp_lu[0][idim]]
                    self.v_[self.sp_l[1][idim]] = self.v.symmetry[0,idim] * u[self.sp_lu[1][idim]]

        elif self.mesh.dim == 2:
            self.v.dataG[2:-2, 2:-2, 0] = u.reshape(self.mesh.N, order='F')
            for idim in range(self.mesh.dim):
                if self.mesh.reduced[idim]:
                    self.v_[self.sp_l[0][idim]] = self.v.symmetry[0,idim] * u[self.sp_lu[0][idim]]
                    self.v_[self.sp_l[1][idim]] = self.v.symmetry[0,idim] * u[self.sp_lu[1][idim]]

        elif self.mesh.dim == 1:
            self.v_[2:-2] = u
            if self.mesh.reduced[0]:
                self.v_[0] = self.v.symmetry[0] * self.v_[3]
                self.v_[1] = self.v.symmetry[0] * self.v_[2]

        return self.Astar_v()

class LLg:
    """
    Lagrange stencil, using the LagrangeMesh for computing the Laplacian and solving
    $(a\Delta + b)u = f$
    """
    def __init__(self, mesh, a, b, symmetry=1, initial_guess=0):
        """
        Args:
            mesh: LagrangeMesh on which to solve the generalized Poisson equation
            a,b: coefficients of the generalized Poisson equation
            symmetry: symmetry of the solution for reduced axes.
            initial_guess: initial guess for the solution, 0 implies u0 is zero everywhere.
                Alternatively, `'rhs'` can be specified which initializes u to the right hand
                side f (during assemble_rhs).
        """
        self.b = b
        self.a = a
        self.u = MeshQuantity(mesh, n_components=1, symmetry=symmetry)
        self.initial_guess = initial_guess
        if initial_guess == 0:
            self.u.data.fill(0)

        # Create LinearOperator instance.
        self.LO = sparse.linalg.LinearOperator(
            (self.u.mesh.linear_size, self.u.mesh.linear_size),
            matvec=self.matvec,
            dtype=np.float64
        )

    def set_symmetry(self, symmetry):
        """set the symmetry that the solution must obey."""
        self.u.set_symmetry(symmetry)

    def assemble_rhs(self, f, symmetry=None, boundary_value=0):
        """
        Args:
            f: MeshQuantity containing the right hand side of the equation
            boundary_value: ignored for LLg stencil
        """
        if symmetry is not None:
            self.set_symmetry(symmetry)

        if not boundary_value == 0:
            raise ValueError(f"LLg stencil is for solving the generalized Poisson equation."
                             f"Boundary conditions are Dirichlet and always 0. ")
        assert f.shape == (self.u.mesh.linear_size,)
        self.f = f
        self.f *= (1/self.a)

        if self.initial_guess == 'rhs':
            self.u.data[:,0] = self.f

    def solve(self, method='cg', **kwargs):
        """Solve the generalized Poisson equation an iterative solver.
        Args:
            method: one of the methods in `iterative_linear_solvers`
            **kwargs: keyword arguments to be passed to the iterative solver.
        """
        iterative_linear_solver = iterative_linear_solvers[method]
        u = iterative_linear_solver(self.LO, self.f, **kwargs)
        return u

    def matvec(self, u):
        self.u.data[:,0] = u
        Lu = self.u.differentiate(axes='xx', recompute=True)
        Lu += (self.b/self.a) * self.u.data
        return Lu[:,0]


class GeneralizedPoissonSolverMF:
    """Generalized Poisson Solver using a matrix-free Finite Differences approach on an
    embedded grid. Since it is matrix-free, it must use iterative solvers for solving
    the linear system. The embedded grid allows for a simple implementation of the symmetry
    boundary conditions, as it suffices to copy the current approximation of the solution
    to the embedding grid across the symmetry boundaries.
    Most of the actual work is done by the Stencil classes.
    """
    def __init__(self, mesh, stencil=L5p, a=1, b=0):
        """
        Args:
            mesh (LagrangeMesh): Lagrange mesh providing the grid on which the Generalized
                Poisson equation will be solved, using Finite Differences.
            stencil: A Stencil class,  Currently, `L3p` or `L5p`, with resp 3-point and 5-point discretization
                in each dimension. For the generalized Poisson equation, the stencil class must be `LLg`.
        """
        if a != 1 and b != 0 and stencil is not LLg:
            raise ValueError(f"The generalized Poisson problem require stencil to be Lagrange, not {stencil}.")

        if stencil is LLg:
            self.stencil = stencil(mesh, a, b)
        else:
            self.stencil = stencil(mesh)

    def assemble_rhs(self, f, boundary_value=0, symmetry=None):
        """Assemble the equations:
          - apply scaling to rhs in case of uniform grid spacing
          - applying the Dirichlet boundary conditions to the rhs.

        Args:
            f (ndarray): right hand side of the generalized Poisson equation. Modified during
                `assemble` and `solve` methods
            boundary_value (Callable): function that yields the Dirichlet boundary value when
                applied to the boundary points of the mesh. Ignored for LLg stencil.
        """
        if symmetry is not None:
            self.stencil.set_symmetry(symmetry)

        self.stencil.assemble_rhs(f, boundary_value)

    def solve(self, method='cg', **kwargs):
        """Solve the generalized Poisson equation with Dirichlet boundary conditions `bcs` using
        an the iterative solvers.
        Args:
            method: one of the methods in `iterative_linear_solvers`
        """
        return self.stencil.solve(method=method, **kwargs)


