"""
Solve the generalized Poisson equation


$$
( a \Delta + b ) F(r) = f(r)
$$

using the Finite Difference Method.
"""

import numpy as np
import scipy.sparse as sparse
from tabulate import tabulate
import pytest
from mocca.mesh import LagrangeMesh, MeshQuantity


def Laplacian1D(n, stencil=3, h=None):
    """Construct 1D Laplacian matrix for a grid with n points.
    Args:
        n (int): dimension of the matrix.
        stencil (list): number of points in the stencil: 3 -> [1,-2,1], 5 -> [-1,16,-30,16,-1]
        h float: lattice spacing, if not None, the matrix is scaled by 1/h**2
    Returns:
        a diagonal matrix of the 1D Laplacian.
    """
    if isinstance(n, list):
        n = n[0]

    if stencil == 3:
        c = [-2,1]
    elif stencil == 5:
        c = [-30,16,-1]
    else:
        raise NotImplementedError(f'{stencil}-point 1D stencil not implemented.')

    if h is not None:
        # scale the coefficients
        inv_h2 = 1./h**2 if (stencil == 3) else \
                 1./(12 * h**2) # stencil == 5
        c = [ ci*inv_h2  for ci in c ]
    else:
        c = [ np.int8(ci) for ci in c]

    for i, ci in enumerate(c):
        if i == 0:
            offsets = [0]
            diagonals = [ci*np.ones(n, dtype=np.int8)]
        else:
            offsets.extend([i,-i])
            di = ci*np.ones(n, dtype=np.int8)
            diagonals.extend([di,di])

    d = sparse.diags_array(diagonals, offsets=offsets)
    return d


def Laplacian2D(n, stencil=3, h=None):
    """Construct 1D Laplacian matrix for a grid with n points.
    Args:
        stencil (list): number of points in the stencil: 3 -> [1,-2,1], 5 -> [-1,16,-30,16,-1]
        n (int): dimension of the matrix.
    Returns:
        a sparse matrix
    """
    if h is None:
        d = sparse.kron( sparse.eye(n[1])         , Laplacian1D(n[0], stencil), format='dia') \
          + sparse.kron(Laplacian1D(n[1], stencil),  sparse.eye(n[0])         , format='dia')
    else:
        d = sparse.kron( sparse.eye(n[1])               , Laplacian1D(n[0], stencil, h[0]), format='dia') \
          + sparse.kron(Laplacian1D(n[1], stencil, h[1]),  sparse.eye(n[0])               , format='dia')
    return d


def Laplacian3D(n, stencil=3, h=None):
    """Construct 1D Laplacian matrix for a grid with n points.
    Args:
        stencil (list): number of points in the stencil: 3 -> [1,-2,1], 5 -> [-1,16,-30,16,-1]
        n (int): dimension of the matrix.
    Returns:
        a sparse matrix.
    """
    if h is None:
        d = sparse.kron( sparse.eye(n[2])               , sparse.kron( sparse.eye(n[1])               , Laplacian1D(n[0], stencil)      , format='dia'), format='dia') \
          + sparse.kron( sparse.eye(n[2])               , sparse.kron(Laplacian1D(n[1], stencil)      ,  sparse.eye(n[0])               , format='dia'), format='dia') \
          + sparse.kron(Laplacian1D(n[2] , stencil)     , sparse.kron( sparse.eye(n[1])               ,  sparse.eye(n[0])               , format='dia'), format='dia')
    else:
        d = sparse.kron( sparse.eye(n[2])               , sparse.kron( sparse.eye(n[1])               , Laplacian1D(n[0], stencil, h[0]), format='dia'), format='dia') \
          + sparse.kron( sparse.eye(n[2])               , sparse.kron(Laplacian1D(n[1], stencil, h[1]),  sparse.eye(n[0])               , format='dia'), format='dia') \
          + sparse.kron(Laplacian1D(n[2], stencil, h[2]), sparse.kron( sparse.eye(n[1])               ,  sparse.eye(n[0])               , format='dia'), format='dia')
    return d

# TODO: dia_entry_set and dia_entry_add are no longer used. discard or improve.
# def dia_entry_set(diag, i_j_val):
#     """Set the [i,j] element of sparse 'dia'-format matrix `diag` to `val`.
#     This method is useful to apply the boundary conditions to a Laplacian matrix.
#
#     Args:
#         diag (scipy.sparse.diag_array): sparse matrix in 'dia' format
#         i_j_val: list of tuples with `(i,j,val)` : `diag[i,j] = val`
#
#     Raises:
#         ValueError if [i,j] is not on one of the diagonals in `diag`
#     """
#     offsets = diag.offsets.tolist()
#     for tpl in i_j_val:
#         i, j, val = tpl
#         id = j-i
#         try:
#             idd = offsets.index(id) # may raise ValueError
#         except ValueError:
#             raise ValueError(f"dia_array has no {id}-th diagonal.")
#         diag.data[idd,j] = val
#         pass
#
#
# def dia_entry_add(diag, i_j_val):
#     """Add val to the [i,j] element of sparse 'dia'-format matrix `diag` to `val`.
#     This method is useful to apply the boundary conditions to a Laplacian matrix.
#
#     Args:
#         diag (scipy.sparse.diag_array): sparse matrix in 'dia' format
#         i_j_val: list of tuples with `(i,j,val)` : `diag[i,j] += val`
#
#     Raises:
#         ValueError if [i,j] is not on one of the diagonals in `diag`
#     """
#     offsets = diag.offsets.tolist()
#     for tpl in i_j_val:
#         i, j, val = tpl
#         id = j - i
#         try:
#             idd = offsets.index(id)  # may raise ValueError
#         except ValueError:
#             raise ValueError(f"dia_array has no {id}-th diagonal.")
#         diag.data[id, j] += val
#         pass
#


class GeneralizedPoissonSolver:
    """
    This class solves the generalized Poisson equation
    """
    # TODO=done: implement reduced axes
    # TODO=done: test symmetric solutions      1D 2D 3D
    # TODO=done: test skew-symmetric solutions 1D 2D 3D
    # TODO=done: test non-uniform spacing         2D 3D
    # TODO: test 5 point stencil
    # TODO: implement generalized poisson
    # TODO: add iterative solution methods
    # TODO: time the methods

    def __init__(self, mesh, stencil=3, bc_scale=None):
        """
        This solver depends only on the mesh, the stencil and the way it treats Dirichlet boundary condition
        (`bc_scale`). The right hand side only appears when calling `assemble`.

        Args:
            f (MeshQuantity): right hand side of the generalized Poisson equation.
            a, b (float) parameters of the generalized Poisson equation.
            bc_scale: if None the boundary conditions are applied by setting all elements in L corresponding to
                boundary points to zero except for the diagonal. The rhs f is set to the boundary value. This
                destroys the symmetry of the matrix.
                Otherwise the diagonal elements of the boundary rows of the matrix are incremented by bc_scale
                (>>1) and the rhs entries is set equal to the boundary value times bc_scale.
        """
        self.mesh = mesh
        self.stencil = stencil

        dim = self.mesh.dim
        N = self.mesh.N
        h = self.mesh.d
        for hi in h[1:dim]:
            if h[0] != hi:
                self.uniform_h = False
                break
        else:
            # Lattice spacing is same in each direction, $h^2$ scaling will be applied to the rhs
            self.uniform_h = h[0]
            h = None # Tells the Laplacian that the h**2 scaling is applied to the rhs

        assert hasattr(self,'uniform_h')

        Laplacian = Laplacian3D if dim==3 else \
                    Laplacian2D if dim==2 else \
                    Laplacian1D if dim==1 else \
                    None

        self.L = Laplacian(N, stencil=stencil, h=h)
        assert self.L.format == 'dia', f"Unexpected format of Laplacian: '{self.L.format}, expected 'dia'."

        self.bc_scale = bc_scale

        # Apply the Dirichlet boundary conditions to L
        self.apply_Dbc_to_matrix()


    def assemble(self, f, boundary_value, symmetry=None):
        """Assemble the equations. This involves treating reduced axes and applying the Dirichlet
        boundary conditions to the rhs.

        Args:
            f (MeshQuantity): right hand side of the generalized Poisson equation.
            boundary_value (Callable): function that yields the Dirichlet boundary value when applied
                to the boundary points of the mesh.
            symmetry (list[int]): the symmetry components of the solution
        """
        assert f.mesh == self.mesh

        # treat reduced axes
        if any(f.mesh.reduced):
            # If we need to solve the Poisson equation for more problems with different symmetry,
            # it is better to keep a copy of L
            if not hasattr(self, 'L0'):
                self.L0 = self.L.copy()
            else:
                self.L = self.L0.copy()

            self.apply_symmetry_to_matrix(symmetry)

        if self.uniform_h:
            f.data[:,0] *= (     self.uniform_h**2) if (self.stencil == 3) else \
                           (12 * self.uniform_h**2)  # (self.stencil == 5)

        self.f = f
        self.apply_Dbc_to_rhs(boundary_value)

        self.L = self.L.tocsr()  # diagonal format raises SparseEfficiencyWarning


    def solve(self, method='direct'):
        """Solve the generalized Poisson equation with Dirichlet boundary conditions `bcs`."""
        if method == 'direct':
            u = sparse.linalg.spsolve(self.L, self.f.data[:,0])
        return u

    def apply_Dbc_to_matrix(self):
        """
        """
        try:
            rows = self.mesh.bp_l
        except AttributeError:
            raise AttributeError(
                f"LagrangeMesh instance has no `bp_l` member. You must call `mesh.collect_boundary_points()` first."
            )

        offsets = self.L.offsets.tolist()
        data = self.L.data
        nd = len(offsets)
        ld = data.shape[1]
        if self.bc_scale is None:
            # set 0-th diagonal to 1 and all others to zero
            for irow in rows:
                for io, o in enumerate(offsets):
                    if o == 0:
                        data[io, irow] = 1
                    elif 0 <= irow + o < ld:
                        data[io, irow + o] = 0
        else:
            id0 = offsets.index(0)
            for irow in rows:
                data[id0, irow] += self.bc_scale

    def apply_Dbc_to_rhs(self, boundary_value):
        """
        Args:
            f: column vector representing the right hand side.
            mesh: LagrangeMesh instance that provides `bpl` member: linear indices of boundary points. You must
                initialize mesh with `collect_boundary_points=True` or explicitly call `mesh.collect_boundary_points()`.
            boundary_value: a function that yields the boundary value at a point in space
            bc_scale: scaling factor for boundary conditions, specifying None yields exact boundary conditions,
                but breaks the symmetry of the matrix L. Specifying `bc_scale` >> 1 yields approximate boundary
                conditions, but keeps the symmetry of the matrix L. (Must be the same value as passed to
                `apply_Dbc_to_matrix`)
        """
        try:
            rows   = self.mesh.bp_l
            bp_xyz = self.mesh.bp_xyz
        except AttributeError:
            raise AttributeError(
                f"LagrangeMesh instance has no `bp_l` or `bp_xyz` member. You must call `mesh.collect_boundary_points()` first."
            )

        if self.mesh.dim == 3:
            g = boundary_value(bp_xyz[:,0], bp_xyz[:,1], bp_xyz[:,2])
        elif self.mesh.dim == 2:
            g = boundary_value(bp_xyz[:,0], bp_xyz[:,1])
        else:
            g = boundary_value(bp_xyz)

        self.f.data[rows,0] = g if (self.bc_scale is None) else \
                              g * self.bc_scale

    def apply_symmetry_to_matrix(self, symmetry):
        """
        Args:
            L: dia_array representing the Laplacian to which to apply Dirichlet boundary conditions,
                or column vector representing the right hand side.
            mesh: LagrangeMesh for which L was constructed
            symmetry (list[int]): the symmetry components of the solution
            stencil: stencil used for L, 3-point or 5-point per dimension
        """
        if symmetry is None:
            raise ValueError("The symmetry components of the solution must be specified.")

        if isinstance(symmetry, int):
            symmetry = [symmetry] * self.mesh.dim

        if isinstance(symmetry, (tuple, list)):
            for s in symmetry:
                if s not in  [1,-1]:
                    raise ValueError(f"The symmetry components of the solution must be 1 or -1, not {s}.")


        if self.mesh.dim > 1:
            # If h is uniform h**2 scaling is applied on the rhs, otherwise it is applied to the
            # matrix and the symmetry action on the matrix must be scaled as well

            sh = symmetry if self.uniform_h else \
                 [si / (     hi**2) for si, hi in zip(symmetry, self.mesh.d)] if (self.stencil == 3) else \
                 [si / (12 * hi**2) for si, hi in zip(symmetry, self.mesh.d)]  # (self.stencil == 5)

            offsets = self.L.offsets.tolist()
            md = len(offsets) // 2
            data = self.L.data
            if self.stencil == 3:  # only changes on the main diagonal
                for idim in range(3):
                    rows = self.mesh.sp_l[idim]
                    if rows is not None:
                        data[md, rows] += sh[idim] # only changes on the main diagonal

            elif self.stencil == 5:  # also changes off the main diagonal
                for idim in range(3):
                    rows = self.mesh.sp_l[idim]
                    n = 1                if (idim == 0) else \
                        n*self.mesh.N[0] if (idim == 1) else \
                        n*self.mesh.N[1]
                    if rows is not None:
                        data[md           , rows                             ] += sh[idim] * 16
                        data[md+(2*idim+1), rows + offsets[md+(2*idim+1)]    ] -= sh[idim]
                        data[md-(2*idim+1), rows + offsets[md+(2*idim+1)] - n] -= sh[idim]

        else: # mesh.dim == 1
            # h**2 scaling always applied to rhs.
            assert self.uniform_h

            s_x = symmetry[0]
            if self.stencil == 3:   # only changes on the main diagonal
                self.L.data[0,0] += s_x

            elif self.stencil == 5:
                self.L.data[0,0] += 16*s_x
                self.L.data[1,1] += -s_x # diagonal  1
                self.L.data[2,0] += -s_x # diagonal -1

            else:
                raise NotImplementedError()
        pass

    def __str__(self):
        """Readable presentation of the linear system to be solved.
        For debugging purposes mainly."""
        Lfull = self.L.toarray()
        f = self.f.data[:,0]
        n = Lfull.shape[0]
        tbl = []
        for i in range(n):
            line = [i]
            line.append(self.mesh.point_info(i))
            line.extend([Lfull[i,j] for j in range(n)])
            line.append(f[i])
            tbl.append(line)

        headers = ["i\j", "point"]
        headers.extend([str(j) for j in range(n)])
        headers.append("f")
        return tabulate(tbl, tablefmt="simple", headers=headers)

class GeneralizedPoissonSolverMF1D:
    """Generalized Poisson Solver using a matrix-free Finite Differences approach on an
    embedded grid. Since it is matrix-free, it must use iterative solvers for solviing
    the linear system. The embedded grid allows for a simple implementation of the symmetry
    boundary conditions, as it suffices to copy the current approximation of the solution
    to the embedding grid acros the symmetry boundaries.
    """
    iterative_linear_solvers = {
        'cg'    : sparse.linalg.cg,
        'cgs'   : sparse.linalg.cgs,
        'gmres' : sparse.linalg.gmres,
        'lgmres': sparse.linalg.lgmres,
        'minres': sparse.linalg.minres,
    }
    def __init__(self, mesh, stencil=3):
        """
        Args:
            mesh (LagrangeMesh): Lagrange mesh providing the grid on which the Generalized
                Poisson equation will be solved, using Finite Differences.
            stencil: stencil to be used for discretization. Currently, 3-point or 5-point in
                each dimension.
        """
        self.mesh = mesh
        self.stencil = stencil

        dim = mesh.dim
        self.boundary_width = \
            1 if (stencil == 3) else\
            2  # (stencil == 5)

        if dim == 3:
            raise NotImplementedError()

        elif dim == 2:
            raise NotImplementedError()

        elif dim == 1:
            # Set up embedding grid:

            self.grid = np.empty(mesh.N[0] + 2*self.boundary_width, dtype=np.float64)
            self.i0 = self.boundary_width
            self.i1 = self.boundary_width + mesh.N[0]
            self.grid[self.i0:self.i1] = self.mesh.grid
            self.h = self.mesh.d[0]
            for i in range(self.boundary_width):
                self.grid[ self.boundary_width-i-1] = self.grid[ self.boundary_width  ] - (i+1)*self.h
                self.grid[-self.boundary_width+i  ] = self.grid[-self.boundary_width-1] + (i+1)*self.h

            # array for embedding the solution u
            self.v = np.empty(mesh.N[0] + 2*self.boundary_width, dtype=np.float64)

            # Create LinearOperator instance.  Apparently, its ctor below calls matvec once
            # and thus any references it uses (c..q. s_x and Au) must be defined beforehand.
            # (It does that to detect the dtype of the result.)
            if self.mesh.reduced[0]:
                self.s_x = 1
            self.Au = np.empty(self.mesh.N[0], dtype=np.float64)
            self.LO = sparse.linalg.LinearOperator(
                (mesh.N[0], mesh.N[0]),
                matvec=self.matvec_3pt_1D if (self.stencil == 3) else \
                       self.matvec_5pt_1D,
            )

    def assemble(self, f, boundary_value, symmetry=None):
        """Assemble the equations:
          - applying the Dirichlet boundary conditions to the rhs.

        Args:
            f (ndarray): right hand side of the generalized Poisson equation. Modified during
                `assemble` and `solve` methosds
            boundary_value (Callable): function that yields the Dirichlet boundary value when
                applied to the boundary points of the mesh.
        """
        assert f.shape == (self.mesh.N[0],)
        self.f = f

        # scale rhs
        if self.stencil == 3:
            self.f *= self.h**2
        else:
            self.f *= 12 * self.h**2

        dim = self.mesh.dim
        if dim == 3:
            raise NotImplementedError()
        elif dim == 2:
            raise NotImplementedError()
        elif dim == 1:
            self.v = boundary_value(self.grid)
            # set all interior nodes to 0
            if self.mesh.reduced[0]:
                i0 = 0
                self.s_x = symmetry[0] if isinstance(symmetry, (tuple, list)) else symmetry
            else:
                i0 = self.i0
            self.v[i0:self.i1] = 0
            if self.stencil == 3:
                b = -2 * self.v[self.i0  :self.i1  ] \
                       + self.v[self.i0-1:self.i1-1] \
                       + self.v[self.i0+1:self.i1+1]
            elif self.stencil == 5:
                b =  -30 *  self.v[self.i0:self.i1    ] \
                    + 16 * (self.v[self.i0+1:self.i1+1] \
                           +self.v[self.i0-1:self.i1-1])\
                         - (self.v[self.i0+2:self.i1+2] \
                           +self.v[self.i0-2:self.i1-2])
            self.f -= b
            self.v[:] = 0 # required by matvec



    def matvec_3pt_1D(self, u):
        # Embed u with boundary region
        self.v[self.i0:self.i1] = u

        if self.mesh.reduced[0]:
            self.v[0] = self.s_x*u[0]

        self.Au[:] = -2 * self.v[self.i0  :self.i1  ] \
                      +   self.v[self.i0-1:self.i1-1] \
                      +   self.v[self.i0+1:self.i1+1]

        return self.Au

    def matvec_5pt_1D(self, u):
        # Embed u with boundary region
        self.v[self.i0:self.i1] = u

        if self.mesh.reduced[0]:
            self.v[0] = self.s_x*u[1]
            self.v[1] = self.s_x*u[0]

        self.Au[:] = -30 * self.v[self.i0  :self.i1  ] \
                     +16 *(self.v[self.i0+1:self.i1+1] \
                          +self.v[self.i0-1:self.i1-1])\
                     -    (self.v[self.i0+2:self.i1+2] \
                          +self.v[self.i0-2:self.i1-2])

        return self.Au

    def solve(self, method='cg'):
        """Solve the generalized Poisson equation with Dirichlet boundary conditions `bcs` using
        an the iterative solvers.
        Args:
            method: one of the methods in `self.iterative_linear_solvers`
        """
        ils = self.iterative_linear_solvers[method]
        u = ils(self.LO, self.f, rtol=1e-7)
        return u


class GeneralizedPoissonSolverMF:
    """Generalized Poisson Solver using a matrix-free Finite Differences approach on an
    embedded grid. Since it is matrix-free, it must use iterative solvers for solviing
    the linear system. The embedded grid allows for a simple implementation of the symmetry
    boundary conditions, as it suffices to copy the current approximation of the solution
    to the embedding grid across the symmetry boundaries.
    Most of the actual work is done by the Stencil classes.
    """
    def __init__(self, mesh, stencil=3):
        """
        Args:
            mesh (LagrangeMesh): Lagrange mesh providing the grid on which the Generalized
                Poisson equation will be solved, using Finite Differences.
            stencil: stencil to be used for discretization. Currently, 3-point or 5-point in
                each dimension.
        """
        # TODO: create an embedding grid in mesh.dim dimensions, in 'Fortran' order.
        # TODO: define indexing arrays for
        #       - interior points,
        #       - boundary points and
        #       - symmetry points for each reduced axis
        self.stencil = L3p(mesh) if (stencil == 3) else \
                       L5p(mesh)

    def assemble_rhs(self, f, boundary_value, symmetry=None):
        """Assemble the equations:
          - apply scaling to rhs in case of uniform grid spacing
          - applying the Dirichlet boundary conditions to the rhs.

        Args:
            f (ndarray): right hand side of the generalized Poisson equation. Modified during
                `assemble` and `solve` methods
            boundary_value (Callable): function that yields the Dirichlet boundary value when
                applied to the boundary points of the mesh.
        """
        self.stencil.assemble_rhs(f, boundary_value)

    def solve(self, method='cg', **kwargs):
        """Solve the generalized Poisson equation with Dirichlet boundary conditions `bcs` using
        an the iterative solvers.
        Args:
            method: one of the methods in `self.iterative_linear_solvers`
        """

        return self.stencil.solve(method=method, **kwargs)


class Stencil:
    """Base class for system of linear equations resulting from a particular stencil"""
    iterative_linear_solvers = {
        'cg'    : sparse.linalg.cg,
        'cgs'   : sparse.linalg.cgs,
        'gmres' : sparse.linalg.gmres,
        'lgmres': sparse.linalg.lgmres,
        'minres': sparse.linalg.minres,
    }

    def __init__(self, mesh, boundary_width):
        """
        Args:
            mesh: LagrangeMesh instance of grid that must be embedded
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

        # Collect boundary points and symmetry points
        def remove_corners(condition, idim):
            """Remove corners from symmetry points"""
            for jdim in range(self.mesh.dim):
                if jdim != idim:
                    L, R = bmesh.g1D[jdim][ self.boundary_width-1], \
                           bmesh.g1D[jdim][-self.boundary_width  ]
                    np.logical_and(condition, bmesh.grid[:,jdim] > L, out=condition)
                    np.logical_and(condition, bmesh.grid[:,jdim] < R, out=condition)
            return condition

        bp = []
        self.sp_l = [[None]*bmesh.dim] * self.boundary_width
        for ib in range(self.boundary_width):
            for idim in range(bmesh.dim):
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
        self.sp_lu = [[None]*bmesh.dim] * self.boundary_width
        for ib in range(self.boundary_width):
            for idim in range(bmesh.dim):
                iu = self.boundary_width - ib -1
                l = self.mesh.g1D[idim][iu]
                self.sp_lu[ib][idim] = np.nonzero(mesh.grid[:, idim] == l)[0]

        # MeshQuantity for embedding the solution u
        self.v = MeshQuantity(bmesh, n_components=1, symmetry=1)
        self.v_ = self.v.data[:,0]

        # Create LinearOperator instance.
        self.LO = sparse.linalg.LinearOperator((mesh.linear_size, mesh.linear_size), matvec=self.matvec, dtype=np.float64)



    def set_symmetry(self, symmetry):
        """set the symmetry that the solution must obey."""
        self.v.set_symmetry(symmetry)

    def denominator(self, idim=0):
        """to be overridden by derived class"""
        raise NotImplementedError()

    def Astar_v(self):
        """to be overridden by derived class"""
        raise NotImplementedError()

    def matvec(self, u):
        """to be overridden by derived class"""
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
        """Solve the generalized Poisson equation with Dirichlet boundary conditions `bcs` using
        an the iterative solvers.
        Args:
            method: one of the methods in `self.iterative_linear_solvers`
        """
        ils = self.iterative_linear_solvers[method]
        u = ils(self.LO, self.f, **kwargs)
        return u


class L3p(Stencil):
    """Classical 3 point stencil for the Laplacian with coefficients[1,-2,1]/h**2"""
    def __init__(self, mesh):
        """
        Args:
            v MeshQuantity on embedding mesh
        """
        super().__init__(mesh, boundary_width=1)

    def denominator(self, idim=0):
        h = self.h_uniform if (self.h_uniform != 0) else self.v.mesh.d[idim]
        h = h**2
        return h

    def Astar_v(self):
        if self.mesh.dim == 3:
            pass
        elif self.mesh.dim == 2:
            bw = self.boundary_width
            # if not hasattr(self, 'u'):
            #     self.uG = np.zeros(self.mesh.N, order='F')
            #     self.u  = self.uG.reshape(self.mesh.linear_size, order='F')
            if self.h_uniform:
                # np.multiply(-4.0, self.v.dataG[bw  :-bw  , bw  :-bw  ], out=self.uG)
                asv  = -4 * self.v.dataG[bw  :-bw  , bw  :-bw  ] \
                          + self.v.dataG[bw-1:-bw-1, bw  :-bw  ] \
                          + self.v.dataG[bw+1:     , bw  :-bw  ] \
                          +(self.v.dataG[bw  :-bw  , bw-1:-bw-1] \
                           +self.v.dataG[bw  :-bw  , bw+1:     ])
                return asv.reshape(self.mesh.linear_size, order='F')
            else:
                raise NotImplementedError()

        elif self.mesh.dim == 1:
            return -2 * self.v_[1:-1] \
                    +   self.v_[0:-2] \
                    +   self.v_[2:  ]

    def matvec(self, u):
        """Embed u with boundary region and copy symmetry points for reduced axes. Then apply Astar_v."""

        if self.mesh.dim >= 2:
            self.v.dataG[1:-1,1:-1,0] = u.reshape(self.mesh.N, order='F')
            for ib in range(self.boundary_width):
                for idim in range(self.mesh.dim):
                    if self.mesh.reduced[idim]:
                        self.v_[self.sp_l[ib][idim]] = self.v.symmetry[0,idim] * u[self.sp_lu[ib][idim]]

        elif self.mesh.dim == 1:
            self.v_[1:-1] = u
            if self.mesh.reduced[0]:
                self.v_[0] = self.v.symmetry[0]*self.v_[1]

        return self.Astar_v()

class L5p(Stencil):
    """Classical 5 point stencil for the Laplacian with coefficients[-1, 16, -30, 16, -1]/12*h**2"""
    def __init__(self, mesh):
        """
        Args:
            v MeshQuantity on embedding mesh
        """
        super().__init__(mesh, boundary_width=2)

    def denominator(self, idim=0):
        h = self.h_uniform if (self.h_uniform != 0) else self.v.mesh.d[idim]
        h = 12 * h**2
        return h

    def Astar_v(self):
        if self.mesh.dim == 3:
            pass
        elif self.mesh.dim == 2:
            pass
        elif self.mesh.dim == 1:
            return  -30 *  self.v_[2:-2] \
                   + 16 * (self.v_[3:-1] \
                          +self.v_[1:-3]) \
                        - (self.v_[0:-4] \
                          +self.v_[4:  ])

    def matvec(self, u):
        """Embed u with boundary region and copy symmetry points for reduced axes. Then apply Astar_v."""
        if self.mesh.dim == 3:
            pass
        elif self.mesh.dim == 2:
            pass
        elif self.mesh.dim == 1:
            self.v_[2:-2] = u
            if self.mesh.reduced[0]:
                self.v_[0] = self.v.symmetry[0] * self.v_[3]
                self.v_[1] = self.v.symmetry[0] * self.v_[2]

        return self.Astar_v()

