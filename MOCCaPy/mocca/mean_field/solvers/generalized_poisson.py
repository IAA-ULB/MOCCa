"""
Solve the generalized Poisson equation

$$
( a \Delta + b ) F(r) = f(r)
$$

using the Finite Difference Method.
"""

import numpy as np
import scipy.sparse as sparse

from src_heph.heph_symmetries import symmetry


def Laplacian1D(n, stencil=3, h=None):
    """Construct 1D Laplacian matrix for a grid with n points.
    Args:
        stencil (list): number of points in the stencil: 3 -> [1,-2,1], 5 -> [-1,16,-30,16,-1]
        n (int): dimension of the matrix.
        h float: lattice spacing, if not None, the matrix is scaled by 1/h**2
    Returns:
        a diagonal matrix of the 1D Laplacian.
    """
    if isinstance(n, (tuple, list)):
        n = n[0]

    if stencil == 3:
        c = [-2,1]
    elif stencil == 5:
        c = [30,16,-1]
    else:
        raise NotImplementedError(f'{stencil}-point 1D stencil not implemented.')

    c = [np.int8(ci) for ci in c]

    if not h is None:
        inv_h2 = 1. / h**2
        c = [ ci*inv_h2 for ci in c ]

    for i, ci in enumerate(c):
        if i == 0:
            offsets = [0]
            diagonals = [ci*np.ones(n, dtype=np.int8)]
        else:
            offsets.extend([1,-1])
            di = ci*np.ones(n-i, dtype=np.int8)
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
        d = sparse.kron( sparse.eye(n[2])          , sparse.kron( sparse.eye(n[1])         , Laplacian1D(n[0], stencil) ), format='dia') \
          + sparse.kron( sparse.eye(n[2])          , sparse.kron(Laplacian1D(n[1], stencil),  sparse.eye(n[0])          ), format='dia') \
          + sparse.kron(Laplacian1D(n[2] , stencil), sparse.kron( sparse.eye(n[1])         ,  sparse.eye(n[0])          ), format='dia')
    else:
        d = sparse.kron( sparse.eye(n[2])               , sparse.kron( sparse.eye(n[1])               , Laplacian1D(n[0], stencil, h[0]) )) \
          + sparse.kron( sparse.eye(n[2])               , sparse.kron(Laplacian1D(n[1], stencil, h[1]),  sparse.eye(n[0])                )) \
          + sparse.kron(Laplacian1D(n[2], stencil, h[2]), sparse.kron( sparse.eye(n[1])               ,  sparse.eye(n[0])                ))
    return d


def dia_entry_set(diag, i_j_val):
    """Set the [i,j] element of sparse 'dia'-format matrix `diag` to `val`.
    This method is useful to apply the boundary conditions to a Laplacian matrix.

    Args:
        diag (scipy.sparse.diag_array): sparse matrix in 'dia' format
        i_j_val: list of tuples with `(i,j,val)` : `diag[i,j] = val`

    Raises:
        ValueError if [i,j] is not on one of the diagonals in `diag`
    """
    offsets = diag.offsets.tolist()
    for tpl in i_j_val:
        i, j, val = tpl
        id = j-i
        try:
            idd = offsets.index(id) # may raise ValueError
        except ValueError:
            raise ValueError(f"dia_array has no {id}-th diagonal.")
        diag.data[idd,j] = val
        pass


def dia_entry_add(diag, i_j_val):
    """Add val to the [i,j] element of sparse 'dia'-format matrix `diag` to `val`.
    This method is useful to apply the boundary conditions to a Laplacian matrix.

    Args:
        diag (scipy.sparse.diag_array): sparse matrix in 'dia' format
        i_j_val: list of tuples with `(i,j,val)` : `diag[i,j] += val`

    Raises:
        ValueError if [i,j] is not on one of the diagonals in `diag`
    """
    offsets = diag.offsets.tolist()
    for tpl in i_j_val:
        i, j, val = tpl
        id = j - i
        try:
            idd = offsets.index(id)  # may raise ValueError
        except ValueError:
            raise ValueError(f"dia_array has no {id}-th diagonal.")
        diag.data[id, j] += val
        pass

def apply_Dbc_to_matrix(L, mesh, bc_scale):
    """
    Args:
        L: dia_array representing the Laplacian to which to apply Dirichlet boundary conditions,
            or column vector representing the right hand side.
        mesh: LagrangeMesh instance that provides `bpl` member: linear indices of boundary points. You must
            initialize mesh with `collect_boundary_points=True` or explicitly call `mesh.collect_boundary_points()`.
        bc_scale: scaling factor for boundary conditions, specifying None yields exact boundary conditions,
            but breaks the symmetry of the matrix L. Specifying bc_scale >> 1 yields approximate boundary
            conditions, but keeps the symmetry of the matrix L.
    """
    try:
        rows = mesh.bp_l
    except AttributeError:
        raise AttributeError(
            f"LagrangeMesh instance has no `bp_l` member. You must call `mesh.collect_boundary_points()` first."
        )

    offsets = L.offsets.tolist()
    data = L.data
    nd = len(offsets)
    ld = data.shape[1]
    if bc_scale is None:
        # set 0-th diagonal to 1 and all others to zero
        for irow in rows:
            for io ,o in enumerate(offsets):
                if o == 0:
                    data[io, irow] = 1
                elif 0 <= irow + o < ld:
                    data[io, irow + o] = 0
    else:
        id0 = offsets.index(0)
        for irow in rows:
            data[id0, irow] += bc_scale


def apply_Dbc_to_rhs(f, mesh, boundary_value, bc_scale):
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
        rows   = mesh.bp_l
        bp_xyz = mesh.bp_xyz
    except AttributeError:
        raise AttributeError(
            f"LagrangeMesh instance has no `bp_l` or `bp_xyz` member. You must call `mesh.collect_boundary_points()` first."
        )

    if mesh.dim == 3:
        g = boundary_value(bp_xyz[:,0], bp_xyz[:,1], bp_xyz[:,2])
    elif mesh.dim == 2:
        g = boundary_value(bp_xyz[:,0], bp_xyz[:,1])
    else:
        g = boundary_value(bp_xyz)

    # TODO: Fix this
    # assert isinstance(f, np.ndarray) and len(f.shape) == 1

    f[rows,0] = g if (bc_scale is None) else \
                g * bc_scale

def apply_symmetry_to_matrix(L, mesh, stencil, symmetry):
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

    if symmetry not in  [1,-1]:
        raise ValueError("The symmetry components of the solution must be 1 or -1.")

    if not isinstance(symmetry, list):
        symmetry = [symmetry]*mesh.dim

    if mesh.dim == 3:
        NotImplementedError

    elif mesh.dim == 2:
        NotImplementedError

    else: # mesh.dim == 1
        # L[0,0] += symmetry[0]
        s_x = symmetry[0]
        if stencil == 3:
            L.data[0,0] += s_x
        elif stencil == 5:
            L.data[0,0] += 16*s_x
            L.data[1,1] += 16-s_x # diagonal  1
            L.data[2,0] += 16-s_x # diagonal -1
        else:
            raise NotImplementedError(f"Stencil {stencil} is not recognized by PyMOCCa.")
    pass

class GeneralizedPoissonSolver:
    """
    This class solves the generalized Poisson equation
    """
    # TODO: test non-uniform spacing
    # TODO: implement reduced axes
    # TODO: test 5 point stencil
    # TODO: implement generalized poisson

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
                break
        else:
            # Lattice spacing is same in each direction, $h^2$ scaling will be applied to the rhs
            self.h = h[0]
            h = None

        Laplacian = Laplacian3D if dim==3 else \
                    Laplacian2D if dim==2 else \
                    Laplacian1D if dim==1 else \
                    None

        self.L = Laplacian(N, stencil=stencil, h=h)

        self.bc_scale = bc_scale

        # Apply the Dirichlet boundary conditions to L
        apply_Dbc_to_matrix(self.L, self.mesh, self.bc_scale)


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

            apply_symmetry_to_matrix(self.L, self.mesh, self.stencil, symmetry)

        if self.h is not None:
            f.data[:,0] *= self.h**2

        apply_Dbc_to_rhs(f.data, f.mesh, boundary_value, self.bc_scale)

        self.L = self.L.tocsr()  # diagonal format raises SparseEfficiencyWarning
        self.f = f


    def solve(self, method='direct'):
        """Solve the generalized Poisson equation with Dirichlet boundary conditions `bcs`."""
        if method == 'direct':
            u = sparse.linalg.spsolve(self.L, self.f.data[:,0])
        return u