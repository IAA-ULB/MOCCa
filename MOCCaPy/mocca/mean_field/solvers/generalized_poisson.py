"""
Solve the generalized Poisson equation

$$
( a \Delta + b ) F(r) = f(r)
$$

using the Finite Difference Method.
"""

import numpy as np
import scipy.sparse as sparse


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
        d = sparse.kron(    sparse.eye(n[2]), sparse.kron(    sparse.eye(n[1]), Laplacian1D(stencil, n[0]) )) \
          + sparse.kron(    sparse.eye(n[2]), sparse.kron(Laplacian1D(stencil, n[1]),     sparse.eye(n[0]) )) \
          + sparse.kron(Laplacian1D(stencil, n[2]), sparse.kron(    sparse.eye(n[1]),     sparse.eye(n[0]) ))
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

def apply_Dbc_L(L, rows, bc_scale):
    """
    Args:
        L: dia_array representing the Laplacian to which to apply Dirichlet boundary conditions,
            or column vector reprenting the right hand side.
    """
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


def apply_Dbc_f(f, F, rows, bc_scale):
    """
    Args:
        f:  column vector reprenting the right hand side.
    """

    assert isinstance(f, np.ndarray) and len(f.shape) == 1
    if bc_scale is None:
        for irow in rows:
            f.data[0] = F.data[0]  # Left  end of the domain
            f.data[-1] = F.data[-1]  # Right end of the domain
    else:
        f.data[0] = F.data[0] * bc_scale  # Left  end of the domain
        f.data[-1] = F.data[-1] * bc_scale  # Right end of the domain


def collect_boundary_conditions(F):
    """
    Args:
        F: MeshQuantity containing the Dirichlet boundary conditions on the boundary points.
    """

    mesh = F.mesh
    if mesh.ndim == 3:
        raise NotImplementedError
    elif mesh.ndim == 2:
        ix0, ix1 = mesh.M[0]//2, mesh.M[0]
        iy0, iy1 = mesh.M[1]//2, mesh.M[2]
        bps2 = []
        for ix in range(ix0, ix1):
            bps2.append([ix, iy1])
        for iy in range(iy0, iy1):
            bps2.append([ix1, iy])
        bps2.append([ix1, iy1])

        halfM0 = mesh.M[0]//2
        if mesh.reduced[0]:
            for ibp, bp in enumerate(bps2):
                bps2[ibp][0] -= halfM0
        else:
            mirror_bps2 = []
            for bp in bps2:
                mirror_bps2.append([bp[0]-halfM0, bp[1]])
            bps2.extend(mirror_bps2)

        halfM1 = mesh.M[1]//2
        if mesh.reduced[1]:
            for ibp, bp in enumerate(bps2):
                bps2[ibp][1] -= halfM1
        else:
            mirror_bps2 = []
            for bp in bps2:
                mirror_bps2.append([bp[0], bp[1]-halfM1])
            bps2.extend(mirror_bps2)
        bcs = [F.dataG[bp[0], bp[1]] for bp in bps2]
        bps = [bp[0] + mesh.N[0]*bp[1] for bp in bps2]

    elif mesh.dim == 1:
        bps = [mesh.N[0] - 1] if mesh.reduced[0] else \
              [0, mesh.N[0]-1]
        bcs = [F.data[bp] for bp in bps]

    return bps, bcs


class GeneralizedPoissonSolver:
    """
    This class solves the generalized Poisson equation
    """
    def __init__(self, f, a=1, b=0, stencil=3, bc_scale=None):
        """
        Args:
            f (MeshQuantity): right hand side of the generalized Poisson equation.
            a, b (float) parameters of the generalized Poisson equation.
            bc_scale: if None the boundary conditions are applied by setting all elements in L corresponding to
                boundary points to zero except for the diagonal. The rhs f is set to the boundary value. This
                destroys the symmetry of the matrix.
                Otherwise the diagonal elements of the boundary rows of the matrix are incremented by bc_scale
                (>>1) and the rhs entries is set equal to the boundary value times bc_scale.
        """
        self.f = f
        self.a = a # TODO: poisson -> generalized poisson
        self.b = b # TODO: poisson -> generalized poisson

        dim = self.f.mesh.dim
        N = self.f.mesh.N
        h = self.f.mesh.d
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
        # Apply the BCs to L
        self.bps = collect_boundary_points(self.f.mesh)
        apply_Dbc(self.L, self.bps, self.bc_scale)

    def assemble(self, F):
        """Solve the generalized Poisson equation with Dirichlet boundary conditions `bcs`.
        Args:
            F (MeshQuantity): on input: the Dirichlet boundary conditions on the boundary points of the mesh.
                on output: the solution of the generalized Poisson equation.
        """
        # Apply the BCs to the rhs f
        if self.f.mesh.dim == 3:
            raise NotImplementedError
        elif self.f.mesh.dim == 2:
            raise NotImplementedError
        elif self.f.mesh.dim == 1:
            if self.f.mesh.reduced[0]:
                raise NotImplementedError
            else:
                self.f.data *= self.h**2
                if self.bc_scale is None:
                    self.f.data[ 0] = F.data[ 0] # Left  end of the domain
                    self.f.data[-1] = F.data[-1] # Right end of the domain
                else:
                    self.f.data[ 0] = F.data[ 0] * self.bc_scale # Left  end of the domain
                    self.f.data[-1] = F.data[-1] * self.bc_scale # Right end of the domain

        else:
            raise NotImplementedError

    def solve(self):
        """Solve the generalized Poisson equation with Dirichlet boundary conditions `bcs`."""
        u = sparse.linalg.spsolve(self.L, self.f.data[:,0])
        return u