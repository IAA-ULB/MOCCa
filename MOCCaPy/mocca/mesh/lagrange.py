import numpy as np
import numpy.typing as npt
from numba import guvectorize, float64

from mocca.mesh.mesh_quantity import MeshQuantity
from mocca.mesh.lagrange_function import lagrange_function

# TODO: consider https://numpy.org/doc/2.1/reference/generated/numpy.meshgrid.html#numpy-meshgrid
#       np.meshgrid returns a sequence with a item for each coordinate axis, whereas
#       create_mesh returns a single array
#       the main benefit it reduction in code size and readability, the performance improvement
#       relative to an *entire* MOCCaPy run is marginal.

def create_mesh(g1D):
    """Create a rectangular mesh (1D, 2D, 3D) from 1D arrays with coordinates.

    Args:
        g1D: list of 1D ndarrays with the grid points on the respective coordinate axes.
            `len(g1D)` is the dimension of the space (1D/2D/3D)
    Returns:
        an array (order='F)
        if `len(g1D)==1`, g1D[0] is returned
        if `len(g1D)==1`, gxy of shape (nx*ny,2) is returned (Fortran order)
        if `len(g1D)==3`, gxzy of shape (nx*ny*nz,3) is returned (Fortran order)
    """
    if len(g1D)==3:
        nx = g1D[0].size
        ny = g1D[1].size
        nz = g1D[2].size
        gxyz = np.empty((nx*ny*nz,3), dtype=np.float64, order='F')
        gxy = create_mesh(g1D[0:2])
        nxy = nx*ny
        for iz in range(nz):
            gxyz[iz * nxy : (iz + 1) * nxy, 0:2] = gxy
            gxyz[iz * nxy : (iz + 1) * nxy, 2] = g1D[2][iz]
        return gxyz

    elif len(g1D) == 2:
        nx = g1D[0].size
        ny = g1D[1].size
        gxy = np.empty((nx*ny,2), dtype=np.float64, order='F')
        for iy in range(ny):
            gxy[iy * nx : (iy + 1) * nx, 0] = g1D[0]
            gxy[iy * nx : (iy + 1) * nx, 1] = g1D[1][iy]
        return gxy
    else:
        return g1D[0]

def _split_axes(axes, Q):
    """Split axes in a part that is already computed and a part that still has to be computed.
    E.g. 'xyz' if 'yz' was already computed yields `0, 1, 'yz'`, where 0 is the axis of the
    derivative that still must be computed, 1 the order of that derivative (thus `0, 1` means
    d/dx), and 'yz' is the derivative that will be recycled.

    Args:
        axes: a sorted string of 'x'|'y'|'z' characters
        Q: the MeshQuantity being differentiated.

    Returns:
        (axis, order, reused_derivative):
            axis (int): differentiation axis
            order (int): differentiation order
            reused_derivative: the reused derivative
    """
    nx = axes.count('x')
    ny = axes.count('y')
    nz = axes.count('z')
    n  = len(axes)

    # How to find the optimal substring that can be reused?
    # - the cost of differentiating in a single direction is the same, irrespective of the order
    # - the cost of differentiating wrt y and z is higher than wrt x due to non-contiguous memory access
    # Hence, only cross derivatives must be checked
    # Since we set up the computation in a way that the lower order derivatives are computed first,

    if nx == n or \
       ny == n or \
       nz == n:
        # Not a cross derivative, nothing to reuse.
        # We need to reshape Q from a linear array over all the grid points to a
        # dimD array over the true grid to leverage np.einsum for computing the
        # differentiation matrix products.
        axis = 0 if nx else \
               1 if ny else \
               2
        order = n
        return axis, order, Q.dataG

    else:
        if nx:
            reused_axes = axes.replace('x', '')
            return 0, nx, Q.derivativesG[reused_axes]
            # we only need the symmetry in the x direction, which hasn't changed.

        else:
            # note that this case implies nx==0, ny!=0 and nz!=0.
            reused_axes = axes.replace('y', '')
            return 1, ny, Q.derivativesG[reused_axes]

class Mesh:
    """Base class for LagrangeMesh and Poisson mesh"""
    def __init__(self, M:int|tuple, d:int|float|tuple, dim:int=0, reduced:tuple|bool=True):
        """Construct a Mesh in 1, 2 or 3 dimensions.

        Args:
            dim (int): dimension of mesh. if not specified, dim is guessed as len(M), where M must be a tuple.
            M: tuple with number of grid points in the respective dimensions. Each component of M=2N must be
                even, yielding N grid points at both sides of the origin. If a single M is specified, M is the
                same in every direction.
                The distinction between M and N is consistent with  Ryssens et al, PHYSICAL REVIEW C 92, 064318
                (2015) section III.B Lagrange-mesh representation.
            d: spacing of points in the respective dimensions. If d is a float or an int, d is the same on each
                coordinate axis.
            reduced: Restrict the corresponding coordinate axis to the positive half-axis. This implies that all
                quantities represented on the mesh are either symmetric or skew-symmetric with respect to that
                axis. A single bool indicates that all coordinate axes are reduced. A tuple of bools can be used
                to reduce some of the axes (True) and others not (False).
        """
        if not (0 <= dim <= 3):
            raise ValueError(f"`dim` parameter violates 0<={dim=}<=3.")

        expected_dim = None
        tuples = [M, d, reduced]
        for tpl in tuples:
            if isinstance(tpl, tuple):
                if expected_dim is None:
                    expected_dim = len(tpl)
                else:
                    if not (expected_dim == len(tpl)):
                        raise ValueError("All tuple arguments must have the same length (=spatial dimensionality).")

        if expected_dim is None:
            # there were no tuples
            self.dim = dim if (dim > 0) else 3
        else:
            if not (0 <= expected_dim <= 3):
                raise ValueError(f"`dim` parameter violates 0<={dim=}<=3.")

            # expected_dim must correspond to dim, unless dim == 0
            if expected_dim != dim and dim != 0:
                raise ValueError(F"Parameter {dim=}, violatingting len(d|M|reduced|shift) != {dim}).")

            self.dim = expected_dim

        if not isinstance(M, tuple):
            if isinstance(M, int):
                M = self.dim * (M,)
            else:
                raise ValueError("M parameter must be an int or a tuple of ints.")

        for mi in M:
            if not (mi > 0) or \
               not (mi % 2 == 0):
                raise ValueError(f"M must be strictly positive and even (got {M=}).")

        self.M = M

        if isinstance(d, tuple):
            d = tuple([float(di) for di in d])
        else:
            if isinstance(d, (int,float)):
                d = self.dim * (float(d),)
            else:
                raise ValueError("d parameter must be an int, float or a tuple of int/floats (got {d=}).")

        for di in d:
            if not(di > 0):
                raise ValueError("d must be strictly positive (got {d=}).")

        self.d = d

        if not isinstance(reduced, tuple):
            if isinstance(reduced, bool):
                reduced = self.dim * (reduced,)
            else:
                raise ValueError(f"`reduced` parameter must be an bool or a tuple of bools (got {reduced=}).")

        self.reduced = reduced

        self.N = [ M//2 if r else M for M,r in zip(self.M,self.reduced)]

        # Compute dv, the integration volume per mesh point.
        self.dv = np.prod(self.d) * 2 ** self.reduced.count(True)

        # Compute the unreduced (!) box widths
        # The full (unreduced box width is needed by the plane wave base functions
        self.box_width = np.array([M*d for (M,d) in zip(self.M, self.d)])

    def _initialize_gridpoints(self, **kwargs):
        raise NotImplementedError

class LagrangeMesh(Mesh):
    """
    Attributes:
        g1D: list of ndarrays containing the grid points in each dimension.
    """
    def __init__(self, M:int|tuple, d:int|float|tuple, reduced:tuple|bool=True, dim:int=0,
                       shift:tuple|float=.0,
                       bc:str='antiperiodic',
                       highest_derivative_order:int=2,
                       collect_boundary_points:bool=False,
                 ) -> None:
        """Construct a Lagrange mesh in 1, 2 or 3 dimensions.

        Args:
            dim, M, d, reduced: see base class Mesh
            bc: boundary condition type. 'antiperiodic' or 'periodic'.
            shift: subtract shift from the grid points. If non-zero, the corresponding `reduced` entry must be
                `False`.
            highest_derivative_order: allow for differentiation up to this order. D-matrices are pre-constructed
                up to this order.
            name: optional name for the mesh.

        Raises:
            ValueError: in case of invalid choices.

        Remark:
            The `reduced` parameter is derived from symmetry considerations and may at some point - when the complexity
             of Hephaestos is taken into account - be replaced with a `Symmetry` object. For the time being, however,
             we content with explicitly indicating which coordinate axes must be 'reduced'
        """
        super().__init__(dim=dim, M=M, d=d, reduced=reduced)

        # Validate shift parameter
        if isinstance(shift, tuple):
            if not len(shift) == self.dim:
                raise ValueError(f"`shift` parameter violates `len(shift)==self.dim`: {len(shift)=} != {self.dim=}.")
        else:
            if isinstance(shift, float):
                shift = self.dim * (shift,)
            else:
                raise ValueError("`shift` parameter must be a float or a tuple of floats (got {shift=}).")

        for (r,s) in zip(self.reduced, shift):
            if (r and (s != 0)):
                raise ValueError("A nonzero shift cannot be applied when reduced is True.")
        self.shift = shift

        # validate boundary condition
        if not (bc in ['antiperiodic', 'periodic']):
            raise ValueError(f"`bc` parameter must be either 'antiperiodic' or 'periodic' (got {bc=}).")

        self.bc = bc
        # Related convenience attributes
        self.antiperiodic = bc == 'antiperiodic'
        self.periodic = not self.antiperiodic # since there are only 2 options.

        self._initalize_gridpoints()

        self._setup_D_matrices(highest_derivative_order)

        if collect_boundary_points:
            self.collect_boundary_points()

    def _initalize_gridpoints(self) -> None:
        # initialize grid points:
        start = self.dim*[.0]
        n_reduced = self.dim*[0]
        g1D = self.dim*[0]
        for idim, (M_i, shift_i, reduced_i, d_i) in enumerate(zip(self.M, self.shift, self.reduced, self.d)):
            if reduced_i:
                n_reduced[idim] = M_i // 2
                start[idim] = 0.5
            else:
                n_reduced[idim] = M_i
                start[idim] = -(M_i - 1) * 0.5 - shift_i / d_i

            g1D[idim] = np.linspace(start[idim], start[idim] + n_reduced[idim] - 1, n_reduced[idim])
            g1D[idim] *= d_i
            # print(f"{idim=} {g1D[idim]=}")

        self.g1D = g1D
        self.grid = create_mesh(self.g1D)

        if self.dim == 3:
            self.gridx = np.empty(n_reduced, order='F')
            self.gridy = np.empty(n_reduced, order='F')
            self.gridz = np.empty(n_reduced, order='F')
            for k in range(n_reduced[2]):
                for j in range(n_reduced[1]):
                    self.gridx[:,j,k] = g1D[0]
            for k in range(n_reduced[2]):
                for i in range(n_reduced[0]):
                    self.gridy[i,:,k] = g1D[1]
            for j in range(n_reduced[1]):
                for i in range(n_reduced[0]):
                    self.gridz[i,j,:] = g1D[2]

        elif self.dim == 2:
            self.gridx = np.empty(n_reduced, order='F')
            self.gridy = np.empty(n_reduced, order='F')
            for j in range(n_reduced[1]):
                self.gridx[:, j] = g1D[0]
            for i in range(n_reduced[0]):
                self.gridy[i, :] = g1D[1]

        else: # self.dim == 1
            self.gridx = g1D[0]
            # print(f"{self.gridx=}")

        self.mesh_shape   = self.gridx.shape
        self.linear_size = int(np.prod(self.mesh_shape))

    def collect_boundary_points(self):
        """Collect and store the boundary points of this mesh. Reduced axes have only boundary points
        on the right side. (Points are ordered in the Fortran way, i.e. left-most indices
        vary faster. Needs to be called only once

        Returns:
            bp_l, bp_ijk, bp_xyz:
                bp_l: linear index of the boundary points, shape (n,)
                bp_ijk: grid indices of  the boundary points, shape (n, self.dim)
                bp_xyz: coordinates of the boundary points, shape (n, self.dim)
        """
        if self.dim == 3:
            # Collect boundary points (as [i,j] indices) on reduced axes, then mirror if not reduced
            lx = self.M[0] // 2 - 1
            ly = self.M[1] // 2 - 1
            lz = self.M[2] // 2 - 1
            x_boundary = [[lx,iy,iz] for iz in range(lz+1) for iy in range(ly+1)]
            y_boundary = [[ix,ly,iz] for iz in range(lz+1) for ix in range(lx)  ]
            z_boundary = [[ix,iy,lz] for iy in range(ly)   for ix in range(lx)  ]
            boundary = []
            boundary.extend(x_boundary)
            boundary.extend(y_boundary)
            boundary.extend(z_boundary)
            if not self.reduced[0]:
                M = self.M[0]
                halfM = M // 2
                boundary = [ [bp[0] + halfM, bp[1], bp[2]] for bp in boundary ]
                # mirror x
                mirror   = [ [M - bp[0] - 1, bp[1], bp[2]] for bp in boundary]
                boundary.extend(mirror)

            if not self.reduced[1]:
                M = self.M[1]
                halfM = M // 2
                boundary = [ [bp[0], bp[1] + halfM, bp[2]] for bp in boundary ]
                # mirror x
                mirror   = [ [bp[0], M - bp[1] - 1, bp[2]] for bp in boundary]
                boundary.extend(mirror)

            if not self.reduced[2]:
                M = self.M[2]
                halfM = M // 2
                boundary = [ [bp[0], bp[1], bp[2] + halfM] for bp in boundary ]
                # mirror x
                mirror   = [ [bp[0], bp[1], M - bp[2] - 1] for bp in boundary]
                boundary.extend(mirror)

            self.bp_ijk = np.array(boundary, order='F')

            self.bp_l   = self.bp_ijk[:, 0] \
                        + self.bp_ijk[:, 1] * self.N[0] \
                        + self.bp_ijk[:, 2] * self.N[0] * self.N[1]
            self.bp_xyz = self.grid[self.bp_l, :]

        elif self.dim == 2:
            # Collect boundary points (as [i,j] indices) on reduced axes, then mirror if not reduced
            lx = self.M[0] // 2 - 1
            ly = self.M[1] // 2 - 1
            x_boundary = [[lx,iy] for iy in range(ly)]   # [lx, ly] not included
            y_boundary = [[ix,ly] for ix in range(lx+1)] # [lx, ly]     included
            boundary = []
            boundary.extend(x_boundary)
            boundary.extend(y_boundary)
            if not self.reduced[0]:
                M = self.M[0]
                halfM = M // 2
                boundary = [ [bp[0] + halfM, bp[1]] for bp in boundary ]
                # mirror x
                mirror   = [ [M - bp[0] - 1, bp[1]] for bp in boundary]
                boundary.extend(mirror)

            if not self.reduced[1]:
                M = self.M[1]
                halfM = M // 2
                boundary = [ [bp[0], bp[1] + halfM] for bp in boundary ]
                # mirror x
                mirror   = [ [bp[0], M - bp[1] - 1] for bp in boundary]
                boundary.extend(mirror)

            self.bp_ijk = np.array(boundary, order='F')
            self.bp_l   = self.bp_ijk[:,0] + self.N[0]*self.bp_ijk[:,1]
            self.bp_xyz = self.grid[self.bp_l, :]

        else: # self.dim == 1
            last = self.N[0] - 1
            self.bp_l   = np.array([   last,]) if self.reduced[0] else \
                          np.array([0, last,])
            self.bp_ijk = self.bp_l
            self.bp_xyz = self.g1D[0][self.bp_l]


    def __repr__(self):
        reduced = ''.join(['T' if r else 'F' for r in self.reduced])
        shape = 'x'.join([str(m) for m in self.mesh_shape ])
        return (f"<LagrangeMesh(Mesh)[{shape}={self.linear_size}, {reduced=}]>")

    # ---------------------------------------------------------------------------
    # Casting arrays between linear and mesh shape
    # (n_gridpoints, n_components) <--> (Nx, Ny, Nz, n_components
    #---------------------------------------------------------------------------
    def cast2grid(self, a):
        """Cast an array `a` to shape (*self.mesh_shape, n_components).
        Underlying data are not copied.

        Note that this casts an array af shape (self.linear_size,) into an array
        of shape (self.linear_size, 1)

        Raises:
            ValueError if a is not commensurate. I.e. its size is not a multiple
            of self.linear_size.
        """
        if a.size % self.linear_size != 0:
            raise ValueError(f"Array is not commensurate with shape {self.mesh_shape}")

        mesh_shape_a = (*self.mesh_shape, a.size // self.linear_size)
        return a if (a.shape == mesh_shape_a) else \
               a.reshape(mesh_shape_a, order='F')


    def cast2linear(self, a):
        """Cast an array a to shape (self.linear_size, n_components). No data are
        copied.

        Raises:
            ValueError if a is not commensurate. I.e. its size is not a multiple
            of self.linear_size.
        """
        if a.size % self.linear_size != 0:
            raise ValueError(f"Array is not commensurate with linear shape ({self.linear_shape},)")

        linear_shape_a = (self.linear_size, a.size // self.linear_size)
        return a if (a.shape == linear_shape_a) else \
               a.reshape(linear_shape_a, order='F')


    def check_observable_shape(self, Q):
        return 'mesh' if not (Q.data.shape == (*self.shape, Q.n_components)) else \
               'linear' if not (Q.data.shape == (int(np.prod(self.shape),), Q.n_components)) else \
               'unknown' # this shouldn't happen


    def is_commensurate(self, a: npt.NDArray) -> bool:
        """Determine if array a is commensurate with the mesh
        Args:
            Q: a quantity discretized on the grid.
        """
        return a.size % self.linear_size == 0

    # ---------------------------------------------------------------------------
    # functions on the mesh
    # ---------------------------------------------------------------------------
    def apply(self, function, kwargs=None, out=None):
        """Apply a function on the mesh, i.e. compute the function value on every grid point.

        Args:
            function: function f(x,y,z) to be applied on the mesh array. If f accepts additional
                parameters, a closure must be provided.
            kwargs: dict with the parameters to be passed to the function.
            out: if provided the result is stored there. (I expect providing `out` makes only
                sense if `function` itself accepts an `out` parameter.)

        Returns:
            a flat mesh array (self.flat_shape).

        Raises:
            ValueError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        self.gridx = self.cast2linear(self.gridx)
        if self.dim > 1:
            self.gridy = self.cast2linear(self.gridy)
            if self.dim > 2:
                self.gridz = self.cast2linear(self.gridz)
                args = (self.gridx, self.gridy, self.gridz)
            else:
                args = (self.gridx, self.gridy)
        else:
            args = (self.gridx, )

        if kwargs is None:
            kwargs = {}

        if out is None:
            return function(*args, **kwargs)
        else:
            try:
                function(*args, **kwargs, out=out)
            except:
                out = function(*args, **kwargs)
        return out


    def integrate(self, Q):
        """Compute the integral of a scalar quantity `q` on the mesh.

        Args:
            Q (MeshQuantity) : discretized on the grid.

        Returns:
            a scalar:

        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        n_components = Q.n_components
        result = np.zeros((n_components,), dtype=np.float64)
        # for ic in range(n_components):
        #     q = Q[ic]
        #     result[ic] = q.sum() * self.dv
        result = q.data.sum(axis=0)
        return result if n_components > 1 else result[0]

    # ---------------------------------------------------------------------------
    # Basis functions
    # ---------------------------------------------------------------------------
    def plane_wave(self, L:float, k:float, r:np.array):
        """Evaluate the plane wave basis function at with wave vector `k` at position `r`.

        According to Ryssens et al, PHYSICAL REVIEW C 92, 064318 (2015), eq(14)
        According to [eq 5.1] in https://github.com/IAA-nuclear/tantalus_full/blob/MOCCaPy/MOCCaPy/mocca/mesh/lagrange.pdf

        Args:
            L: box width
            k: wave vector, shape is `(self.dim,)`. values must be odd half integers ±1/2, ±3/2, ...
            r: position vector, shape is `(self.dim, nr), nr being the number of evaluation points.
            as_complex: if True, return a complex array. otherwise return two arrays with the real and imaginary parts,
                resp.

        Returns:
            array of shape `(nr,2)` with complex and imaginary parts
        """
        arg = (2*np.pi * k / L) * r
        result = np.empty((len(r),2), dtype=np.float64, order='F')
        result[:,0] = np.cos(arg).ravel()
        result[:,1] = np.sin(arg).ravel()
        result *= np.sqrt(1/L)
        return result

    def basis_function(self, ijk, r, ijk_sign=1):
        """Evaluate the plane wave basis function corresponding to grid point `x_ijk = [x_i,y_j,z_k]` 
        at position `r`. The wave vector `k` is related to `x_i` as `k = x_i/dx`.

        According to [eq 5.1] in https://github.com/IAA-nuclear/tantalus_full/blob/MOCCaPy/MOCCaPy/mocca/mesh/lagrange.pdf

        Args:
            r: array of coordinates at which to evaluate the plane wave basis function.
            ijk: grid point index. A D-tuple of grid indices with `D == self.dim'. Positive integers
                are preferred. For reduced axes mesh points on the negative axis are selected by setting the
                corresponding entry in ijk_sign to -1.
            ijk_sign: optional sign of the mesh point coordinates. Only necessary if an index is zero. (-0 would select
                the first mesh point on the positive axis, rather than on the negative axis. This is fixed by setting
                the corresponding entry in ijk_sign to -1.)
                (ignored for non-reduced axes).
        Returns:
            (data, symmetry)
            data: Array of values of the basis function at r.
        Raises:
            ValueError:
        """
        if not isinstance(r, np.ndarray):
            raise ValueError(f"Expecting a numpy array for `r parameter, got {type(r)}")

        if isinstance(ijk_sign,int):
            ijk_sign = self.dim * (ijk_sign,)

        if self.dim == 3:
            assert r.shape[1] == 3
            i = ijk[0]
            j = ijk[1]
            k = ijk[2]
            two_pi_K = np.array([ 2. * np.pi * ( self.g1D[0][i] if not self.reduced[0] else
                                                 self.g1D[0][i] * ijk_sign[0]
                                               ) / (self.box_width[0] * self.d[0])
                                , 2. * np.pi * ( self.g1D[1][j] if not self.reduced[0] else
                                                 self.g1D[1][j] * ijk_sign[1]
                                               ) / (self.box_width[1] * self.d[1])
                                , 2. * np.pi * ( self.g1D[2][k] if not self.reduced[0] else
                                                 self.g1D[2][k] * ijk_sign[2]
                                               ) / (self.box_width[2] * self.d[2])
                                ], dtype=np.float64, order='F')
            factor = np.sqrt(1/(self.box_width[0] * self.box_width[1] * self.box_width[2]))

        elif self.dim == 2:
            assert r.shape[1] == 2
            i = ijk[0]
            j = ijk[1]
            two_pi_K = np.array([ 2. * np.pi * ( self.g1D[0][i] if not self.reduced[0] else
                                                 self.g1D[0][i] * ijk_sign[0]
                                               ) / (self.box_width[0] * self.d[0])
                                , 2. * np.pi * ( self.g1D[1][j] if not self.reduced[0] else
                                                 self.g1D[1][j] * ijk_sign[1]
                                               ) / (self.box_width[1] * self.d[1])
                                ], dtype=np.float64, order='F')
            factor = np.sqrt(1/(self.box_width[0] * self.box_width[1]))

        elif self.dim == 1:
            if len(r.shape) == 1:
                r = r.reshape((r.shape[0], 1), order='F')
            assert r.shape[1] == 1
            # ijk == i
            i = ijk
            two_pi_K = np.array([ 2. * np.pi * ( self.g1D[0][i] if not self.reduced[0] else
                                                self.g1D[0][i] * ijk_sign[0]
                                               ) / (self.box_width[0] * self.d[0])
                                ], dtype=np.float64, order='F')
            factor = np.sqrt(1 / self.box_width[0])

        arg = r @ two_pi_K
        result = np.empty((r.shape[0],2), dtype=np.float64, order='F')
        np.cos(arg, out=result[:,0])
        np.sin(arg, out=result[:,1])
        result *= factor

        return result


    # ---------------------------------------------------------------------------
    # Differentiation
    # ---------------------------------------------------------------------------
    def _compute_D1(self, axis) -> None:
        """Compute the full matrix D1 for `axis` (whether it is reduced or not).

        Returns:
            The full matrix D1

        Remark:
            see eq 10.1 in lagrange.md. Note that the paper by Ryssens et al 2015 has a sign error.
        """
        # We first compute the full D1, also for the reduced case. Then we comput $E^{ll}$ and $E^{lr}$
        twoN = self.M[axis]
        D1 = np.empty((twoN,twoN), dtype=np.float64) # deliberately not order='F' for performance reasons
        alternating_sign = np.empty((twoN+1), dtype=np.float64)
        alternating_sign[::2] = 1
        alternating_sign[1::2] = -1
        # alternating_sign = [1, -1, 1, -1, ...]
        # In agreement with the Warning following equations 10.1 and 10.2 in lagrange.md, which demands
        # an extra sign flip because they use (i-j) in the sine arguments instead of (j-i)
        # i = np.linspace(0, twoN - 1, twoN)
        minus_i = -np.linspace(0, twoN - 1, twoN)
        # In agreement with the formulas in Ryssens et al 2015 (eq 18),
        # j is the row index and i the column index (D_ji)
        for j in range(twoN):
            # row = j-i
            row = minus_i+j # j-i
            row *= (np.pi / twoN) # twoN = 2N
            row[j] = .1 # avoid division by 0 in row[j] below
            row = (np.pi / (twoN * self.d[axis])) / np.sin(row)
            # division by 0 in row[j] yields `inf`, to be replaced by 0, (l'Hopitals rule)
            row[j] = .0
            if j % 2 == 0: # j is even
                # flip the sign of entries 1 3 5 ... i.e. those with odd i
                # since j is even, i-j is odd iff i is odd: (-1)^(i-j) is -1 for odd i
                row = row * alternating_sign[:twoN]
                # assert alternating_sign[0] == 1 # True by construction
            else:          # j is odd
                # flip the sign of entries 0 2 4 ... i.e. those with even i
                # since j is odd, i-j is odd iff i is even  and (-1)^(i-j) is -1 for even i
                row = row * alternating_sign[1:]  # flip sign of entries 0 2 4 ...
                # assert alternating_sign[1] == -1  # True by construction
            D1[j, :] = row

        return D1


    @property
    def highest_derivative_order(self):
        """Highest differentiation order achievable by this Mesh. This can be chosen when constructing a
        LagrangeMesh object."""
        return self.D.shape[1]


    def _setup_D_matrices(self, highest_derivative_order):
        """setup D matrices infrastructure"""
        self.D = np.empty((3, highest_derivative_order), dtype=np.ndarray)
        #   Stores the full D matrix if the axis is not reduced,
        #   and the reduced matrix sum D+E, otherwise

        self.DmE = np.empty((3, highest_derivative_order), dtype=np.ndarray)
        #   Stores nothing if the axis is not reduced,
        #   and the reduced matrix difference D-E otherwise.

        # Compute D1 for all the axes
        self.D[0, 0] = self._compute_D1(0)
        if self.dim > 1:
            if (self.M[1] == self.M[0]) and \
               (self.d[1] == self.d[0]) and \
               (self.reduced[1] == self.reduced[0]):
                # If M, d and reduced are the same for the y-axis as for the x-axis, it has the same D1 matrix as the x-axis
                self.D[1, 0] = self.D[0, 0]
            else:
                # otherwise, compute its D1-matrix
                self.D[1, 0] = self._compute_D1(1)

            if self.dim > 2:
                if (self.M[2] == self.M[0]) and \
                   (self.d[2] == self.d[0]) and \
                   (self.reduced[2] == self.reduced[0]):
                    # If M, d and reduced are the same for the z-axis as for the x-axis, it has the same D1 matrix as the x-axis
                    self.D[2, 0] = self.D[0, 0]
                elif (self.M[2] == self.M[1]) and \
                     (self.d[2] == self.d[1]) and \
                     (self.reduced[2] == self.reduced[1]):
                    # If M, d and reduced are the same for the z-axis as for the y-axis, it has the same D1 matrix as the xy-axis
                    self.D[2, 0] = self.D[1, 0]
                else:
                    # otherwise, compute its D1-matrix
                    self.D[2, 0] = self._compute_D1(2)

        # Higher order D matrices
        for axis in range(self.dim):
            # compute the D2, D3,... matrices
            for order in range(1, highest_derivative_order):
                self.D[axis, order] = self.D[axis, 0] @ self.D[axis, order - 1]

            # for reduced axes derive Dlr and Ell from the higher order D matrices
            if self.reduced[axis]:
                N = self.M[axis]//2
                for order in range(highest_derivative_order):
                    D = np.empty((N, N), dtype=np.float64)
                    E = np.empty((N, N), dtype=np.float64)

                    # Copy the columns of Dlr into D
                    # Reverse the columns of the lower left quadrant Dll and copy into E
                    D[:,:] = self.D[axis, order][N:,N:] # copy Dlr (lower right quadrant)
                    for icol in range(N):
                        E[:,N-1-icol] = self.D[axis, order][N:,  icol] # the 1st column of Dll becomes the last column of E
                                                                       # the 2nd column of Dll becomes the second last column of E
                                                                       # ...
                    self.D  [axis, order] = D + E # replaces the full D matrix
                    self.DmE[axis, order] = D - E


    def _get_D(self, axis:int, order:int):
        """Get the pre-computed D matrix for (non-reduced) coordinate axis `axis`,
        for differentiation order `order`.
        """
        return self.D[axis, order - 1]


    def _get_DpE(self, axis:int, order:int):
        """Get the pre-computed D+E matrix for (reduced) coordinate axis `axis`,
        for differentiation order `order`.
        """
        return self.D[axis, order - 1] # `self.D` is Not a typo: the original full D matrix was replaced by DpE


    def _get_DmE(self, axis:int, order:int):
        """Get the pre-computed D-E matrix for (reduced) coordinate axis `axis`,
        for differentiation order `order`.
        """
        return self.DmE[axis, order - 1]


    def differentiate(self, Q, axes, out):
        """Differentiate MeshQuantity Q wrt axes.

        Args:
            Q: MeshQuantity to be differentiated.
            axes: Axes to be differentiated.
                A single derivative can be requested as a `str` combining the characters 'x', 'y', 'z'. E.g. the
                 `str` 'xyz' requests d^3/dxdydz. Multiple derivatives can be requested too:
                 In principle the order of the 'xyz' characters in axes is immaterial. Axes is sorted and repeating
                 characters are combined. E.g. 'xyx' -> 'xxy' is computed as D2x x D1y x Q.
        Returns:
            A `ndarray` of the same shape as Q.data containing the requested derivative.

        Note that
        - Q.data must have shape (*self.mesh_shape, n_components) be in mesh_shape and that the caller is responsible for this.
        - this is a low-level method, to be called by MeshQuantity.differentiate()

        """
# Handle single entries:

        # Find out if we can reuse a previously computed derivative as a starting point. E.g.:
        # - axes = 'xyz', if 'z' is already computed, use it as a starting point and apply D1x*D1y to dQdz
        # - axes = 'xyz', if 'y' is already computed, use it as a starting point and apply D1x*D1z to dQdy
        # - axes = 'xyz', if 'yz' is already computed, use it as a starting point and apply D1x to d2Qdydz
        # - axes = 'xx',  if 'x' is already computed, do NOT it as a starting point and apply D2x  to Q
        # We prefer to reuse component with 'y' and 'z', because the einsum operations imply non-contiguous

        axis, order, op2 = _split_axes(axes, Q)

        # Symmetry considerations.
        # The scheme for reusing derivatives is as follows:
        # - if nothing is reused the derivative is applies as Dnx * Dny * Dnz * Q
        # - if some precomputed derivative is reused, this is
        #   - either: Dnx * Dny * Q_nz with Q_nz = d^nzQ/dz^nz
        #   - or    : Dnx * Q_nynz  with Q_nynz = d^(ny+nz)Q/dy^ny.dz^nz
        # Hence none of the matrix multiplications acts on a precomputed axis. As a consequence the symmetry
        # sign of the part acted on is simply given by Q.symmetry.
        # The symmetry of the result is given by :
        #   (-1)^nx * Q.symmetry[0,:]
        #   (-1)^ny * Q.symmetry[1,:]
        #   (-1)^nz * Q.symmetry[2,:]

        if axis == 0: # x-axis
            subscripts = 'il,ljkq' if (self.dim == 3) else \
                         'il,ljq'  if (self.dim == 2) else \
                         'il,lq'   #  (self.dim == 1)
            if not self.reduced[axis]:
                D = self._get_D(axis=0, order=order)
                np.einsum(subscripts, D, op2, out=out, order='F', optimize=True)
            else:
                subscripts = subscripts[:-1] # drop the trailing 'q`, we're updating one component at a time.
                for iq in range(Q.n_components):
                    DE = self._get_DpE(axis=0, order=order) if Q.symmetry[iq,0] == 1 else \
                         self._get_DmE(axis=0, order=order)
                    op2_iq, out_ = (op2[:, :, :, iq], out[:,:,:,iq]) if (self.dim == 3) else \
                                   (op2[:, :, iq]   , out[:,:,iq]  ) if (self.dim == 2) else \
                                   (op2[:, iq]      , out[:,iq]    )
                    np.einsum(subscripts, DE, op2_iq, out=out_, order='F', optimize=True)

        elif self.dim >=1 and axis == 1: # y-axis
            subscripts = 'jl,ilkq' if (self.dim == 3) else \
                         'il,jlq'  #  (self.dim == 2)
            if not self.reduced[axis]:
                D = self._get_D(axis=1, order=order)
                np.einsum(subscripts, D, op2, out=out, order='F', optimize=True)
            else:
                subscripts = subscripts[:-1] # drop the trailing 'q`, we're updating one component at a time.
                for iq in range(Q.n_components):
                    DE = self._get_DpE(axis=1, order=order) if Q.symmetry[iq,1] == 1 else \
                         self._get_DmE(axis=1, order=order)
                    op2_iq, out_ = (op2[:, :, :, iq], out[:,:,:,iq]) if (self.dim == 3) else \
                                   (op2[:, :, iq]   , out[:,:,iq]) #  (self.dim == 2)
                    np.einsum(subscripts, DE, op2_iq, out=out_, order='F', optimize=True)

        elif axis == 2 : # z-axis (self.dim == 3 obviously)
            subscripts = 'kl,ijlq'
            if not self.reduced[axis]:
                D = self._get_D(axis=2, order=order)
                np.einsum(subscripts, D, op2, out=out, order='F', optimize=True)
            else:
                subscripts = subscripts[:-1] # drop the trailing 'q`, we're updating one component at a time.
                for iq in range(Q.n_components):
                    DE = self._get_DpE(axis=2, order=order) if Q.symmetry[iq,2] == 1 else \
                         self._get_DmE(axis=2, order=order)
                    np.einsum(subscripts, DE, op2[:,:,:,iq], out=out[:,:,:,iq], order='F')

    # ---------------------------------------------------------------------------
    # Interpolation methods
    #----------------------------------------------------------------------------
    # TODO: provide out= parameter?
    def interpolate(self, Q:MeshQuantity, r:npt.NDArray, algo='new', out=None) -> npt.NDArray:
        """Interpolate a quantity `Q` on the mesh.

        Args:
            Q: representation of a scalar quantity on the grid.
            r: array of `p` points at which to interpolate the scalar quantity `q`. 'r.shape == (p, self.dim)'
        """
        # TODO (?) speed this stuff up... Performance may be wrecked by numerous nested python loops.

        if algo == 'new':
            return self._interpolateND_new(Q, r, out=out)
        else:
            # TODO: obsolete? _interpolateND_new is generic and much faster
            if self.dim == 3:
                self.gridx = self.cast2linear(self.gridx)
                self.gridy = self.cast2linear(self.gridy)
                self.gridz = self.cast2linear(self.gridz)
                return self._interpolate3D(Q, r)
            elif self.dim == 2:
                self.gridx = self.cast2linear(self.gridx)
                self.gridy = self.cast2linear(self.gridy)
                return self._interpolate2D(Q, r)
            else:
                self.gridx = self.cast2linear(self.gridx)
                return self._interpolate1D(Q, r)


    def _interpolate1D(self, Q:MeshQuantity, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 1D  mesh."""
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.zeros(shape=(nr,nq), dtype=np.float64, order='F')
        for ig in range(self.linear_size):
            for iq in range(nq):
                sign = Q.sign(iq)
                fx = self.lagrange_function(r, ig, sign)
                Qr[:,iq] += Q.data[ig, iq] * fx

        return Qr


    def _interpolate2D(self, Q:MeshQuantity, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 2D mesh.

        Args:
            Q: A MeshQuantity with values specified on all grid points
            r: array of points at which to interpolate Q. 'r.shape == (nr, self.dim)'
        Returns:
            An array with the interpolated values.
        """
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.zeros(shape=(nr,nq), dtype=np.float64, order='F')
        ig = 0
        Nx = len(self.g1D[0])
        Ny = len(self.g1D[1])
        for iy in range(Ny):
            for ix in range(Nx):
                for iq in range(nq):
                    sign = Q.sign(iq)
                    fxy = self.lagrange_function(r, (ix,iy), sign)
                    Qr[:,iq] += Q.data[ig, iq] * fxy
                ig += 1

        return Qr

    def _interpolateND_new(self, Q:MeshQuantity, r:npt.NDArray, out):
        """Interpolate `Q` on a ND mesh.

        Args:
            Q: A MeshQuantity with values specified on all grid points
            r: array of points at which to interpolate Q. 'r.shape == (nr, self.dim)'
        """
        pi_over_Delta = [np.pi / d for d in self.d]
        inv2N         = [1. / M    for M in self.M]
        g1D = self.g1D
        f = [np.empty(self.mesh_shape[idim], dtype=np.float64) for idim in range(self.dim)]
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.empty((nr,nq), dtype=np.float64, order='F') if (out is None) else \
             out
        R = r if (len(r.shape) > 1) else \
            r.reshape((nr,1))

        for iq in range(nq):
            h_ijk = Q.dataG[:,:,:,iq] if (self.dim == 3) else \
                    Q.dataG[:,:,  iq] if (self.dim == 2) else \
                    Q.dataG[:,    iq]
            lifs = [lif if not self.reduced[idim] else
                    lif_symm if Q.symmetry[iq,idim] == 1 else
                    lif_skew for idim in range(self.dim)]
            for p in range(nr): # loop over all points to be interpolated
                rp = R[p,:]
                for idim in range(self.dim):
                    lifs[idim](g1D[idim], rp[idim], pi_over_Delta[idim], inv2N[idim], out=f[idim])

                Qr[p,iq] = np.einsum('ijk,i,j,k', h_ijk, f[0], f[1], f[2]) if (self.dim == 3) else \
                           np.einsum('ij,i,j'   , h_ijk, f[0], f[1]      ) if (self.dim == 2) else \
                           np.einsum('i,i'      , h_ijk, f[0]            )

        return Qr

    def _interpolate3D(self, Q:MeshQuantity, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 3D mesh.

        Args:
            Q: A MeshQuantity with values specified on all grid points
            r: array of points at which to interpolate Q. 'r.shape == (nr, self.dim)'
        """
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.zeros(shape=(nr, nq), dtype=np.float64, order='F')
        ig = 0
        Nx = len(self.g1D[0])
        Ny = len(self.g1D[1])
        Nz = len(self.g1D[2])
        for iz in range(Nz):
            for iy in range(Ny):
                for ix in range(Nx):
                    for iq in range(nq):
                        sign = Q.sign(iq)
                        fxy = self.lagrange_function(r, (ix, iy,iz), sign)
                        Qr[:, iq] += Q.data[ig, iq] * fxy
                    ig += 1
        return Qr

    # TODO: formulate in terms of lif/lif_symm/lif_skew?
    def lagrange_function(self, x:npt.NDArray|float, ijk:tuple|int, sign=1):
        """Evaluate the Lagrange function corresponding to the `ijk` grid point at `x`.
        Args:
            x : array of points at which to evaluate the Lagrange function.
            ijk : index of the grid point selecting the gridpoint of the Lagrange function.
            sign : the symmetry ±1. For reduced axes the lagrange function of negative grid points
                is multplied by `sign`.
        Returns:
            1D Array with values of the Lagrange function corresponding to the `ijk` grid point at `x`.
            If an axis is reduced, several arrays are returned, one for each combination of positeve
            and negative grid points.
            reduced == [True,True,True] -> + + +   [True,True] -> + +
                                           + + -                  + -
                                           + - +                  - +
                                           + - -                  - -
                                           - + +        [True] -> +
                                           - + -                  -
                                           - - +
                                           - - -
            If one of the axes is not reduced the rows with a - on that axes are removed
        """
        dim = len(ijk)
        assert dim == self.dim

        if self.dim == 3:

            i,j,k = ijk
            x_i, y_j, z_k = self.g1D[0][i], self.g1D[1][j], self.g1D[2][k]

            if self.reduced[0]:
                if sign[0] == 1:
                    fx = ( lagrange_function(x[:, 0],  x_i, self.d[0], self.M[0])
                         + lagrange_function(x[:, 0], -x_i, self.d[0], self.M[0]) )
                else:
                    fx = ( lagrange_function(x[:, 0],  x_i, self.d[0], self.M[0])
                         - lagrange_function(x[:, 0], -x_i, self.d[0], self.M[0]) )

                if self.reduced[1]:
                    if sign[1] == 1:
                        fy = ( lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                             + lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) )
                    else:
                        fy = ( lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                             - lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) )

                    if self.reduced[2]:
                        if sign[2] == 1:
                            fz = ( lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2])
                                 + lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) )
                        else:
                            fz = ( lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2])
                                 - lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) )

                    else:
                        fz = lagrange_function(x[:, 2], z_k, self.d[2], self.M[2])
                else:
                    fy = lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                    if self.reduced[2]:
                        if sign[2] == 1:
                            fz = ( lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2])
                                 + lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) )
                        else:
                            fz = ( lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2])
                                 - lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) )

                    else:
                        fz = lagrange_function(x[:, 2], z_k, self.d[2], self.M[2])
            else:
                fx = lagrange_function(x[:, 0],  x_i, self.d[0], self.M[0])
                
                if self.reduced[1]:
                    if sign[1] == 1:
                        fy = ( lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                             + lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) )
                    else:
                        fy = ( lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                             - lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) )

                    if self.reduced[2]:
                        if sign[2] == 1:
                            fz = ( lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2])
                                 + lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) )
                        else:
                            fz = ( lagrange_function(x[:, 2],  z_k, self.d[2], self.M[2])
                                 - lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]) )

                    else:
                        fz = lagrange_function(x[:, 2], z_k, self.d[2], self.M[2])
                else:
                    fy = lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                    if self.reduced[2]:
                        if sign[2] == 1:
                            fz = (lagrange_function(x[:, 2], z_k, self.d[2], self.M[2])
                                  + lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]))
                        else:
                            fz = (lagrange_function(x[:, 2], z_k, self.d[2], self.M[2])
                                  - lagrange_function(x[:, 2], -z_k, self.d[2], self.M[2]))

                    else:
                        fz = lagrange_function(x[:, 2], z_k, self.d[2], self.M[2])

            return fx * fy * fz 

        elif self.dim == 2:
            i, j = ijk
            x_i, y_j = self.g1D[0][i], self.g1D[1][j]

            if self.reduced[0]:
                if sign[0] == 1:
                    fx = ( lagrange_function(x[:, 0],  x_i, self.d[0], self.M[0])
                         + lagrange_function(x[:, 0], -x_i, self.d[0], self.M[0]) )
                else:
                    fx = ( lagrange_function(x[:, 0],  x_i, self.d[0], self.M[0])
                         - lagrange_function(x[:, 0], -x_i, self.d[0], self.M[0]) )

                if self.reduced[1]:
                    if sign[1] == 1:
                        fy = ( lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                             + lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) )
                    else:
                        fy = ( lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                             - lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) )
                else:
                    fy = lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
            else:
                fx = lagrange_function(x[:, 0],  x_i, self.d[0], self.M[0])
                if self.reduced[1]:
                    if sign[1] == 1:
                        fy = ( lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                             + lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) )
                    else:
                         fy = ( lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])
                              - lagrange_function(x[:, 1], -y_j, self.d[1], self.M[1]) )
                else:
                    fy = lagrange_function(x[:, 1],  y_j, self.d[1], self.M[1])

            return fx * fy

        else:
            i = ijk
            x_i = self.gridx[i]

            if self.reduced[0]:
                fx = ( lagrange_function(x,  x_i, self.d[0], self.M[0])
                     + lagrange_function(x, -x_i, self.d[0], self.M[0]) ) if sign[0] == 1 else \
                     ( lagrange_function(x,  x_i, self.d[0], self.M[0])
                     - lagrange_function(x, -x_i, self.d[0], self.M[0]) )
            else:
                fx =   lagrange_function(x,  x_i, self.d[0], self.M[0])

            return fx

@guvectorize([(float64[:], float64, float64, float64, float64[:])], '(n),(),(),()->(n)')
def lif(ui, u, pi_over_Delta, inv2N, out):
    """Lagrange interpolation function on not-reduced mesh. Used by `LagrangeMesh._interpolateND_new()`
    Args:
         ui: array with $\frac{\pi}{\Delta}u_i$, i=1..N corresponding to the grid points on the u-axis
         u : $\frac{\pi}{\Delta}u$ interpolation point on u-axis
         inv2N: $\frac{1}{2N}$
    """
    out[:] = inv2N * np.sin(pi_over_Delta*(u-ui))/np.sin(inv2N*pi_over_Delta*(u-ui))

@guvectorize([(float64[:], float64, float64, float64, float64[:])], '(n),(),(),()->(n)')
def lif_symm(ui, u, pi_over_Delta, inv2N, out):
    """Lagrange interpolation function on reduced mesh, for a symmetry quantity.
    Used by `LagrangeMesh._interpolateND_new()`

    Args:
        as for `lif`.
    """
    out[:] = inv2N * ( np.sin(pi_over_Delta * (u - ui)) / np.sin(inv2N * pi_over_Delta * (u - ui)) +
                       np.sin(pi_over_Delta * (u + ui)) / np.sin(inv2N * pi_over_Delta * (u + ui))
                     )


@guvectorize([(float64[:], float64, float64, float64, float64[:])], '(n),(),(),()->(n)')
def lif_skew(ui, u, pi_over_Delta, inv2N, out):
    """Lagrange interpolation function on reduced mesh, for skew-symmetry quantity.
    Used by `LagrangeMesh._interpolateND_new()`

    Args:
        as for `lif`.
    """
    out[:] = inv2N * ( np.sin(pi_over_Delta * (u - ui)) / np.sin(inv2N * pi_over_Delta * (u - ui)) -
                       np.sin(pi_over_Delta * (u + ui)) / np.sin(inv2N * pi_over_Delta * (u + ui))
                     )

