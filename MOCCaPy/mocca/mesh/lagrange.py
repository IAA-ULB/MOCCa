import numpy as np
import numpy.typing as npt
# from mocca.mesh import ijk

from mocca.mesh.observable import Observable
from mocca.mesh.lagrange_function import lagrange_function

from src_heph.heph_symmetries import symmetry


def create_mesh(gx, gy=None, gz=None):
    """Create a rectangular mesh from 1D arrays with coordinates.

    Args:
        gx (ndarray): 1D array
        gy (ndarray): 1D array
        gz (ndarray): 1D array
            gx, gy, gz may refer to the same array.
    Returns:
        an array (order='F)
        if only gx is provided, gx is returned
        if only gx and gy are provided, gxy of shape (nx*ny,2) is returned
        if gx, gy and gz are provided, gxzy of shape (nx*ny*nz,3) is returned
    """
    if gz is not None:
        assert gy is not None, "Argument gy is required"
        nx = len(gx)
        ny = len(gy)
        nz = len(gz)
        gxyz = np.empty((nx*ny*nz,3), dtype=float, order='F')
        gxy = create_mesh(gx, gy)
        nxy = nx*ny
        for iz in range(nz):
            gxyz[iz * nxy:(iz + 1) * nxy, 0:2] = gxy
            gxyz[iz * nxy:(iz + 1) * nxy, 2] = gz[iz]
        return gxyz

    elif gy is not None:
        nx = len(gx)
        ny = len(gy)
        gxy = np.empty((nx*ny,2), dtype=float, order='F')
        for iy in range(ny):
            gxy[iy * nx:(iy + 1) * nx, 0] = gx
            gxy[iy * nx:(iy + 1) * nx, 1] = gy[iy]
        return gxy
    else:
        return gx


def _can_reuse(axes, Q):
    """Return a derivative that can be reused or None."""
    if Q.derivative_is_uptodate(axes):
        return Q.derivatives[axes], _get_symmetry(axes,Q)
    else:
        return None, None

def _get_symmetry(axes, Q):
    """Return the symmetry of the derivative of Q wrt axes components.
    """
    symmetry = Q.symmetry.copy()

    # an even number of differentiations keeps the symmetry sign
    # an odd number of differentiations flips the symmetry sign
    # x-axis
    if Q.mesh.reduced[0]:
        if axes.count('x') % 2:
            symmetry[0,:] *= -1
    # y-axis
    if Q.mesh.dim > 1 and Q.mesh.reduced[1]:
        if axes.count('y') % 2:
            symmetry[1,:] *= -1
    # z-axis
    if Q.mesh.dim > 2 and Q.mesh.reduced[2]:
        if axes.count('z') % 2:
            symmetry[2,:] *= -1

    return symmetry

def _split_axes(axes, Q):
    """split axes in a part that is already computed and a part that still has to be computed.
    Args:
        axes: a sorted string of 'x'|'y'|'z' characters
        Q: the observable being differentiated.
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
        return axis, order, Q.data.reshape((*Q.mesh.shape, Q.n_components), order='F')

    else:
        if nx:
            reused_axes = axes.replace('x', '')
            return 0, nx, Q.derivatives[reused_axes]
            # we only need the symmetry in the x direction, which hasn't changed.

        else:
            # note that this case implies nx==0, ny!=0 and nz!=0.
            reused_axes = axes.replace('y', '')
            return 1, ny, Q.derivatives[reused_axes]


class LagrangeMesh:
    def __init__(self, M: int|tuple,
                       d: int|float|tuple,
                       dim: int=0,
                       bc:str='antiperiodic',
                       reduced:tuple|bool=True,
                       shift:tuple|float=.0,
                       highest_derivative_order:int=2
                 ) -> None:
        """Construct a Lagrange mesh in 1, 2 or 3 dimensions.

        Args:
            dim (int): dimension of mesh. if not specified, dim is guessed as len(M), where M must be a tuple.
            M: tuple with number of grid points in the respective dimensions. Each component of M=2N must be even,
                yielding N grid points at both sides of the origin. If a single M is specified, M is the same in every
                direction.
                The distinction between M and N is consistent with  Ryssens et al, PHYSICAL REVIEW C 92, 064318 (2015)
                section III.B Lagrange-mesh representation.
            d: spacing of points in the respective dimensions. If d is a float or an int d is the same in each direction.
            bc: boundary condition type. 'antiperiodic' or 'periodic'.
            reduced: restrict the mesh to the positive half-axis.
            shift: subtract shift from the grid points. If non-zero, corresponding `reduced` entry must be `False`.

        Raises:
            AssertionError: in case of invalid choices

        Remark:
            The `reduced` parameter is derived from symmetry considerations and may at some point - when the complexity
             of Hephaestos is taken into account - be replaced with a `Symmetry` object. For the time being, however,
             we content with explicitly indicating which coordinate axes must be 'reduced'
        """
        # initialize M and d as dim-tuples
        if dim == 0 and isinstance(M, tuple):
            self.dim = len(M)
        else:
            self.dim = dim

        assert 1 <= self.dim <= 3

        if isinstance(M, int):
            assert M > 0, "M must be strictly positive."
            assert M % 2 == 0, "M must be an even number."
            self.M = tuple(self.dim*[M])
        else:
            assert isinstance(M, tuple)
            assert len(M) == self.dim
            for mi in M:
                assert mi >= 0, "M must be strictly positive."
                assert mi % 2 == 0, "M must be an even number."
            for mi in M:
                assert isinstance(mi, int)
            self.M = M

        # box spacing (ints are converted to floats)
        if isinstance(d,tuple):
            assert len(d) == self.dim
            for di in d:
                assert di > 0, "d must be strictly positive."
            self.d = np.array([float(di) for di in d])
        else:
            assert d > 0, "d must be strictly positive."
            self.d = np.array(self.dim*[float(d)])

        # validate boundary condition
        assert bc in ['antiperiodic', 'periodic']
        self.bc = bc
        # convenience attributes
        self.antiperiodic = bc == 'antiperiodic'
        self.periodic = not self.antiperiodic # since there are only 2 options.

        # validate reduced
        if isinstance(reduced, bool):
            self.reduced = tuple(self.dim*[reduced])
        else:
            assert isinstance(reduced, tuple)
            assert len(reduced) == self.dim
            self.reduced = reduced

        # Compute dv (for integration)
        self.dv = np.prod(self.d) * 2 ** self.reduced.count(True)

        # Compute the unreduced (!) box widths
        # The full (unreduced box width is needed by the plane wave base functions
        self.box_width = np.array([M*d for (M,d) in zip(self.M, self.d)])

        # validate shift
        if isinstance(shift, float):
            self.shift = tuple(self.dim*[shift])
        else:
            assert isinstance(shift, tuple)
            assert len(shift) == self.dim
            self.shift = shift
        for (r,s) in zip(self.reduced, self.shift):
            if r:
                assert s == 0., "A nonzero shift cannot be applied when reduced is True."

        # initialize grid points:
        start = self.dim*[.0]
        n_reduced = self.dim*[0]
        g1D = self.dim*[0]
        for i, (n_i, shift_i, reduced_i, d_i) in enumerate(zip(self.M, self.shift, self.reduced, self.d)):
            if reduced_i:
                n_reduced[i] = n_i//2
                start[i] = 0.5
            else:
                n_reduced[i] = n_i
                start[i] = -(n_i - 1)*0.5 - shift_i/d_i

            g1D[i] = np.linspace(start[i], start[i] + n_reduced[i] - 1, n_reduced[i])
            g1D[i] *= d_i
            # print(f"{i=} {g1D[i]=}")

        self.gx = g1D[0]
        if self.dim > 1:
            self.gy = g1D[1]
        if self.dim > 2:
            self.gz = g1D[2]

        if self.dim == 3:
            self.grid = create_mesh(self.gx, self.gy, self.gz)

            self.gridx = np.empty(n_reduced, order='F')
            self.gridy = np.empty(n_reduced, order='F')
            self.gridz = np.empty(n_reduced, order='F')
            # for k in range(n_reduced[2]):
            #     for j in range(n_reduced[1]):
            #         for i in range(n_reduced[0]):
            #             self.gridx[i, j, k] = g1D[0][i]
            #             self.gridy[i, j, k] = g1D[1][j]
            #             self.gridz[i, j, k] = g1D[2][k]
            for k in range(n_reduced[2]):
                for j in range(n_reduced[1]):
                    self.gridx[:,j,k] = g1D[0]
            for k in range(n_reduced[2]):
                for i in range(n_reduced[0]):
                    self.gridy[i,:,k] = g1D[1]
            for j in range(n_reduced[1]):
                for i in range(n_reduced[0]):
                    self.gridz[i,j,:] = g1D[2]
            # print(f"{self.gridx=}")
            # print(f"{self.gridy=}")
            # print(f"{self.gridz=}")
            # print(f"{self.gridx.ravel(order='F')=}")
            # print(f"{self.gridy.ravel(order='F')=}")
            # print(f"{self.gridz.ravel(order='F')=}")

        elif self.dim == 2:
            self.grid = create_mesh(self.gx, self.gy)

            self.gridx = np.empty(n_reduced, order='F')
            self.gridy = np.empty(n_reduced, order='F')
            for j in range(n_reduced[1]):
                self.gridx[:, j] = g1D[0]
            for i in range(n_reduced[0]):
                self.gridy[i, :] = g1D[1]
            # print(f"{self.gridx=}")
            # print(f"{self.gridy=}")
            # print(f"{self.gridx.ravel(order='F')=}")
            # print(f"{self.gridy.ravel(order='F')=}")

        else: # self.dim == 1
            self.grid = create_mesh(self.gx)
            self.gridx = g1D[0]
            # print(f"{self.gridx=}")

        self.shape = self.gridx.shape
        self.flat_shape = (int(np.prod(self.shape)),)

        self._setup_D_matrices(highest_derivative_order)

    def flatten(self, array=None):
        """Flatten the mesh array.

        Args:
            array: if None, self.gridx|y|z are flattened. Otherwise, array is flattened.

        Returns:
            None or the flattened array.

        Raises:
            AssertionError: if array.shape != self.shape.
        """
        if array is None:
            # self._reshape(flat=True)
            if len(self.gridx.shape) == 1:
                # gridx/gridy/gridz already flat
                # (We assume that either all or none of them are flat or unflattened at the same time)
                return None
            else:
                try:
                    self.gridx = self.flatten(self.gridx)
                    # AttributeError raised if one of gridy/z does not exist
                    self.gridy = self.flatten(self.gridy)
                    self.gridz = self.flatten(self.gridz)
                except AttributeError:
                    pass

        else:
            if array.shape[0] == self.flat_shape[0]:
                # already flat
                return array
            else:
                # collapse the first dim dimensions into one and copy the remaining dimensions:
                flat_shape = list(self.flat_shape) # a single element
                for d in array.shape[self.dim:]:
                    flat_shape.append(d)
                flat_shape = tuple(flat_shape)
                return array.reshape(flat_shape, order='F')

    def unflatten(self, array=None) -> None:
        """Unflatten the mesh array. For most computations

        Args:
            array: if None, self.gridx|y|z are unflattened. Otherwise, array is unflattened.

        Returns:
            None or the unflattened array.

        Raises:
            AssertionError: if array.shape != self.flat_shape.
        """
        if array is None:
            if len(self.gridx.shape) == self.dim:
                # gridx/gridy/gridz already unflattened
                # (We assume that either all or none of them are flat or unflattened at the same time)
                return None
            else:
                try:
                    self.gridx = self.unflatten(self.gridx)
                    # AttributeError raised if one of gridy/z does not exist
                    self.gridy = self.unflatten(self.gridy)
                    self.gridz = self.unflatten(self.gridz)
                except AttributeError:
                    pass
        else:
            if array.shape[0:self.dim] == self.shape[0]:
                # already unflattened
                return array
            else:
                # Take the unflattened shape and append the multicomponent dimensions of array
                shape = list(self.shape)
                for d in array.shape[1:]:
                    shape.append(d)
                shape = tuple(shape)
                return array.reshape(shape, order='F')

    def is_flat(self, Q: Observable) -> bool:
        """Determine if `q` is flat.
        Args:
            Q: a quantity discretised on the grid.
        """
        return Q.data.shape[0] == self.flat_shape[0]

    def is_unflattened(self, Q: npt.NDArray) -> bool:
        """Determine if `q` is unflattened.
        Args:
            Q: a quantity discretised on the grid.
        """
        return Q.shape[0:self.dim] == self.shape

    def is_grid_quantity(self, Q: npt.NDArray) -> bool:
        """Determine if `Q` is a grid quantity.
        Args:
            Q: a quantity discretised on the grid.
        """
        return self.is_flat(Q) or self.is_unflattened(Q)

    def n_gridpoints(self):
        return int(np.prod(self.shape))

    def get_number_of_components(self, Q: npt.NDArray):
        """compute the number of components in agrid quantity `Q`.
        Args:
            Q: a quantity discretised on the grid.
        """
        # TODO: does this belong here? It probably dates from before the Observable class
        return len(Q) // self.flat_shape[0]

    def get_component(self, Q: npt.NDArray, i: int):
        """Get the i-th component of `Q`."""
        # TODO: does this belong here? It probably dates from before the Observable class
        return Q[i * self.flat_shape[0]:(i + 1) * self.flat_shape[0]]

    def apply(self, function):
        """Apply a function on the mesh, i.e. compute the function value on every grid point.

        Args:
            function: function f(x,y,z) to be applied on the mesh array. If f accepts additional parameters, a closure
                must be defined.

        Returns:
            a flat mesh array (self.flat_shape).

        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        # TODO: use arguments for output variables? as in numpy out=
        self.flatten()
        return function(self.gridx, self.gridy, self.gridz)

    def integrate(self, Q):
        """Compute the integral of a scalar quantity `q` on the mesh.

        Args:
            Q (Observable) : discretised on the grid.

        Returns:
            a scalar:

        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        n_components = Q.n_components
        result = np.zeros((n_components,), dtype=float)
        for ic in range(n_components):
            q = Q[ic]
            result[ic] = q.sum() * self.dv
        return result if n_components > 1 else result[0]

    def plane_wave_1D(self, L:float, k:float, r:np.array):
        """Evaluate the plane wave basis function at with wave vector `k` at position `r`.

        According to Ryssens et al, PHYSICAL REVIEW C 92, 064318 (2015), eq(14)
        According to [eq 5.1] in https://github.com/IAA-nuclear/tantalus_full/blob/MOCCaPy/MOCCaPy/mocca/mesh/lagrange.pdf

        Args:
            k: wave vector, shape is `(self.dim,)`. values must be odd half integers ±1/2, ±3/2, ...
            r: position vector, shape is `(self.dim, nr), nr being the number of evaluation points.
            as_complex: if True, return a complex array. otherwise return two arrays with the real and imaginary parts,
                resp.

        Returns:
            array of shape `(nr,2)` with complex and imaginary parts
        """
        arg = (2*np.pi * k / L) * r
        result = np.empty((len(r),2), dtype=float, order='F')
        result[:,0] = np.cos(arg)
        result[:,1] = np.sin(arg)
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
            Array of values of the basis function at r.
        """
        assert isinstance(r, np.ndarray)

        if isinstance(ijk_sign,int):
            ijk_sign = self.dim * (ijk_sign,)

        if self.dim == 3:
            assert r.shape[1] == 3
            i = ijk[0]
            j = ijk[1]
            k = ijk[2]
            two_pi_K = np.array([ 2. * np.pi * ( self.gx[i] if not self.reduced[0] else
                                                 self.gx[i] * ijk_sign[0]
                                               ) / (self.box_width[0] * self.d[0])
                                , 2. * np.pi * ( self.gy[j] if not self.reduced[0] else
                                                 self.gy[j] * ijk_sign[1]
                                               ) / (self.box_width[1] * self.d[1])
                                , 2. * np.pi * ( self.gz[k] if not self.reduced[0] else
                                                 self.gz[k] * ijk_sign[2]
                                               ) / (self.box_width[2] * self.d[2])
                                ], dtype=float, order='F')
            factor = np.sqrt(1/(self.box_width[0] * self.box_width[1] * self.box_width[2]))

        elif self.dim == 2:
            assert r.shape[1] == 2
            i = ijk[0]
            j = ijk[1]
            two_pi_K = np.array([ 2. * np.pi * ( self.gx[i] if not self.reduced[0] else
                                                 self.gx[i] * ijk_sign[0]
                                               ) / (self.box_width[0] * self.d[0])
                                , 2. * np.pi * ( self.gy[j] if not self.reduced[0] else
                                                 self.gy[j] * ijk_sign[1]
                                               ) / (self.box_width[1] * self.d[1])
                                ], dtype=float, order='F')
            factor = np.sqrt(1/(self.box_width[0] * self.box_width[1]))

        elif self.dim == 1:
            if len(r.shape) == 1:
                r = r.reshape((r.shape[0], 1), order='F')
            assert r.shape[1] == 1
            # ijk == i
            i = ijk
            two_pi_K = np.array([ 2. * np.pi * ( self.gx[i] if not self.reduced[0] else
                                                self.gx[i] * ijk_sign[0]
                                               ) / (self.box_width[0] * self.d[0])
                                ], dtype=float, order='F')
            factor = np.sqrt(1 / self.box_width[0])

        arg = r @ two_pi_K
        result = np.empty((r.shape[0],2), dtype=float, order='F')
        np.cos(arg, out=result[:,0])
        np.sin(arg, out=result[:,1])
        result *= factor

        return result

    # Differentiation methods
    # ---------------------------------------------------------------------------
    def _compute_D1(self, axis) -> None:
        """Compute the full matrix D1 (see eq 10.1 in lagrange.md) for axis
        (wether it is reduced or not).

        Returns:
            The full matrix D1 (see eq 10.1 in lagrange.md).
        """
        # We first compute the full D1, also for the reduced case. Then we comput $E^{ll}$ and $E^{lr}$
        twoN = self.M[axis]
        D1 = np.empty((twoN,twoN), dtype=float) # deliberately not order='F' for performance reasons
        alternating_sign = np.empty((twoN+1), dtype=float)
        alternating_sign[::2] = 1
        alternating_sign[1::2] = -1
        # alternating_sign = [1, -1, 1, -1, ...]
        # In agreement with the Warning following equations 10.1 and 10.2 in lagrange.md, which demands
        # an extra sign flip because they use (i-j) in the sine argements instead of (j-i)
        # i = np.linspace(0, twoN - 1, twoN)
        minus_i = -np.linspace(0, twoN - 1, twoN)
        # In agreement with the formulas in Ryssens et al 2015 (eq 18),
        # j is the row index and i the column index (D_ji)
        for j in range(twoN):
            # row = i-j
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
        return self.D.shape[1]

    def _setup_D_matrices(self, highest_derivative_order):
        """setup D matrices infrastructure"""
        self.D = np.empty((3, highest_derivative_order), dtype=np.ndarray)
        #   stores the full D matrix if the axis is not reduced,
        #   and the reduced matrix sum D+E, otherwise
        self.DmE = np.empty((3, highest_derivative_order), dtype=np.ndarray)
        #   stores nothing if the axis is not reduced,
        #   and the reduced matrix difference D-E otherwise.
        self.D[0, 0] = self._compute_D1(0)
        if self.dim > 1:
            if (self.M[1] == self.M[0]) and \
                    (self.d[1] == self.d[0]) and \
                    (self.reduced[1] == self.reduced[0]):
                self.D[1, 0] = self.D[0, 0]
            else:
                self.D[1, 0] = self._compute_D1(1)

            if self.dim > 2:
                if (self.M[2] == self.M[0]) and \
                        (self.d[2] == self.d[0]) and \
                        (self.reduced[2] == self.reduced[0]):
                    self.D[2, 0] = self.D[0, 0]
                elif (self.M[2] == self.M[1]) and \
                        (self.d[2] == self.d[1]) and \
                        (self.reduced[2] == self.reduced[1]):
                    self.D[2, 0] = self.D[1, 0]
                else:
                    self.D[2, 0] = self._compute_D1(2)

        for axis in range(self.dim):
            for order in range(1, highest_derivative_order):
                self.D[axis, order] = self.D[axis, 0] * self.D[axis, order - 1]

            if self.reduced[axis]:
                N = self.M[axis]//2
                for order in range(1, highest_derivative_order):
                    # Storing D
                    E = np.empty((N, N), dtype=float)
                    D = np.empty((N, N), dtype=float)

                    # Reverse the columns of the lower left quadrant Dll and copy into ED
                    # Copy the first column of Dll

                    for icol in range(N):
                        E[:,N-1-icol] = self.D[axis, order][N:,icol] # reverse the order of Dll (lower left quadrant)
                        D[:,icol] = self.D[axis, order][N:,N+icol]   # copy Dlr (lower right quadrant)
                    self.D[axis, order] = D + E # replaces the full D matrix
                    self.DmE[axis, order] = D - E

    def _get_D(self, axis:int, order:int):
        return self.D[axis, order - 1]

    def _get_DpE(self, axis:int, order:int):
        return self.D[axis, order - 1]

    def _get_DmE(self, axis:int, order:int):
        return self.DmE[axis, order - 1]

    def differentiate(self, Q, axes, accumulate_in=None):
        """Differentiate Q wrt axes.
        Args:
            Q: Observable to be differentiated.
            axes: Axes to be differentiated.
                A single derivative can be requested as a `str` combining the characters 'x', 'y', 'z'. E.g. the
                 `str` 'xyz' requests d^3/dxdydz. Multiple derivatives can be requested too:
                   - 'Grad': all first order derivatives
                   - 'Hessian': all second order derivatives
                   - 'Laplacian' d^2/dx^2 + d^2/dy^2 + d^2/dz^2)
                   - 'Tensor3': all third order derivatives
                   - 'Tensor4': all fourth order derivatives
                 In principle the order of the 'xyz' characters in axes is immaterial. Axes is sorted and repeating
                 characters are combined. E.g. 'xyx' -> 'xxy' is computed as D2x x D1y x Q.
                 Furthermore, axes can be a `list` of such strings, requesting the computation of many different
                 derivatives. The list is sorted from low to high differentiation order, and higher order derivatives
                 will reuse previously computed derivatives. E.g., if `xx` has been computed and
        Returns:
            If a single derivative is requested, returns a `ndarray` containing the requested derivative. (The Laplacian
            counts as a single derivative too). Otherwise a tensor is returned containing `ndarray`s with the individual
            derivatives.
        """
        if isinstance(axes, str):
            if axes[0].isupper():
                # All multi-component derivatives, i.e. those requiring at least two simple derivatives
                assert self.dim >= 2
                if axes == 'Grad':
                    result = np.empty((self.dim,), dtype=np.ndarray)
                    result[0] = self.differentiate(Q, axes='x')
                    result[1] = self.differentiate(Q, axes='y')
                    if self.dim == 3:
                        result[2] = self.differentiate(Q, axes='z')

                    return result

                elif axes == 'Hessian':
                    result = np.empty((self.dim,self.dim), dtype=np.ndarray)
                    if self.dim == 2:
                        self.differentiate(Q, axes=['x', 'y'], accumulate_in=result)
                    else:
                        self.differentiate(Q, axes=['x', 'y', 'z'], accumulate_in=result)
                    result[0, 0] = self.differentiate(Q, axes='xx')
                    result[0, 1] = self.differentiate(Q, axes='xy')
                    result[1, 0] = self.differentiate(Q, axes='yx')
                    result[1, 1] = self.differentiate(Q, axes='yy')
                    if self.dim == 3:
                        result[0, 2] = self.differentiate(Q, axes='xz')
                        result[1, 2] = self.differentiate(Q, axes='yz')
                        result[2, 0] = self.differentiate(Q, axes='zx')
                        result[2, 1] = self.differentiate(Q, axes='zy')
                        result[2, 2] = self.differentiate(Q, axes='zz')

                    return result

                elif axes == 'Laplacian':
                    assert self.dim >= 2
                    if axes in Q.derivatives:
                        out = Q.derivatives[axes]
                    else:
                        out = np.zeros_like(self.data)
                        Q.derivatives[axes] = out

                    self.differentiate(Q, axes='xx', accumulate_in=out)
                    self.differentiate(Q, axes='yy', accumulate_in=out)
                    if self.dim >= 2:
                        self.differentiate(Q, axes='zz', accumulate_in=out)

                    return out

                elif axes == 'Tensor3':
                    raise NotImplementedError(f"{axes=} is not implemented.")

                if axes == 'Tensor4':
                    raise NotImplementedError(f"{axes=} is not implemented.")

                else:
                    raise NotImplementedError(f"{axes=} is not implemented.")

            else:
                # TODO: The D matrix products for differentiation require Q in grid shape, instead of flattened.
                #       The resulting derivatives are also in grid shape. Devise a mechanism to keep track of the shape.

                # All simple derivatives. `axes` is composed as a sequence of 'x'|'y'|'z' characters.
                # sort the `axes` str, as the order of differentiation is immaterial
                axes = ''.join(sorted(axes)) # E.g. 'xyzx' -> 'xxyz', which is  evaluated as Dx2*Dy*Dz*Q
                if Q.derivative_is_uptodate(axes):
                    return Q.derivatives[axes]

                # Find out if we can reuse a previously computed derivative as a starting point. E.g.:
                # - axes = 'xyz', if 'z' is already computed, use it as a starting point and apply D1x*D1y to dQdz
                # - axes = 'xyz', if 'y' is already computed, use it as a starting point and apply D1x*D1z to dQdy
                # - axes = 'xyz', if 'yz' is already computed, use it as a starting point and apply D1x to d2Qdydz
                # - axes = 'xx',  if 'x' is already computed, do NOT it as a starting point and apply D2x  to Q
                # We prefer to reuse component with 'y' and 'z', because the einsum operations imply non-contiguous

                # memory access
                # Reuse previously allocated memory, or allocate
                if axes in Q.derivatives:
                    # Reuse previously allocated memory
                    out = Q.derivatives[axes]
                else:
                    # allocate memory for the derivative.
                    out = np.empty((*self.shape,Q.n_components), dtype=float, order='F')
                    Q.derivatives[axes] = out

                axis, order, op2 = _split_axes(axes, Q)

                if axis == 0: # x-axis
                    subscripts = 'il,lq'  if (self.dim == 1) else \
                                 'il,ljq' if (self.dim == 2) else \
                                 'il,ljkq'
                    if not self.reduced[axis]:
                        D = self._get_D(axis=0, order=order)
                        np.einsum(subscripts, D, op2, out=out)
                    else:
                        for iq in range(Q.n_components):
                            DE = self._get_DpE(axis=0, order=order) if symmetry[0, iq] == 1 else \
                                 self._get_DmE(axis=0, order=order)
                            if self.dim == 2:
                                np.einsum(subscripts, DE, op2[:, :, :, iq], out=out)
                            elif self.dim == 1:
                                np.einsum(subscripts, DE, op2[:, :, iq], out=out)
                            else:
                                np.einsum(subscripts, DE, op2[:, iq], out=out)
                            if nx % 2:
                                symmetry[0, iq] *= -1

                elif axis == 1: # y-axis
                    subscripts = 'il,jlq' if (self.dim == 2) else \
                                 'il,jlkq'
                    if not self.reduced[axis]:
                        D = self._get_D(axis=1, order=order)
                        np.einsum(subscripts, D, op2, out=out)
                    else:
                        for iq in range(Q.n_components):
                            DE = self._get_DpE(axis=1, order=order) if symmetry[1, iq] == 1 else \
                                 self._get_DmE(axis=1, order=order)
                            if self.dim == 2:
                                np.einsum(subscripts, DE, op2[:, :, :, iq], out=out)
                            else:
                                np.einsum(subscripts, DE, op2[:, :, iq], out=out)
                            if ny % 2:
                                symmetry[1, iq] *= -1

                else: # z-axis
                    subscripts = 'il,jklq'
                    if not self.reduced[axis]:
                        D = self._get_D(axis=2, order=order)
                        np.einsum(subscripts, D, op2, out=out)
                    else:
                        for iq in range(Q.n_components):
                            DE = self._get_DpE(axis=2, order=order) if symmetry[2,iq] == 1 else \
                                 self._get_DmE(axis=2, order=order)
                            np.einsum(subscripts, DE, op2[:,:,:,iq], out=out)
                            if nz % 2:
                                symmetry[2, iq] *= -1

            Q.derivative_set_uptodate(axes)
            return out

        elif isinstance(axes, list):
            # Add components to allow reuse of derivatives:
            # every compose
            if 'Grad' in axes:
                axes = ['x', 'y', 'z'] + axes

            if 'Laplacian' in axes:
                axes = ['xx', 'yy', 'zz'] + axes

            if 'Hessian' in axes:
                axes = ['y', 'z',
                        'xx', 'xy', 'xz',
                              'yy', 'yz',
                                    'zz'
                       ] + axes
                # As D2x*Q is computed equally efficient as D1x*dQdx and the the same holds for y and z,
                # there is no point in reusing 'x' for 'xx'. However, reusing 'y' for 'xy' and 'z' for 'xz'
                # and 'yz' is beneficial.
                # The genaral rule is that everything starting with 'x' can be dropped, except the one needed
                # by the composite, c.q 'Hessian'.

            if 'Tensor3' in axes:
                axes = ['y', 'z',
                        'yy', 'yz', 'zz',
                        'xxx', 'xxy', 'xxz',
                               'xyy', 'xyz',
                                      'xzz',
                        'yyy', 'yyz',
                               'yzz',
                        'zzz',
                       ] + axes
                # 'x', 'xx' and and 'xy' are dropped for the same reason as above.

            if 'Tensor4' in axes:
                axes = ['y', 'z',
                        'yy', 'xz', 'yz', 'zz',
                        'yyy', 'yyz', 'yzz', 'zzz',
                        'xxxx', 'xxxy', 'xxxz',
                                'xxyy', 'xxyz',
                                        'xxzz',
                                'xyyy', 'xyyz',
                                        'xyzz',
                                        'xzzz',
                                'yyyy', 'yyyz',
                                        'yyzz',
                                        'yzzz',
                                        'zzzz',
                        ] + axes
                # all starting with 'x' and order < 4 are dropped.

            axes = list(set(axes)) # now every entry occurs only once
            # sort the list in-place:
            # - composite derivatives come first, simple derivatives in alphabetical order and 'x' < 'xx'
            axes.sort()
            # move the composite derivatives (starting with a capital) to the back.
            while axes[0][0].isupper():
                composite = axes.pop(0)
                axes.append(composite)

            for ax in axes:
                self.differentiate(Q, axes=ax)

            return None  # returning a list would make no sense, the user must access the requested derivatives via
                         # `Q.derivatives`

        raise ValueError(f"Axes {axes} not supported.")

    # Interpolation methods
    #---------------------------------------------------------------------------
    def interpolate(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate a quantity `Q` on the mesh.

        Args:
            Q: representation of a scalar quantity on the grid.
            r: array of `p` points at which to interpolate the scalar quantity `q`. 'r.shape == (p, self.dim)'
        """
        # TODO (?) speed this stuff up... Performance may be wrecked by numerous nested python loops.
        self.flatten()
        
        if self.dim == 3:
            return self._interpolate3D(Q, r)
        elif self.dim == 2:
            return self._interpolate2D(Q, r)
        else:
            return self._interpolate1D(Q, r)


    def _interpolate1D(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 1D  mesh."""
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.zeros(shape=(nr,nq), dtype=float, order='F')
        for ig in range(self.n_gridpoints()):
            for iq in range(nq):
                q = Q[iq]
                sign = Q.sign(iq)
                fx = self.lagrange_function(r, ig, sign)
                Qr[:,iq] += q[ig] * fx

        return Qr


    def _interpolate2D(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 2D mesh.

        Args:
            Q: An observable with values specified on all grid points
            r: array of points at which to interpolate Q. 'r.shape == (nr, self.dim)'
        Returns:
            An array with the interpolated values.
        """
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.zeros(shape=(nr,nq), dtype=float, order='F')
        ig = 0
        Nx = len(self.gx)
        Ny = len(self.gy)
        for iy in range(Ny):
            for ix in range(Nx):
                for iq in range(nq):
                    q = Q[iq]
                    sign = Q.sign(iq)
                    fxy = self.lagrange_function(r, (ix,iy), sign)
                    Qr[:,iq] += q[ig] * fxy
                ig += 1

        return Qr


    def _interpolate3D(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 2D mesh.

        Args:
            Q: An observable with values specified on all grid points
            r: array of points at which to interpolate Q. 'r.shape == (nr, self.dim)'
        """
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.zeros(shape=(nr, nq), dtype=float, order='F')
        ig = 0
        Nx = len(self.gx)
        Ny = len(self.gy)
        Nz = len(self.gz)
        for iz in range(Nz):
            for iy in range(Ny):
                for ix in range(Nx):
                    for iq in range(nq):
                        q = Q[iq]
                        sign = Q.sign(iq)
                        fxy = self.lagrange_function(r, (ix, iy,iz), sign)
                        Qr[:, iq] += q[ig] * fxy
                    ig += 1
        return Qr


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

        if self.dim == 3:

            i,j,k = ijk
            x_i, y_j, z_k = self.gx[i], self.gy[j], self.gz[k]

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
            x_i, y_j = self.gx[i], self.gy[j]

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

