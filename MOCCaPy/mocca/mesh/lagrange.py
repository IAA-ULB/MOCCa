import numpy as np
import numpy.typing as npt
from mocca.mesh.observable import Observable
from numpy.ma.core import shape


def lagrange_function(x:npt.NDArray|float, x_i:float, d:float, N:int):
    """Compute the 1D Lagrange function on `x`.

    Args:
        x: array of points at which to compute the Lagrange function.
        x_i: a grid point : ±1/2 dx, ±3/2 dx, ±5/2 dx, ...
        d: grid spacing
        N: number of grid points on the positive axis
        reduced: whether the axis is reduced or not.
    """
    one_over_2N = 1/(2*N)
    A_i = (np.pi/d)*(x - x_i)
    result = one_over_2N*np.sin(A_i)/np.sin(one_over_2N*A_i)
    # if reduced:
    #     A_i = (np.pi / d) * (x + x_i)
    #     result += one_over_2N * np.sin(A_i) / np.sin(one_over_2N * A_i)
    return result

class LagrangeMesh:
    def __init__(self, N: int|tuple, d: int|float|tuple, dim: int=0,
                       bc='antiperiodic',
                       reduced:tuple|bool = True,
                       shift:tuple|float = .0,
                 ) -> None:
        """Construct a Lagrange mesh in 1, 2 or 3 dimensions.

        Args:
            dim (int): dimension of mesh. if not specified, dim is guessed as len(N), where N must be a tuple.
            N: number of points in the respective dimensions. If N is an int N is the same in each direction.
            d: spacing of points in the respective dimensions. If N is a float or an int d is the same in each direction.
            bc: boundary condition type. 'antiperiodic' or 'periodic'.
            reduced: restrict the mesh to the positive half-axis. The corresponding `N` entry is halved.
            shift: subtract shift from the grid points. If non-zero, corresponding `reduced` entry must be `False`.

        Raises:
            AssertionError: in case of invalid choices

        Remark:
            The `reduced` parameter is derived from symmetry considerations and may at some point - when the complexity
             of Hephaestos is taken into account - be replaced with a `Symmetry` object. For the time being, however,
             we content with explicitly indicating which coordinate axes must be 'reduced'
        """
        # initialize N and d as dim-tuples
        if dim == 0 and isinstance(N, tuple):
            self.dim = len(N)
        else:
            self.dim = dim

        assert 1 <= self.dim <= 3

        if isinstance(N, int):
            assert N > 0, "N must be strictly positive."
            assert N % 2 == 0, "N must be an even number."
            self.N = tuple(self.dim*[N])
        else:
            assert isinstance(N, tuple)
            assert len(N) == self.dim
            for ni in N:
                assert ni >= 0, "N must be strictly positive."
                assert ni % 2 == 0, "N must be an even number."
            for ni in N:
                assert isinstance(ni, int)
            self.N = N

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
        self.box_width = np.array([N*d for (N,d) in zip(self.N, self.d)])

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
        for i, (n_i, shift_i, reduced_i, d_i) in enumerate(zip(self.N, self.shift, self.reduced, self.d)):
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
            self.gridx = g1D[0]
            # print(f"{self.gridx=}")

        self.shape = self.gridx.shape
        self.flat_shape = (int(np.prod(self.shape)),)

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
        # TODO use arguments for output variables?
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


    def basis_function(self, ijk, r):
        """Evaluate the plane wave basis function corresponding to grid point `x_ijk = [x_i,y_j,z_k]` 
        at position `r`. The wave vector `k` is related to `x_i` as `k = x_i/dx`.

        According to [eq 5.1] in https://github.com/IAA-nuclear/tantalus_full/blob/MOCCaPy/MOCCaPy/mocca/mesh/lagrange.pdf

        Args:
            ijk: grid point index. A D-tuple of grid indices with `D == self.dim'.  If a grid axis
                is reduced, grid indices run from 0, ±1, ±2, ..., ±N. Negative values select basis f
                unctions corresponding to grid points on the negative axis. Otherwise, grid indices
                run from 0 to `2N-1`.
        """
        if self.dim == 1:
            # ijk == i
            i = ijk
            two_pi_K = 2. * np.pi * ( self.gx[i] if not self.reduced[0] else
                                      self.gx[i] if i>0 else
                                     -self.gx[i] ) / (self.box_width[0] * self.d[0])
            arg = r * two_pi_K
            factor = np.sqrt(1/self.box_width[0])
        elif self.dim == 2:
            i = ijk[0]
            j = ijk[1]
            two_pi_K = np.array([ 2. * np.pi * ( self.gx[i] if not self.reduced[0] else
                                                 self.gx[i] if i>0 else
                                                -self.gx[i] ) / (self.box_width[0] * self.d[0])
                                , 2. * np.pi * ( self.gy[j] if not self.reduced[0] else
                                                 self.gy[j] if i>0 else
                                                -self.gy[j] ) / (self.box_width[1] * self.d[1])
                                ])
            arg = r @ two_pi_K
            factor = np.sqrt(1/(self.box_width[0] * self.box_width[1]))

        else:
            i = ijk[0]
            j = ijk[1]
            k = ijk[2]
            two_pi_K = np.array([ 2. * np.pi * ( self.gx[i] if not self.reduced[0] else
                                                 self.gx[i] if i>0 else
                                                -self.gx[i] ) / (self.box_width[0] * self.d[0])
                                , 2. * np.pi * ( self.gy[j] if not self.reduced[0] else
                                                 self.gy[j] if i>0 else
                                                -self.gy[j] ) / (self.box_width[1] * self.d[1])
                                , 2. * np.pi * ( self.gz[k] if not self.reduced[0] else
                                                 self.gz[k] if i>0 else
                                                -self.gz[k] ) / (self.box_width[2] * self.d[2])
                                ])
            arg = r @ two_pi_K
            factor = np.sqrt(1/(self.box_width[0] * self.box_width[1] * self.box_width[2]))

        nr = r.shape[0] 
        result = np.empty((nr,2), dtype=float, order='F')
        result[:,0] = np.cos(arg)
        result[:,1] = np.sin(arg)
        result *= factor
        return result


    def derive1(self, Q:npt.NDArray):
        """Compute the 1st order derivative of a scalar quantity `q` on the mesh.

        Args:
            q: scalar quantity discretised on the grid. Thus `q.shape in [self.shape, self.flat_shape]` evaluates
                to True

        Returns:


        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        # TODO : implement

    def derive2(self, Q:npt.NDArray):
        """Compute the 2nd order derivative of a scalar quantity `q` on the mesh.

        Args:
            q: scalar quantity discretised on the grid. Thus `q.shape in [self.shape, self.flat_shape]` evaluates
                to True

        Returns:


        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        # TODO : implement
        
    
    def is_flat(self, Q:Observable) -> bool:
        """Determine if `q` is flat.
        Args:
            Q: a quantity discretised on the grid.
        """
        return Q.data.shape[0] == self.flat_shape[0]
    
    
    def is_unflattened(self, Q:npt.NDArray) -> bool:
        """Determine if `q` is unflattened.
        Args:
            Q: a quantity discretised on the grid.
        """
        return Q.shape[0:self.dim] == self.shape        


    def is_grid_quantity(self, Q:npt.NDArray) -> bool:
        """Determine if `Q` is a grid quantity.
        Args:
            Q: a quantity discretised on the grid.
        """
        return self.is_flat(Q) or self.is_unflattened(Q)


    def n_gridpoints(self):
        return np.prod(self.shape)
    
    def get_number_of_components(self, Q:npt.NDArray):
        """compute the number of components in agrid quantity `Q`.
        Args:
            Q: a quantity discretised on the grid.
        """
        return len(Q)//self.flat_shape[0]
    
    
    def get_component(self, Q:npt.NDArray, i:int):
        """Get the i-th component of `Q`."""
        return Q[i*self.flat_shape[0]:(i+1)*self.flat_shape[0]]


    def interpolate(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate a quantity `Q` on the mesh.

        Args:
            Q: representation of a scalar quantity on the grid.
            r: array of `p` points at which to interpolate the scalar quantity `q`. 'r.shape == (p, self.dim)'
        """
        self.flatten()
        self.flatten(Q.data)
        assert self.is_flat(Q)

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
        Qr = np.empty(shape=(nr,nq), dtype=float, order='F')
        N = len(self.gridx) if self.reduced[0] else len(self.gridx)//2
        for ig in range(self.n_gridpoints()):
            for iq in range(nq):
                q = Q[iq]
                sign = Q.symmetry[iq]
                if self.reduced[0]:
                    Qr[:,iq] = q[ig] * (        lagrange_function(r,  self.gridx[ig], self.d[0], N)
                                       + sign * lagrange_function(r, -self.gridx[ig], self.d[0], N)
                                       )
                else:
                    Qr[:,iq] = q[ig] * lagrange_function(r, self.gridx[ig], self.d[0], N)
        return Qr


    def _interpolate2D(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 2D  mesh."""
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.empty(shape=(nr,nq), dtype=float, order='F')
        N = len(self.gridx) if self.reduced[0] else len(self.gridx)//2
        for ig in range(self.n_gridpoints()):
            for iq in range(nq):
                q = Q[iq]
                sign = Q.symmetry[iq]
                if self.reduced[0]:
                    Qr[:,iq] = q[ig] * (        lagrange_function(r,  self.gridx[ig], self.d[0], N)
                                       + sign * lagrange_function(r, -self.gridx[ig], self.d[0], N)
                                       ) \
                                     * (        lagrange_function(r,  self.gridy[ig], self.d[1], N)
                                       + sign * lagrange_function(r, -self.gridy[ig], self.d[1], N)
                                       )
                else:
                    Qr[:,iq] = q[ig] * lagrange_function(r, self.gridx[ig], self.d[0], N) \
                                     * lagrange_function(r, self.gridy[ig], self.d[1], N)
        return Qr


    def _interpolate3D(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 2D  mesh."""
        nr = r.shape[0]
        nq = Q.n_components
        Qr = np.empty(shape=(nr,nq), dtype=float, order='F')
        N = len(self.gridx) if self.reduced[0] else len(self.gridx)//2
        for ig in range(self.n_gridpoints()):
            for iq in range(nq):
                q = Q[iq]
                sign = Q.symmetry[iq]
                if self.reduced[0]:
                    Qr[:,iq] = q[ig] * (        lagrange_function(r,  self.gridx[ig], self.d[0], N)
                                       + sign * lagrange_function(r, -self.gridx[ig], self.d[0], N)
                                       ) \
                                     * (        lagrange_function(r,  self.gridy[ig], self.d[1], N)
                                       + sign * lagrange_function(r, -self.gridy[ig], self.d[1], N)
                                       ) \
                                     * (        lagrange_function(r,  self.gridz[ig], self.d[2], N)
                                       + sign * lagrange_function(r, -self.gridz[ig], self.d[2], N)
                                       )
                else:
                    Qr[:,iq] = q[ig] * lagrange_function(r, self.gridx[ig], self.d[0], N) \
                                     * lagrange_function(r, self.gridy[ig], self.d[1], N) \
                                     * lagrange_function(r, self.gridz[ig], self.d[2], N)
        return Qr


def lagrange_function(x:npt.NDArray, x_i:float, dx:float, N:int):
    """Compute the 1D Lagrange function on `x`.

    Args:
        x: array of points at which to compute the Lagrange function.
        x_i: grid point: ±1/2, ±3/2, ±5/2, ...
        dx: grid spacing
        N: number of grid points on the positive axis
    Returns:
        An array of values of Lagrange function corresponding to `x`.
    Caveat:
        The function returns nan when `x == x_i` because it evaluates 0/0. Checking for this corner case implies giving
        up numpy array functions. As we assume that the user is not interested in finding the value at x_i because it
        is known to be 1. The user must avoid this. However, for `x = x_i + 1e-9` the result is very close to 1,
    """
    one_over_2N = 1/(2*N)
    A_i = np.pi/dx * (x - x_i)
    return one_over_2N*np.sin(A_i)/np.sin(one_over_2N*A_i)
