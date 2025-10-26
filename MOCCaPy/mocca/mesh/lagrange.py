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
    def __init__(self, n: int|tuple, d: int|float|tuple, dim: int=0,
                       bc='antiperiodic',
                       reduced:tuple|bool = True,
                       shift:tuple|float = .0,
                 ) -> None:
        """Construct a Lagrange mesh in 1, 2 or 3 dimensions.

        Args:
            dim (int): dimension of mesh. if not specified, dim is guessed as len(n), where n must be a tuple.
            n: number of points in the respective dimensions. If n is an int n is the same in each direction.
            d: spacing of points in the respective dimensions. If n is a float or an int d is the same in each direction.
            bc: boundary condition type. 'antiperiodic' or 'periodic'.
            reduced: restrict the mesh to the positive half-axis. The corresponding `n` entry is halved.
            shift: subtract shift from the grid points. If non-zero, corresponding `reduced` entry must be `False`.

        Raises:
            AssertionError: in case of invalid choices

        Remark:
            The `reduced` parameter is derived from symmetry considerations and may at some point - when the complexity
             of Hephaestos is taken into account - be replaced with a `Symmetry` object. For the time being, however,
             we content with explicitly indicating which coordinate axes must be 'reduced'
        """
        # initialize n and d as dim-tuples
        if dim == 0 and isinstance(n, tuple):
            self.dim = len(n)
        else:
            self.dim = dim

        assert 1 <= self.dim <= 3

        if isinstance(n, int):
            assert n > 0, "n must be strictly positive."
            assert n % 2 == 0, "n must be an even number."
            self.n = tuple(self.dim*[n])
        else:
            assert isinstance(n, tuple)
            assert len(n) == self.dim
            for ni in n:
                assert ni >= 0, "n must be strictly positive."
                assert ni % 2 == 0, "n must be an even number."
            for ni in n:
                assert isinstance(ni, int)
            self.n = n

        # box spacing (ints are converted to floats)
        if isinstance(d,tuple):
            assert len(d) == self.dim
            for di in d:
                assert di > 0, "d must be strictly positive."
            self.d = tuple(float(di) for di in d)
        else:
            assert d > 0, "d must be strictly positive."
            self.d = tuple(self.dim*[float(d)])

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
        self.box_width = np.array([n*d for (n,d) in zip(self.n, self.d)])

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
        for i, (n_i, shift_i, reduced_i, d_i) in enumerate(zip(self.n, self.shift, self.reduced, self.d)):
            if reduced_i:
                n_reduced[i] = n_i//2
                start[i] = 0.5
            else:
                n_reduced[i] = n_i
                start[i] = -(n_i - 1)*0.5 - shift_i/d_i

            g1D[i] = np.linspace(start[i], start[i] + n_reduced[i] - 1, n_reduced[i])
            g1D[i] *= d_i
            # print(f"{i=} {g1D[i]=}")

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


    def plane_wave(self, k:np.array, r:np.array):
        """Evaluate the plane wave basis function at with wave vector `k` at position `r`.

        Args:
            k: wave vector, shape is `(self.dim,)`. values must be odd half integers +/-1/2, +/-3/2, ...
            r: position vector, shape is `(self.dim, nr), nr being the number of evaluation points.
            as_complex: if True, return a complex array. otherwise return two arrays with the real and imaginary parts,
                resp.

        Returns:
            array of shape `(2,nr)` with complex and imaginary parts

        Remarks:
            It would be nice to return a matrix of shape `(nr,2)` with the real and imaginary parts as contiguous
            columns. However, this require a copy and moving the data. As this is inherently inefficient we do not
            facilitate this. If necessary, we could delegate this to Fortran code.
        """

        one_over_sqrt_bw = np.sqrt(1 / np.prod(self.box_width))
        two_pi_j = np.pi * 2j
        k = np.array([0.5,1.5,2.5])
        k /= self.box_width
        pw = one_over_sqrt_bw*np.exp(two_pi_j*r@k)
        return pw


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
            pass # TODO implement
        elif self.dim == 2:
            pass # TODO implement
        else:
            return self._interpolate1D(Q, r)

        # ijk = IJK(self.shape)
        # n_gridpoints = ijk.ig_end
        # n_components = self.get_number_of_components(Q) # TODO n_components should be an attribute of Q, rather than self
        # i,j,k = ijk.index
        # while not ijk.done():  # loop over gridpoints
        #     for iq in range(n_components):
        #         q = self.get_component(Q, iq) # TODO should be method of Q, rather than self
        #         # loop over interpolation points r
        #     ijk.inc()


    def _interpolate1D(self, Q:Observable, r:npt.NDArray) -> npt.NDArray:
        """Interpolate `Q` on a 1D  mesh."""
        nr = r.shape[0]
        nc = Q.n_components
        Qr = np.empty(shape=(nr,nc), dtype=float, order='F')
        N = len(self.gridx) if self.reduced[0] else len(self.gridx)//2
        for i in range(self.n_gridpoints()):
            for j in range(nc):
                q = Q[j]
                sign = Q.symmetry[j]
                if self.reduced[0]:
                    Qr[:,j] = q[i] * (        lagrange_function(r,  self.gridx[i], self.d[0], N)
                                     + sign * lagrange_function(r, -self.gridx[i], self.d[0], N)
                                     )
                else:
                    Qr[:,j] = q[i] * lagrange_function(r, self.gridx[i], self.d[0], N)
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
