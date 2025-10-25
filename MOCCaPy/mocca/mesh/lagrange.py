import numpy as np
import numpy.typing as npt

class IJK:
    """IJK is an index object to iterate over the grid points of a rectangular mesh, like the LagrangeMesh. It
    provides a `ig` attribute running from 1 to n_gridpoints, and attributes `i`, `j` and `k`, referring to the
    spatial indices of the grid points. The `ig` attribute corresponds to flattened arrays, and the `i`, `j` and `k`
    attributes correspond to the unflattened arrays.
    """
    def __init__(self, shape):
        """IJK constructor.
        Args:
            shape (tuple): shape of the mesh, can be 1D, 2D, or 3D.
        """
        self.dim = len(shape)
        # self.shape is a 3D tuple, even if self.dim<3, self.shape[i] == 1 for i>self.dim-1
        self.shape = list(shape)
        self.shape.extend([1,1])
        self.shape = self.shape[0:3]
        self.shape = tuple(self.shape)
        self.ig_end = int(np.prod(np.array(shape)))
        self.reset()


    def inc(self):
        """Increment the index to point to the next grid point in memory.
        Returns:
            The spatial index of the grid point pointed to after being incremented. Three values (`i`, `j` and `k`) are
            returned even if the mesh is 1D or 2D. In the latter case, only `i`, resp. `i` and `j` are relevant.
        """
        self.ig += 1
        self.index[0] += 1
        if self.index[0] == self.shape[0]:
            self.index[0] = 0
            if self.dim > 1:
                self.index[1] += 1
                if self.index[1] == self.shape[1]:
                    self.index[1] = 0
                    if self.dim > 2:
                        self.index[2] += 1
        return self.index

    def done(self) -> bool:
        """Returns `True` if the last `inc` call incremented past the end of the grid points, `False` otherwise."""
        return self.ig >= self.ig_end

    def ijk2ig(self, ijk):
        """ijk->ig"""
        return ijk[0] + ijk[1] * self.shape[0] + ijk[2] * (self.shape[0] * self.shape[1])

    def ig2ijk(self, ig):
        """ig->ijk
        returns:
            [i,j,k] corresponding to ig. Always returns 3 indices, even if `self.dim < 3`. Indices for non-existing
            dimensions are zero. This behavior is consistent with `self.index`.
        """
        ijk =[0,0,0]
        if self.dim > 2:
            ijk[2] = ig // (self.shape[0]*self.shape[1])
            ig %= (self.shape[0]*self.shape[1])
        if self.dim > 1:
            ijk[1] = ig // self.shape[0]
            ig %= self.shape[0]
        ijk[0] = ig
        return ijk


    def reset(self, ig: int = None, ijk=None):
        """Resets the index to point to the first grid point in the mesh."""
        # TODO implement ig and ijk

        if ig is None and ijk is None:
            self.index = [0,0,0]
            self.ig = 0
        elif ig is not None:
            if self.dim == 3:
                self.index[2] = ig//self.shape[0]*self.shape[1]
        elif ijk is not None:
            if self.dim == 3:
                self.ig = ijk[0] + ijk[1] * self.shape[0] + ijk[2] * (self.shape[0] * self.shape[1])


class LagrangeMesh:
    def __init__(self, n: int|tuple, d: int|float|tuple, dim: int=0,
                       bc='antiperiodic',
                       reduce:tuple|bool = True,
                       shift:tuple|float = .0,
                 ) -> None:
        """Construct a Lagrange mesh in 1, 2 or 3 dimensions.

        Args:
            dim (int): dimension of mesh. if not specified, dim is guessed as len(n), where n must be a tuple.
            n: number of points in the respective dimensions. If n is an int n is the same in each direction.
            d: spacing of points in the respective dimensions. If n is a float or an int d is the same in each direction.
            bc: boundary condition type. 'antiperiodic' or 'periodic'.
            reduce: restrict the mesh to the positive half-axis. The corresponding `n` entry is halved.
            shift: subtract shift from the grid points. If non-zero, corresponding `reduce` entry must be `False`.

        Raises:
            AssertionError: in case of invalid choices

        Remark:
            The `reduce` parameter is derived from symmetry considerations and may at some point - when the complexity
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

        # validate reduce
        if isinstance(reduce, bool):
            self.reduce = tuple(self.dim*[reduce])
        else:
            assert isinstance(reduce, tuple)
            assert len(reduce) == self.dim
            self.reduce = reduce

        # Compute dv (for integration)
        self.dv = np.prod(self.d) * 2 ** self.reduce.count(True)

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
        for (r,s) in zip(self.reduce, self.shift):
            if r:
                assert s == 0., "A nonzero shift cannot be applied when reduce is True."

        # initialize grid points:
        start = self.dim*[.0]
        n_reduced = self.dim*[0]
        g1D = self.dim*[0]
        for i, (n_i, shift_i, reduce_i, d_i) in enumerate(zip(self.n, self.shift, self.reduce, self.d)):
            if reduce_i:
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


    def integrate(self, Q:npt.NDArray):
        """Compute the integral of a scalar quantity `q` on the mesh.

        Args:
            q: scalar quantity discretised on the grid. Thus `q.shape in [self.shape, self.flat_shape]` evaluates
                to True
        Returns:
            a scalar:

        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        # TODO extend to multicomponent arrays, e.g. HFPsi
        self.flatten(q)
        return q.sum() * self.dv


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
        
    
    def is_flat(self, Q:npt.NDArray) -> bool:
        """Determine if `q` is flat.
        Args:
            Q: a quantity discretised on the grid.
        """
        return Q.shape[0] == self.flat_shape[0]        
    
    
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

    
    def get_number_of_components(self, Q:npt.NDArray):
        """compute the number of components in agrid quantity `Q`.
        Args:
            Q: a quantity discretised on the grid.
        """
        return len(Q)//self.flat_shape[0]
    
    
    def get_component(self, Q:npt.NDArray, i:int):
        """Get the i-th component of `Q`."""
        return Q[i*self.flat_shape[0]:(i+1)*self.flat_shape[0]]


    def interpolate(self, Q:npt.NDArray, r):
        """Interpolate a quantity `Q` on the mesh.

        Args:
            Q: representation of a scalar quantity on the grid.
            r: array of `p` points at which to interpolate the scalar quantity `q`. 'r.shape == (p,self.dim)'
        """
        self.flatten()
        self.flatten(Q)
        assert self.is_flat(Q)
        ijk = IJK(self.shape)
        n_gridpoints = ijk.ig_end
        n_components = self.get_number_of_components(Q) # TODO n_components should be an attribute of Q, rather than self
        i,j,k = ijk.index
        while not ijk.done():  # loop over gridpoints
            for iq in range(n_components):
                q = self.get_component(Q, iq) # TODO should be method of Q, rather than self
                # loop over interpolation points r
            ijk.inc()