import numpy as np


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

        # compute dv (for integration)
        self.dv = np.prod(self.d) * 2 ** self.reduce.count(True)

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
            self._reshape(flat=True)
        else:
            if len(array.shape) == 1:
                # already flat
                return array
            else:
                assert array.shape == self.shape
                return array.reshape(self.flat_shape, order='F')

    def unflatten(self, array=None) -> None:
        """Unflatten the mesh array.

        Args:
            array: if None, self.gridx|y|z are unflattened. Otherwise, array is unflattened.

        Returns:
            None or the unflattened array.

        Raises:
            AssertionError: if array.shape != self.flat_shape.
        """
        if array is None:
            self._reshape(flat=False)
        else:
            if len(array.shape) == self.dim:
                # already unflattened
                return array
            else:
                assert array.shape == self.flat_shape
                return array.reshape(self.shape, order='F')

    def _reshape(self, flat:bool) -> None:
        """Reshape `self.gridx`, `self.gridy` (if existing) and `self.gridz` (if existing).

        Args:
            flat: if True the arrays are flattened. otherwise they are reshaped into self.shape.
        """
        if self.dim > 1:
            # (if self.dim == 1 flat and unflattened are equivalent, conversion is not necessary)
            # below self.dim is at least 2, so gridx and gridy exist
            if flat and len(self.gridx.shape) == self.dim:
                # gridx/y/z are not flat, and we must flatten them
                self.gridx = self.gridx.reshape(self.flat_shape, order='F')
                self.gridy = self.gridy.reshape(self.flat_shape, order='F')
                if self.dim == 3:
                    self.gridz = self.gridz.reshape(self.flat_shape, order='F')
            elif not flat and len(self.gridx.shape) == 1:
                # gridx/y/z are flat, and we must unflatten them
                self.gridx = self.gridx.reshape(self.shape, order='F')
                self.gridy = self.gridy.reshape(self.shape, order='F')
                if self.dim == 3:
                    self.gridz = self.gridz.reshape(self.shape, order='F')

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
        self.flatten()
        return function(self.gridx, self.gridy, self.gridz)

    def integrate(self, q):
        """Compute the integral of a scalar quantity `q` on the mesh.

        Args:
            q: scalar quantity discretised on the grid. Thus `q.shape in [self.shape, self.flat_shape]` evaluates
                to True
        Returns:
            a scalar:

        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        self.flatten(q)
        return q.sum() * self.dv

    def derive1(self, q):
        """Compute the 1st order derivative of a scalar quantity `q` on the mesh.

        Args:
            q: scalar quantity discretised on the grid. Thus `q.shape in [self.shape, self.flat_shape]` evaluates
                to True

        Returns:


        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        # TODO : implement

    def derive2(self, q):
        """Compute the 2nd order derivative of a scalar quantity `q` on the mesh.

        Args:
            q: scalar quantity discretised on the grid. Thus `q.shape in [self.shape, self.flat_shape]` evaluates
                to True

        Returns:


        Raises:
            AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
        """
        # TODO : implement

    # def derive3(self, q):
    #     """Compute the 3d order derivative of a scalar quantity `q` on the mesh.
    #
    #     Args:
    #         q: scalar quantity discretised on the grid. Thus `q.shape in [self.shape, self.flat_shape]` evaluates
    #             to True
    #
    #     Returns:
    #
    #
    #     Raises:
    #         AssertionError: if `not q.shape in [self.shape, self.flat_shape]`.
    #     """
    #     # TODO : implement

    def interpolate(self, q, ):
        """Interpolate a scalar quantity `q` on the mesh."""
        # TODO : implement
