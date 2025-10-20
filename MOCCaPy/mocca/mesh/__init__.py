class LagrangeMesh:
    def __init__(self, n: int|tuple, d: int|float|tuple, dim: int=0,
                       bc='antiperiodic',

                 ) -> None:
        """Construct a Lagrange mesh in 1, 2 or 3 dimensions.

        Args:
            dim (int): dimension of mesh. if 0, dim is guessed as len(n)
            n: number of points in the respective dimensions. If n is an int n is the same in each direction.
            d: spacing of points in the respective dimensions. If n is a float or an int d is the same in each direction.
            bc: boundary condition type. 'antiperiodic' or 'periodic'.
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

        # initialize grid points:
