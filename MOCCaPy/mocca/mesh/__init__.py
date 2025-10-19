class LagrangeMesh:
    def __init__(self, n: int|tuple, d: int|float|tuple, dim: int=3, bc='antiperiodic') -> None:

        assert 1 <= dim <= 3
        self.dim = dim

        # dimension of the grid (nx,ny,nz) (which applies to 1/8 of the box?)
        if isinstance(n, int):
            self.n = tuple(dim*[n])
        else:
            assert isinstance(n, tuple)
            assert len(n) == dim
            for ni in n:
                assert isinstance(ni, int)
            self.n = n

        # box spacing (ints are converted to floats)
        if isinstance(d,tuple):
            assert len(d) == 3
            self.d = tuple(float(di) for di in d)
        else:
            df = float(d)
            self.d = tuple(dim*[float(d)])