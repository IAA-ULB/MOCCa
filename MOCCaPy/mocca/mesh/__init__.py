class Cartesian3D:
    def __init__(self, n: int|tuple[int,int,int], d: int|float|tuple[float,float,float]):
        
        # dimension of the grid (nx,ny,nz) (which applies to 1/8 of the box?)
        if isinstance(n, int):
            self.n = (n,n,n)
        else:
            assert len(n) == 3
            self.n = n

        # box spacing (ints are converted to floats)
        if isinstance(d,tuple):
            assert len(d) == 3
            self.d = tuple(float(di) for di in d)
        else:
            df = float(d)
            self.d = (df,df,df)