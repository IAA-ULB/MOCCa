class Cartesian3D:
    def __init__(self, n:tuple[int,int,int]|int, d: float):
        if isinstance(n, int):
            self.n = (n,n,n)
        else:
            assert len(n) == 3
            self.n = n
        self.d = d