from .operator import Operator

class Overlap(Operator):
    """Compute the matrix representation of the overlap operator. <hfpsi|hfpsi>."""
    def __init__(self, hfpsi):
        super().__init__(hfpsi)
        self.O_ket = self.ket
