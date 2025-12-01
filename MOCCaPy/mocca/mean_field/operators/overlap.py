from .operator import Operator

class Overlap(Operator):
    """Compute the matrix representation of the overlap operator. <hfpsi|hfpsi>."""
    def __init__(self, mfs):
        super().__init__(mfs)

    def __repr__(self):
        return f"<Overlap(Operator)[{self.operand_data.mfs}]>"

    def apply_base_operator(self):
        self.operand_data.O_ket = self.operand_data.ket
