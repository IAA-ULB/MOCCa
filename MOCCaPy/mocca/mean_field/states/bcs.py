from .bogoliubov import BogoliubovState


class BCSState(BogoliubovState):
    """This class implements a BCS mean-field state.

    BCS mean-field states are characterised by occupancies ($\rho$) being a
    diagonal matrix and thus _isa_ `BogoliubovState`.
    """


    def __init__(self, hfpsi) -> None:
        """
        Args:
            hfpsi: HFPsi object, containing the single-particle wave functions.
        """
        super().__init__(hfpsi)
