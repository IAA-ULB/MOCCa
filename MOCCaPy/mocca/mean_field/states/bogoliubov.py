import numpy as np

from .hfpsi import HFPsi

class BogoliubovState:
    """This class implements a Bogoliubov mean-field state. This is the most general
    mean-field state and therefor serves as the base class for all other mean-field
    states.

    Bogoliubov mean-field states are characterised by occupancies ($\rho$) being a
    full matrix.
    """
    def __init__(self, hfpsi) -> None:
        """
        Args:
             kwargs_hfpsi: dict with all arguments needed to intialize a HFPsi object.
        """
        # Separate the HFPsi arguments in kwargs:
        self.hfpsi = hfpsi
        self.sp_energies = self.hfpsi.sp_energies
        self.occupancies = np.empty((self.hfpsi.n_total_wf, self.hfpsi.n_total_wf), dtype=float, order='F')
        #   or           = None?
