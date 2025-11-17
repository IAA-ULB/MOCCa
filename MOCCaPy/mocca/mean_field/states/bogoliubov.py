import numpy as np

from .hfpsi import HFPsi

class BogoliubovState:
    """This class implements a Bogoliubov mean-field state. This is the most general
    mean-field state and therefor serves as the base class for all other mean-field
    states.

    Bogoliubov mean-field states are characterised by occupancies ($\rho$) being a
    full matrix.
    """
    # Rather than accepting the arguments of HFPsi.__init__ and forwarding them,
    # We prefer to create a HFPsi object and hand it to BogoliubovState. This
    # simplifies handling the arguments, and separates responsibilities.

    def __init__(self, hfpsi) -> None:
        """
        Args:
             hfpsi: HFPsi object, containing the single-particle wave functions.
        """
        # Separate the HFPsi arguments in kwargs:
        self.hfpsi = hfpsi
        self.sp_energies = self.hfpsi.sp_energies
        self.occupancies = np.empty((self.hfpsi.n_total_wf, self.hfpsi.n_total_wf), dtype=float, order='F')
        #   or           = None?

    @property
    def n_neutrons(self):
        return hfpsi.n_neutrons

    @property
    def n_protons(self):
        return hfpsi.n_protons

    @property
    def n_total_wf(self):
        return hfpsi.n_total_wf

    @property
    def n_neutron_wf(self):
        return hfpsi.n_neutron_wf

    @property
    def n_proton_wf(self):
        return hfpsi.n_proton_wf
