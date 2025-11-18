import numpy as np


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

    def __init__(self, hfpsi, _is_base_class=False) -> None:
        """
        Args:
             hfpsi: HFPsi object, containing the single-particle wave functions.
             _is_base_class: True if called by super().__init__(), False otherwise.
        """
        # Separate the HFPsi arguments in kwargs:
        self.hfpsi = hfpsi
        self.sp_energies = self.hfpsi.sp_energies
        if not _is_base_class:
            self.rho = np.empty((self.hfpsi.n_total_wf, self.hfpsi.n_total_wf), dtype=np.float64, order='F')

    def validate(self):
        """Verify that all conditions for representing a Bogoliubov mean-field State
        are satisfied.

        Returns:
            None
        Raises:
            AssertionError: if any of the conditions are not satisfied.
        """
        assert self.rho.shape == (self.n_total_wf,self.n_total_wf)


    @property
    def n_neutrons(self):
        return self.hfpsi.n_neutrons

    @property
    def n_protons(self):
        return self.hfpsi.n_protons

    @property
    def n_total_wf(self):
        return self.hfpsi.n_total_wf

    @property
    def n_neutron_wf(self):
        return self.hfpsi.n_neutron_wf

    @property
    def n_proton_wf(self):
        return self.hfpsi.n_proton_wf
