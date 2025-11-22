import numpy as np

from .bcs import BCSState

class SlaterDeterminant(BCSState):
    """This class implements a Slater determinant mean-field state.

    BCS mean-field states are characterized by occupancies ($\rho$) being a
    diagonal matrix whose elements are 0 or 1, and thus _isa_ `BCSState`.
    """

    def __init__(self, hfpsi, _is_base_class=False):
        """Mean-field state based on a Slater determinant.

        Args:
            hfpsi: HFPsi object, containing the single-particle wave functions.
            _is_base_class: True if called by super().__init__(), False otherwise.
        """
        super().__init__(hfpsi, _is_base_class=True)
        if not _is_base_class:
            self.rho = np.zeros(self.n_total_wf, dtype=np.float64)

    def validate(self):
        """Verify that all conditions for representing a BCS mean-field State
        are satisfied.

        Returns:
            None
        Raises:
            AssertionError: If any of the conditions are not satisfied.
        """
        assert self.rho.shape == (self.n_total_wf, )
        assert (self.rho == 0).sum() + (self.rho == 1).sum() == self.n_total_wf


    @property
    def n_neutrons(self):
        """Return the number of neutrons in this mean-field state."""
        return self.hfpsi.n_neutrons

    @property
    def n_protons(self):
        """Return the number of protons in this mean-field state."""
        return self.hfpsi.n_protons

    @property
    def n_neutron_wf(self):
        """Return the number of neutron single particle wave functions in this mean-field state."""
        return self.hfpsi.n_neutron_wf

    @property
    def n_proton_wf(self):
        """Return the number of proton single particle wave functions in this mean-field state."""
        return self.hfpsi.n_proton_wf

    @property
    def n_total_wf(self):
        """Return the total number of single particle wave functions in this mean-field state."""
        return self.hfpsi.n_proton_wf

    @property
    def mesh(self):
        """Return the mesh on which this mean-field state is discretized."""
        return self.hfpsi.mesh

    def energy(self, edf=None) -> float:
        """"""
        return .0

    def spwf_energy(self, potentials, edf=None) -> float:
        """"""
        return .0

    def multipole_moments(self, edf=None) -> dict[str, float]:
        """"""
        return {'?': .0}