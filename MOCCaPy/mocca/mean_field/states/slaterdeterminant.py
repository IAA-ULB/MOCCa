from .bcs import BCSState

class SlaterDeterminant(BCSState):
    """This class implements a Slater determinant mean-field state.

    BCS mean-field states are characterised by occupancies ($\rho$) being a
    diagonal matrix whose elements are 0 or 1, and thus _isa_ `BCSState`.
    """

    def __init__(self, hfpsi):
        """Mean-field state based on a Slater determinant.

        Args:
            hfpsi: HFPsi object, containing the single-particle wave functions.
        """
        super().__init__(hfpsi)


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
        """Return the mesh on which this mean-field state is discretised."""
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