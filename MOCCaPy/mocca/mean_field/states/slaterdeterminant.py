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
            self.rho = np.zeros(self.hfpsi.n_total_wf, dtype=np.float64)

    def __repr__(self):
        return f"<SlaterDeterminant(BCSState)[{self.hfpsi}]>"

    def dbg_assert(self):
        """Assert some conditions that may indicate bugs when `False`."""

        print(f"\nSlaterDeterminant.dbg_assert() called on instance `{self}`")
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
        return self.hfpsi.n_total_wf

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

    def occupancies(self, h_ii) -> np.ndarray:
        """Compute the occupancies of this mean-field state.

        Args:
            h_ii: diagonal elements of the Hamiltonian.

        Returns:
            rho, rho2, order_neutrons, order_protons
            rho: the occupancies
            rho2: number of neutrons/protons in the current level and all levels below.
            order_neutrons: array of indices that sort the neutron part of h_diag
            order_protons: array of indices that sort the proton part of h_diag
            Args:
        """
        return _sd_occupancies(self.hfpsi.n_neutrons, self.hfpsi.n_protons, h_ii, self.hfpsi.n_neutron_wf)


def _sd_occupancies(n_neutrons, n_protons, h_diag, nwn):
    """Compute SlaterDeterminant occupancies, given these arguments:

    Args:
        n_neutrons: number of neutrons
        n_protons: number of protons
        h_diag: diagonal elements of the hamiltonian
        nwn: number of neutron wave functions

    Returns:
        rho, rho2, order_neutrons, order_protons
        rho: the occupancies
        rho2: number of neutrons/protons in the current level and all levels below.
        order_neutrons: array of indices that sort the neutron part of h_diag
        order_protons: array of indices that sort the proton part of h_diag
    """
    rho  = np.zeros_like(h_diag)
    rho2 = np.zeros_like(h_diag)
    order_neutrons = np.argsort(h_diag[:nwn])
    order_protons  = np.argsort(h_diag[nwn:]) + nwn

    n = n_neutrons // 2
    rho [order_neutrons[:n]] = 1.0
    rho2[order_neutrons[:n]] = 1.0
    if n_neutrons % 2 == 1:
        rho[order_neutrons[n+1]] = 1.0
    p = n_protons // 2
    rho [order_protons [:p]] = 1.0
    rho2[order_protons [:p]] = 1.0
    if n_protons % 2 == 1:
        rho[order_protons [p+1]] = 1.0
    rho2 += rho
    _sum = 0.
    for i in range(nwn):
        _sum += rho2[order_neutrons[i]]
        rho2[order_neutrons[i]] = _sum
    _sum = .0
    for i in range(nwn,h_diag.size):
        _sum += rho2[order_protons[i-nwn]]
        rho2[order_protons[i-nwn]] = _sum
    return rho, rho2, order_neutrons, order_protons
