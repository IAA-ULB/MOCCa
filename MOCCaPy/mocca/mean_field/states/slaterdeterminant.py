import numpy as np

from .bcs import BCSState

class SlaterDeterminant(BCSState):
    """This class implements a Slater determinant mean-field state.

    BCS mean-field states are characterised by occupancies ($\rho$) being a
    diagonal matrix whose elements are 0 or 1, and thus _isa_ `BCSState`.
    """
    # TODO: derive from BCSState, which must derive from BogoliubovState. Move data members as needed.
    #       ?? do we need three different classes? Or do their properties emerge from the way they are
    #       computed?
    #       ?? Check issue 'mean-field state classes' #42
    #       https://realpython.com/inheritance-composition-python/
    #       Inheritance models an `is a` relationship: a square is also a rectangle, a
    #       SlaterDeterminant is also a BCSState, and a BCSState  is also BogoliubovState.
    #       A derived class is said to derive, inherit, or extend a base class
    #       So far, SlaterDeterminant is very little more than a wrapper for HFPsi...


    def __init__(self, hfpsi):
        """Mean-field state based on a Slater determinant.

        Args:
            hfpsi: HFPsi object, containing the single-particle wave functions.
        """
        super().__init__(hfpsi)


    @property
    def n_neutrons(self):
        return self.hfpsi.n_neutrons

    @property
    def n_protons(self):
        return self.hfpsi.n_protons

    @property
    def n_neutron_wf(self):
        return self.hfpsi.n_neutron_wf

    @property
    def n_proton_wf(self):
        return self.hfpsi.n_proton_wf

    @property
    def mesh(self):
        return self.hfpsi.mesh

    #

    def energy(self, edf=None) -> float:
        """"""
        return .0

    def spwf_energy(self, potentials, edf=None) -> float:
        """"""
        return .0

    def multipole_moments(self, edf=None) -> dict[str, float]:
        """"""
        return {'?': .0}