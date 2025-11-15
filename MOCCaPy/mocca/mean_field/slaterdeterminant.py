import numpy as np

# from .bcs import BCSState
from .hfpsi import HFPsi

class SlaterDeterminant:
    """
    Een SlaterDeterminant object heeft - from the top of my head - slechts drie cruciale datastructuren nodig.

     1. een set golffuncties (en alle meta-data die er mee geassocieerd zijn)
     2. een set occupations, i.e. een vector van reële getallen met lengte len(golffuncties)die aangeven welke
        van de golffuncties bezet zijn door deeltjes en welke leeg zijn.
     3. een set single-particle energies, i.e. een vector met lengte len(golffuncties)die aangeeft wat de
        energieën van de golffuncties zijn.
    ?
    """

    def __init__(self,
                 n_neutrons:int, n_protons:int,
                 n_neutron_wf:int, n_proton_wf:int,
                 mesh,
                 init:str='nilsson',
                 osc_freq:tuple[float]=None,
                 ):
        """Mean-field state based on a Slater determinant.

        All arguments are forwarded to HFPsi._init__().

        Args:
            n_neutrons: number of neutrons.
            n_protons: number of protons.
            n_neutron_wf: Number of neutron wave functions.
            n_proton_wf: Number of proton wave functions.
            mesh: Mesh on which the wave functions are represented.
                Currently only LagrangeMesh objects are supported.
            init: initialisation strategy for the single particle wave functions:
                'nilsson' or 'random'.
            osc_freq: optional, Oscillation frequencies for Nilsson initialisation.
        """
        self.hfpsi = HFPsi(n_neutrons=n_neutrons, n_protons=n_protons,
                           n_neutron_wf=n_neutron_wf, n_proton_wf=n_proton_wf,
                           mesh=mesh,
                           init=init, osc_freq=osc_freq
                           )
        self.occupancies = np.empty(self.hfpsi.n_total_wf, dtype=float)
        self.sp_energies = np.empty(self.hfpsi.n_total_wf, dtype=float)


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

    def energy(self, edf=None) -> float:
        """"""
        return .0

    def spwf_energy(self, potentials, edf=None) -> float:
        """"""
        return .0

    def multipole_moments(self, edf=None) -> dict[str, float]:
        """"""
        return {'?': .0}