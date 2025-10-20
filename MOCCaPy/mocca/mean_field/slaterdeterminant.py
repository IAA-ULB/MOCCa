from .bcs import BCSState


class SlaterDeterminant(BCSState):
    """

    """
    wf_initializers = {
        'Nilsson' : create_wf_nilsson,
        'random'  : create_wf_random,
    }
    def __init__(self, Z:int, N:int, mesh, nwp:int, nwn:int, wf_init='Nilsson'):
        super().__init__()
        self.Z = Z
        self.N = N
        self.mesh = mesh
        self.nwp = nwp
        self.nwn = nwn
        self.wf_init_function = self.__class__.wf_initializers[wf_init]
        self.wf = self.wf_init_function()


    def energy(self, edf=None) -> float:
        """"""
        return .0

    def spwf_energy(self, potentials, edf=None) -> float:
        """"""
        return .0

    def multipole_moments(self, edf=None) -> dict[str, float]:
        """"""
        return {'?': .0}