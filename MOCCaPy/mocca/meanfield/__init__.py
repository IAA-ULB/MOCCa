def init_nilsson(wf) -> None:
    """"""


def init_random(wf) -> None:
    """"""


class BogoliubovState:
    """

    """
    def __init__(self) -> None:
        """"""


class BCSState(BogoliubovState):
    """

    """
    def __init__(self):
        super().__init__()

# so far we are putting all functionality in SlaterDeterminant and non in the base classes
class SlaterDeterminant(BCSState):
    """

    """
    wf_initializers = {
        'Nilsson' : init_nilsson,
        'random' : init_random,
    }
    def __init__(self, Z:int, N:int, mesh, nwp:int, nwn:int, wf_init='Nilsson'):
        super().__init__()
        self.Z = Z
        self.N = N
        self.mesh = mesh
        self.nwp = nwp
        self.nwn = nwn
        self.wf_init_function = self.__class__.wf_initializers[wf_init]


    def energy(self, edf=None) -> float:
        """"""
        return .0

    def spwf_energy(self, potentials, edf=None) -> float:
        """"""
        return .0

    def multipole_moments(self, edf=None) -> dict[str, float]:
        """"""
        return {'?': .0}