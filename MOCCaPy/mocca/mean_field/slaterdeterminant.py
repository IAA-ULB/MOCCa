from .bcs import BCSState


class SlaterDeterminant(BCSState):
    """

    """
    # Registry for wave function initializers
    wf_initializers = {
        'Nilsson' : create_wf_nilsson,
        'random'  : create_wf_random,
    }
    def __init__(self, Z:int, N:int, mesh, nwp:int, nwn:int, wf_init='Nilsson'):
        """Mean-field state base on a Slater determinant.

        Args:
            Z (int): Number of protons.
            N (int): Number of neutrons.
            mesh: mesh on which the mean-field state is discretised. Ccurrently, only Lagrange meshes (LagrangeMesh)
                are supported.
            nwp (int): Number of proton single particle wave functions. (?)
            nwn (int): Number of neutron single particle wave functions. (?)
            wf_init (str): Initialization method for the wave funtion. 'Nilsson' or 'random'.
        """
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