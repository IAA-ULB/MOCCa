def heavy_ball():
    """"""
    
def linear_mix():
    """"""
    
class MFESolver:
    # registry for spwf algorithms
    spwf_algos = {
        'heavy_ball' : heavy_ball,
    }
    # registry for scf algorithms
    scf_algos = {
        'linear_mix' : linear_mix,
    }
    def __init__(self,
        wf0, edf,
        spwf_algo:str = 'heavy_ball',
        scf_algo:str = 'linear_mix',
        energy_tol:float = 1e-9,
        moment_tol:float = 1e-3,
    ):
        """

        Raises:
            KeyError if spwf_algo|scf_algo is not found in the corresponding registries
        """
        self.wf0 = wf0
        self.edf = edf
        self.spwf_algo = self.__class__.spwf_algos[spwf_algo]
        self.scf_algo = self.__class__.scf_algos[scf_algo]
        self.energy_tol = energy_tol
        self.moment_tol = moment_tol

    def solve(self):
        """Solve the mean-field equations.

        Returns:
            wf
            densities
            potentials
        """
        wf = self.wf0
        densities = None
        potentials = None
        return wf, densities, potentials

    def write_hdf5(self, hdf5_file) -> None:
        """"""