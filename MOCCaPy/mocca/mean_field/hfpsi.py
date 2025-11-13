import numpy as np

from MOCCaPy.mocca.mean_field.bcs import BCSState
from .nil8 import nilsson

class HFPsi:
    """Hartree-Fock wave function."""

    def __init__(self, n_neutron_wf:int, n_proton_wf:int,
                 mesh,
                 init:str='nilsson',
                 osc_freq:tuple[float]=None,
                 ):
        """Constructor for HFPsi (using the same layout as MOCCa).

        Args:
            n_neutron_wf: Number of neutron wave functions.
            n_proton_wf: Number of proton wave functions.
            mesh: Mesh on which the wave functions are represented.
                Currently only LagrangeMesh objects are supported.
            init: initialisation strategy for the single particle wave functions:
                'nilsson' or 'random'.
            osc_freq: optional, Oscillation frequencies for Nilsson initialisation.
        """
        # modified after  subroutine iniwavefunctions in src/wavefunctions.f90 line 854
        self.mesh = mesh
        self.n_neutron_wf = n_neutron_wf
        self.n_proton_wf = n_proton_wf
        self.hfpsi_data = None

        nwt = n_neutron_wf + n_proton_wf
        kparz = np.empty(nwt, dtype=int)
        esp1  = np.empty(nwt, dtype=float)
        meven = 0
        modd  = 0
        nwp   = self.n_proton_wf
        nwn   = self.n_neutron_wf
        npp   = 0
        npn   = 0
        mx,my,mz = self.mesh.shape
        dx = self.mesh.d[0]
        assert (self.mesh[1] == self.mesh[0]) and \
               (self.mesh[2] == self.mesh[0]), \
               f"Nilsson initialisation strategy does not support " \
               f"meshes with different spacing on the coordinates axes."
        spwf_map = None
        init_kwargs = {
            'wfs' : self.hfpsi_data,
            'kparz' : kparz,
            'esp1'  : esp1
            'meven' : meven,
            'modd'  : modd,
            'nwt'   : nwt,
            'nwp'   : nwp,
            'nwn'   : nwn,
            'npp'   : npp,
            'npn'   : npn,
            'mx'    : mx,
            'my'    : my,
            'mz'    : mz,
            'dx'    : dx,
            'osc_freq': osc_freq, 
            'spwf_map': spwf_map,
        }
        assert init in ['nilsson', 'random'], f"Initialisation strategy {init} is not supported."
        self.init = init
        if init == 'nilsson':
            self._init_wf_nilsson(init_kwargs)
        else:
            self._init_wf_random(init_kwargs)

        # wavefunctions.f90::914
        # allocate(HFPSI(ININX * ININY * ININZ, 4, sum(HFblocks))); hfpsi = 0.0d0

    def _init_wf_nilsson(self, init_kwargs:dict):
        """
        Remark: cfr src/nil8.f90 subroutine nilsson
        """
        nilsson(**init_kwargs)


    def _init_wf_random(self, init_kwargs:dict):
        """
        Remark: cfr subroutine randomspwfs in src/wavefunctions.f90 line 1006
        """
        raise NotImplementedError
