import numpy as np

# from MOCCaPy.mocca.mean_field.bcs import BCSState
from .nil8_f90 import nilsson
#   Note that Fortran `integer`s are 32-bit, which corresponds to `dtype=np.int32`.
#   The standard Python `int`s are 64-bit

class HFPsi:
    """Hartree-Fock wave function."""
    init_registry = {
        'nilsson' : nilsson,
        'random'  : random,
    }
    def __init__(self,
                 n_neutrons:int, n_protons:int,
                 n_neutron_wf:int, n_proton_wf:int,
                 mesh,
                 init:str='nilsson',
                 osc_freq:tuple[float]=None,
                 ):
        """Constructor for HFPsi (using the same layout as MOCCa).

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
        # modified after  subroutine iniwavefunctions in src/wavefunctions.f90 line 854
        self.mesh = mesh
        self.n_neutron_wf = n_neutron_wf
        self.n_proton_wf = n_proton_wf

        # TODO: allocate
        self.hfpsi_data = np.array([[[]]], dtype=float)
        nwt = n_neutron_wf + n_proton_wf
        #   see wavefunctions.f90 lne 883
        kparz = np.empty(nwt, dtype=np.int32)
        #   see nil8.f90 line 155
        esp1  = np.zeros(nwt, dtype=float)
        #   see nil8.f90 line 155
        meven = max(11,int(1.5*max(n_neutron_wf,n_proton_wf)**(1./3.)))
        #   see wavefunctions.f90 lne 891
        modd  = meven - 1
        #   see wavefunctions.f90 lne 89
        nwp   = self.n_proton_wf
        nwn   = self.n_neutron_wf
        npp   = n_protons
        npn   = n_neutrons
        mx,my,mz = self.mesh.mesh_shape
        #   see wavefunctions.f90 lne 869
        dx = self.mesh.d[0]
        # TODO: allocate
        spwf_map = np.array([], dtype=np.int32)

        # Argumenten voor _init_wf_nilsson
        init_kwargs = {
            'wfs'   : self.hfpsi_data,
            'kparz' : kparz,
            'esp1'  : esp1,
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
        assert init in ['nilsson', 'random'], f"Initialisation strategy '{init}' is not supported."
        if init == 'nilsson':
            assert osc_freq is not None, f"Initialisation strategy '{init}' requires sscillation frequencies."
            assert (self.mesh.d[1] == self.mesh.d[0]) and \
                   (self.mesh.d[2] == self.mesh.d[0]), \
                   f"Nilsson initialisation strategy does not support " \
                   f"meshes with different spacing on the coordinates axes."
        self.init = init
        init = self.init_registry[init]
        init(init_kwargs)

        # integer, parameter   :: Blocks                  = 8  ! This can always be fixed
        # integer              :: HFBlocks(Blocks)        = 0
        # integer              :: HFBlocks_global(Blocks) = 0
        # hfblocks_global = 0
        n_blocks = 8
        hfblocks_global = np.zeros((n_blocks,), dtype=np.int32)
        # do i=1,ININWN
        #     if(kparz(i) .gt. 0) HFBlocks_global(1) = HFBlocks_global(1) +1
        #     if(kparz(i) .lt. 0) HFBlocks_global(3) = HFBlocks_global(3) +1
        hfblocks_global[0] = (kparz[:nwn] > 0).sum() # count the 1s
        hfblocks_global[2] = (kparz[:nwn] < 0).sum() # count the -1s
        # enddo
        # do i=ININWN+1,ININWT
        #     if(kparz(i) .gt. 0) HFBlocks_global(5) = HFBlocks_global(5) +1
        #     if(kparz(i) .lt. 0) HFBlocks_global(7) = HFBlocks_global(7) +1
        # enddo
        hfblocks_global[4] = (kparz[nwn:] > 0).sum()
        hfblocks_global[6] = (kparz[nwn:] < 0).sum()
        # ! then, we are capable of figuring out the way to balance the spwfs among
        # ! the different MPI ranks.
        # call loadbalance(HFblocks_global,&
        # &                        HFblocks,spwf_map,rank_map, spwf_inverse)
        # ! now each MPI rank knows which spwfs it should grab and can allocate
        # ! the required space.
        # allocate(HFPSI(ININX*ININY*ININZ,4,sum(HFblocks))); hfpsi = 0.0d0
        #
        # ! and finally, we can call initialise_wavefunctions a second time in order
        # ! actually the requested spwfs in this partical process.
        init(init_kwargs)

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
