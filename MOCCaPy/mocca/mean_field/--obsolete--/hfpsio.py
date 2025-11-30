import numpy as np

from mocca.mesh.mesh_quantity import MeshQuantity
# from mocca.mean_field.bcs import BCSState
from .nil8_f90 import nilsson
from .randomspwfs_f90 import randomspwfs
# TODO: both need adaptations if MPI is used

#   Note that Fortran `integer`s are 32-bit, which corresponds to `dtype=np.int32`.
#   The standard Python `int`s are 64-bit

class HFPsiO(MeshQuantity):
    """Hartree-Fock wave function."""

    # TODO: Maybe HFPsi needs to derive from MeshQuantity, so that we can take lagrange derivatives
    #       This may require an extension of MeshQuantity as HFPsi.data has shape (mesh.linear_size,
    #       4, n_spwf), and MeshQuantity does only provide one dimension for the components.

    init_f90_registry = {
        'nilsson' : nilsson,
        'random'  : randomspwfs,
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
            init: initialization strategy for the single particle wave functions:
                'nilsson' or 'random'.
            osc_freq: optional, Oscillation frequencies for Nilsson initialization.
        """
        # This code is messy because it relies on F90 modules nil8_f90 and randomspwfs_f90, bad design
        # decisions in there percolate upwards
        super.__init__()

        self.n_neutrons = n_neutrons
        self.n_protons  = n_protons
        self.n_neutron_wf = n_neutron_wf
        self.n_proton_wf  = n_proton_wf
        self.mesh = mesh
        self.init = init
        self.osc_freq = osc_freq
        if init in self.init_f90_registry:
            self._init_f90()
        else:
            raise ValueError(f"Currently only f90 initialization strategies ('nilsson', 'random') are supported. "
                             f"Initialization strategy '{init}' is not recognized."
                             )
        self.n_spwf = self.data.shape[2]

    def _init_f90(self):
        """
        Do the actual initialization of self.data.
        """
        # Modified after  subroutine iniwavefunctions in src/wavefunctions.f90 line 854
        # This code is probably far too complicated, but following the original MOCCa implementation.

        self.data = np.array([[[]]], dtype=np.float64)
        nwt = self.n_neutron_wf + self.n_proton_wf
        #   see wavefunctions.f90 lne 883
        kparz = np.empty(nwt, dtype=np.int32)
        #   see nil8.f90 line 155
        esp1  = np.zeros(nwt, dtype=np.float64)
        #   see nil8.f90 line 155
        meven = max(11,int(1.5*max(self.n_neutron_wf, self.n_proton_wf)**(1./3.)))
        #   see wavefunctions.f90 lne 891
        modd  = meven - 1
        #   see wavefunctions.f90 lne 89
        nwp   = self.n_proton_wf
        nwn   = self.n_neutron_wf
        npp   = self.n_protons
        npn   = self.n_neutrons
        mx,my,mz = self.mesh.mesh_shape
        #   see wavefunctions.f90 lne 869
        dx = self.mesh.d[0]
        spwf_map = np.array([], dtype=np.int32)
        if self.osc_freq is None:
            self.osc_freq = np.empty(self.mesh.dim, dtype=np.int32) # dummy argument for randomspwfs
        # Arguments for nilsson and randomspwfs
        # TODO: it is a bad design decision to pass the same arguments
        #       to both functions. Fix?
        #       randomspwfs only needs: (input) nwt, nwp, nwn, spwf_map,
        #       (output) psi, spe, par
        init_args = [
            self.data,
            kparz,
            esp1,
            meven,
            modd,
            nwt,
            nwp,
            nwn,
            npp,
            npn,
            mx,
            my,
            mz,
            dx,
            self.osc_freq,
            spwf_map,
        ]

        assert self.init in ['nilsson', 'random'], f"Initialization strategy '{init}' is not supported."
        if self.init == 'nilsson':
            assert self.osc_freq is not None, f"Initialization strategy '{init}' requires sscillation frequencies."
            assert (self.mesh.d[1] == self.mesh.d[0]) and \
                   (self.mesh.d[2] == self.mesh.d[0]), \
                   f"Nilsson initialization strategy does not support " \
                   f"meshes with different spacing on the coordinates axes."
        init = self.init_f90_registry[self.init] # fetch the method from the init_registry
        init(*init_args)

        # integer, parameter   :: Blocks                  = 8  ! This can always be fixed
        # integer              :: HFBlocks(Blocks)        = 0
        # integer              :: HFBlocks_global(Blocks) = 0
        # hfblocks_global = 0
        hfblocks_global = np.zeros(8, dtype=np.int32)
        hfblocks_local  = np.zeros(8, dtype=np.int32)
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
        ##< following code adapted from subroutine loadbalance in wavfunctions.f90
        # call loadbalance(
        #    , HFblocks_global  ! intent(in)
        #    , HFblocks         ! intent(out)
        #    , spwf_map         ! intent(out), allocatable
        #    , rank_map         ! intent(out), allocatable
        #    , spwf_inverse     ! intent(out), allocatable
        #    )
        # integer              :: ranks_per_block(Blocks)
        ranks_per_block = np.zeros(8, dtype=np.int32)

        # Nspwf  = sum(blocks_global)
        n_spwf = hfblocks_global.sum()

        # if(.not.allocated(rank_map)) allocate(rank_map(Nspwf), spwf_inverse(Nspwf))
        # rank_map     = 0 ; spwf_inverse = 0 ; blocks_local = 0; ranks_per_block = 0
        rank_map     = np.zeros(n_spwf, dtype=np.int32)
        spwf_inverse = np.zeros(n_spwf, dtype=np.int32)

        max_spwf_per_rank = np.int32(100_000_000) ### see geninfo.f90 line 151 (default value)
        nprocs = 1                                ### see geninfo.f90 line 146 (default value)
        mpi_rank = 0                              ### see geninfo.f90 line 146 (default value)

        activeblocks = 0
        # do B=1,8
        for ib in range(8):
        #   if(blocks_global(B) .ne. 0) then
            if hfblocks_global[ib] != 0:
        #     activeblocks = activeblocks + 1
                activeblocks += 1
        #     ! ensure that every active block gets
        #     !  (1) sufficient processes such that no process should go above max_spwf_per_rank
        #     !  (2) at least one attributed process
        #     !
        #     ! Important note: the code does not strictly enforce max_spwf_per_rank because
        #     ! of the rounding to nearest integer below; max_spwf_per_rank should be understood
        #     ! more as a rough guideline.
        #     ranks_per_block(B) = max(1,ceiling(blocks_global(B)/(1.0d0*max_spwf_per_rank)))
                ranks_per_block[ib] = np.max(1,np.ceil(hfblocks_global[ib] // max_spwf_per_rank))
        #   endif
        # enddo
        # already_assigned = sum(ranks_per_block)
        already_assigned = ranks_per_block.sum()

        # if(mod(activeblocks, NPROCS) .ne. 0) then
        #   call stp('Incompatible number of MPI ranks.')
        # endif
        assert activeblocks % nprocs == 0, "Incompatible number of MPI ranks"

        # blocks_per_rank = activeblocks / NPROCS
        blocks_per_rank = activeblocks // nprocs

        # block_count = -1 ! unintuitive starting point: first block will be '0'
        block_count = -1
        # do B=1,8
        for ib in range(8):
        #   if(blocks_global(B) .eq. 0) cycle
            if hfblocks_global[ib] == 0: continue
        #   block_count = block_count + 1
            block_count += 1
        #   if( block_count / blocks_per_rank .eq. MPI_rank) then
            if block_count // blocks_per_rank == mpi_rank:
        #     ! attention, INTEGER division in the line above
        #     blocks_local(B) = blocks_global(B)
                hfblocks_local[ib] = hfblocks_global[ib]
        #   endif
        # enddo

        # ! Find the first non-zero size in blocks_local
        # do B=1,8
        for ib in range(8):
        #   if(blocks_local(B) .ne. 0) exit
            if hfblocks_local[ib] != 0: break
        # enddo
        # offset = sum(blocks_global(1:B-1))
        offset = hfblocks_global[:ib].sum()
        # ! Calculate the spwf <-> rank mapping and its inverse
        # allocate(spwf_map(sum(blocks_local)))
        n_spwf_local = hfblocks_local.sum()
        spwf_map = np.empty(n_spwf_local, dtype=np.int32)
        # do i=1,sum(blocks_local)
        for i in range(n_spwf_local):
        #   spwf_map(i)            = offset + i
            spwf_map[i] = offset+i
        #   spwf_inverse(offset+i) = i
            spwf_inverse[offset + i] = i
        #   rank_map(offset+i)     = MPI_RANK
            rank_map[offset + i] = mpi_rank
        # enddo
        ##>

        # ! now each MPI rank knows which spwfs it should grab and can allocate
        # ! the required space.
        # allocate(HFPSI(ININX*ININY*ININZ,4,sum(HFblocks))); hfpsi = 0.0d0
        self.data = np.zeros((self.mesh.linear_size,4,n_spwf_local), dtype=np.float64, order='F')
        #
        # ! and finally, we can call initialise_wavefunctions a second time in order
        init(*init_args)
        # ! actually the requested spwfs in this partical process.
