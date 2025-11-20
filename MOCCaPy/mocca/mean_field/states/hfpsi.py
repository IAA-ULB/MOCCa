import numpy as np

# from MOCCaPy.mocca.mean_field.bcs import BCSState
from mocca.mesh.observable import Observable
from mocca.mean_field.operators import Overlap

from mocca.f90.nil8_f90 import nilsson
from mocca.f90.randomspwfs_f90 import randomspwfs
#   these routines are pretty much black box to me...
# TODO: both need adaptations if MPI is used

#   Note that Fortran `integer`s are 32-bit, which corresponds to `dtype=np.int32`,
#   while the standard Python `int`s are 64-bit

# TODO: test
class NumpyWfInitializer:
    def __init__(self,
                 n_neutrons: int, n_protons: int,
                 n_neutron_wf: int, n_proton_wf: int,
                 mesh,
                 ):
        """
        A Python implementation of randomspwfs.f90.
        """

        self.n_neutrons = n_neutrons
        self.n_protons = n_protons
        self.n_neutron_wf = n_neutron_wf
        self.n_proton_wf = n_proton_wf
        self.mesh = mesh

    def __call__(self):
        """Do the initialization.
        Returns:
            hfpsi, hfblocks, sp_energies
        """
        nwn = self.n_neutron_wf
        nwp = self.n_proton_wf
        nwt = nwn + nwp
        par = np.empty(nwt, dtype=np.float64)
        par[          :nwn//2    ] = +1
        par[nwn//2    :nwn       ] = -1
        par[nwn       :nwn+nwp//2] = +1
        par[nwn+nwp//2:          ] = -1

        hfblocks  = np.zeros(8, dtype=np.int32)
        hfblocks[0] = (par[:nwn] > 0).sum() # count the 1s
        hfblocks[2] = (par[:nwn] < 0).sum() # count the -1s
        hfblocks[4] = (par[nwn:] > 0).sum()
        hfblocks[6] = (par[nwn:] < 0).sum()

        spwf_map = np.arange(0, nwt, dtype=np.int32)

        rng = np.random.default_rng()
        hfpsi = rng.random(self.mesh.linear_size * 4 * nwt)
        hfpsi = hfpsi.reshape((self.mesh.linear_size, 4, nwt), order='F')

        sp_energies = np.empty(nwt, dtype=np.float64)
        sp_energies.fill(100)

        return hfpsi, hfblocks, sp_energies


class F90WfInitializer:
    def __init__(self,
                 n_neutrons:int, n_protons:int,
                 n_neutron_wf:int, n_proton_wf:int,
                 mesh,
                 init:str,
                 osc_freq:tuple[float]=None,
                 ):
        """

        """
        assert init in ('nilsson', 'randomspwfs'), f"F90 initialisation strategy '{init}' is not recognizeded."

        assert mesh.dim == 3, \
            f"F90 initialisation strategies support only 3D meshes, not {mesh.dim}D."
        assert (mesh.d[1] == mesh.d[0]) and \
               (mesh.d[2] == mesh.d[0]), \
            f"Nilsson initialisation strategy does not support " \
            f"meshes with different spacing on the coordinates axes."

        if init == 'nilsson':
            assert osc_freq is not None
            if not isinstance(osc_freq, np.ndarray):
                osc_freq = np.array(osc_freq, dtype=np.float64)
            assert osc_freq.shape == (3,), "Expecting ndarray of shape (3,) and dtype=np.float64"
            assert osc_freq.dtype == float, "Expecting ndarray of shape (3,) and dtype=np.float64"

        if init == 'randomspwfs' and osc_freq is None:
            osc_freq = np.empty(mesh.dim, dtype=np.float64)

        self.n_neutrons = n_neutrons
        self.n_protons  = n_protons
        self.n_neutron_wf = n_neutron_wf
        self.n_proton_wf  = n_proton_wf
        self.mesh = mesh
        self.osc_freq = osc_freq
        self.init = nilsson if (init == 'nilsson') else \
            randomspwfs

    def __call__(self):
        """Do the initialization.
        Returns:
            hfpsi, hfblocks, esp1
        """
        hfpsi = np.array([[[]]], dtype=np.float64)
        self.n_total_wf = self.n_neutron_wf + self.n_proton_wf
        #   see wavefunctions.f90 lne 883
        kparz = np.empty(self.n_total_wf, dtype=np.int32)
        #   see nil8.f90 line 155
        esp1  = np.zeros(self.n_total_wf, dtype=np.float64)
        #   see nil8.f90 line 155
        meven = max(11,int(1.5*max(self.n_neutron_wf, self.n_proton_wf)**(1./3.)))
        #   see wavefunctions.f90 lne 891
        modd  = meven - 1
        #   see wavefunctions.f90 lne 89
        spwf_map = np.array([], dtype=np.int32)

        # Arguments for nilsson and randomspwfs
        # TODO: it is a bad design decision to pass the same arguments
        #       to both functions. Fix?
        #       randomspwfs only needs: (input) nwt, nwp, self.n_neutron_wf, spwf_map,
        #       (output) psi, spe, par
        init_args = [
            hfpsi,
            kparz,
            esp1,
            meven,
            modd,
            self.n_total_wf,
            self.n_proton_wf,
            self.n_neutron_wf,
            self.n_protons,
            self.n_neutrons,
            self.mesh.mesh_shape[0],
            self.mesh.mesh_shape[1],
            self.mesh.mesh_shape[2],
            self.mesh.d[0],
            self.osc_freq,
            spwf_map,
        ]


        self.init(*init_args)

        hfblocks  = np.zeros(8, dtype=np.int32)
        hfblocks[0] = (kparz[:self.n_neutron_wf] > 0).sum() # count the 1s
        hfblocks[2] = (kparz[:self.n_neutron_wf] < 0).sum() # count the -1s
        hfblocks[4] = (kparz[self.n_neutron_wf:] > 0).sum()
        hfblocks[6] = (kparz[self.n_neutron_wf:] < 0).sum()

        spwf_map = np.arange(0, self.n_total_wf, dtype=np.int32)

        hfpsi = np.zeros((self.mesh.linear_size, 4, self.n_total_wf), dtype=np.float64, order='F')

        self.init(*init_args)

        return hfpsi, hfblocks, esp1


class HFPsi(Observable):
    """ Data structure for Hartree-Fock wave function."""

    @classmethod
    def like(cls, hfpsi):
        return cls(
            init=hfpsi,
            n_neutrons=None, n_protons=None,
            n_neutron_wf=None, n_proton_wf=None,
            mesh=None,
        )

    def __init__(self,
                 n_neutrons:int, n_protons:int,
                 n_neutron_wf:int, n_proton_wf:int,
                 mesh,
                 init:str='nilsson',
                 osc_freq:tuple[float]=None,
                 ):
        """Constructor for HFPsi (using the same data structure as in MOCCa).

        Args:
            n_neutrons: number of neutrons.
            n_protons: number of protons.
            n_neutron_wf: Number of neutron wave functions.
            n_proton_wf: Number of proton wave functions.
            mesh: Mesh on which the wave functions are represented.
                Currently only LagrangeMesh objects are supported.
            init: initialisation strategy for the single particle wave functions:
                'nilsson', 'randomspwfs' (which recycle MOCCa f90 code), or
                'np.random' (which is entirely based on Numpy).
            osc_freq: Oscillation frequencies, optional, only required for Nilsson initialisation.
        """
        init_strategies = ['nilsson', 'randomspwfs', 'np.random']
        if init not in init_strategies:
            if isinstance(init, HFPsi):
                # init is another HFPsi instance. copy its variables and allocate empty memory
                self.init = init

                self.n_neutrons   = init.n_neutrons
                self.n_protons    = init.n_protons
                self.n_total_wf   = init.n_total_wf
                self.n_neutron_wf = init.n_neutron_wf
                self.n_proton_wf  = init.n_proton_wf
                self.hfblocks     = init.hfblocks
                self.hfblockrange = init.hfblockrange
                super().__init__(mesh=init.mesh, data=np.empty_like(init.data), symmetry=init.symmetry)
                return

            elif init is None:
                # only for test purposes
                self.n_neutrons   = n_neutrons
                self.n_protons    = n_protons
                self.n_total_wf   = n_neutron_wf + n_proton_wf
                self.n_neutron_wf = n_neutron_wf
                self.n_proton_wf  = n_proton_wf
                self.hfblocks     = 8*[None]
                self.hfblockrange = 8*[None]
                self.data = np.empty((mesh.linear_size,4,self.n_total_wf), dtype=np.float64, order='F')
                self.symmetry = np.ones((4*self.n_total_wf,mesh.dim), dtype=np.int32, order='F')
                super().__init__(
                    mesh=mesh,
                    data=np.empty((mesh.linear_size,4,self.n_total_wf), dtype=np.float64, order='F'),
                    symmetry=np.ones((4*self.n_total_wf,mesh.dim), dtype=np.int32, order='F'),
                    name='hfpsi'
                )
                print(f"\nWARNING: Initialisation strategy {None} is only for testing purposes.\n")
                return
            else:
                raise ValueError(f"Unknown initialisation strategy '{init}'. Expecting one of {init_strategies} or HFPsi instance.")
        initializer = \
            F90WfInitializer(
                n_neutrons=n_neutrons, n_protons=n_protons, n_neutron_wf=n_neutron_wf, n_proton_wf=n_proton_wf,
                mesh=mesh, init=init, osc_freq=osc_freq
            ) if init in ['nilsson', 'randomspwfs'] else\
            NumpyWfInitializer(
                n_neutrons=n_neutrons, n_protons=n_protons, n_neutron_wf=n_neutron_wf, n_proton_wf=n_proton_wf, mesh=mesh
            )
        self.init = init

        self.n_neutrons = n_neutrons
        self.n_protons  = n_protons
        self.n_total_wf = n_neutron_wf + n_proton_wf
        self.n_neutron_wf = n_neutron_wf
        self.n_proton_wf  = n_proton_wf


        data, self.hfblocks, self.sp_energies = initializer()
        self.hfblockrange = 8*[None]

        istart = 0
        istop = self.hfblocks[0]
        for i in range(7):
            self.hfblockrange[i] = (istart, istop)
            istart = istop
            istop += self.hfblocks[i+1]
        self.hfblockrange[7] = (istart, self.n_total_wf)

        # TODO: adapt for 1D and 2D meshes
        # allocate(sx(4,sum(hfblocks)), sy(4,sum(hfblocks)), sz(4,sum(hfblocks)))
        symmetry = np.empty((4 * self.n_total_wf, mesh.dim), dtype=np.int32)

        # do i=1, HFBlocks(1)
        #     sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        #     sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1
        #     sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        #     sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
        # enddo
        iblock = 0
        for i_spwf in range(*self.hfblockrange[iblock]):
            symmetry[4*i_spwf:4*i_spwf + 4, 0] = ( 1,-1,-1, 1 ) # sx
            symmetry[4*i_spwf:4*i_spwf + 4, 1] = ( 1,-1, 1,-1 ) # sy
            symmetry[4*i_spwf:4*i_spwf + 4, 2] = ( 1, 1,-1,-1 ) # sz

        # do i=HFBlocks(1) + 1,HFBlocks(1) + HFBlocks(3)
        #     sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        #     sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1
        #     sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        #     sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
        # enddo
        iblock += 2
        for i_spwf in range(*self.hfblockrange[iblock]):
            symmetry[4*i_spwf:4*i_spwf + 4, 0] = ( 1,-1,-1, 1 ) # sx
            symmetry[4*i_spwf:4*i_spwf + 4, 1] = ( 1,-1, 1,-1 ) # sy
            symmetry[4*i_spwf:4*i_spwf + 4, 2] = (-1,-1, 1, 1 ) # sz

        # do i=HFBlocks(1) + HFBlocks(3)+1,HFBlocks(1) + HFBlocks(3) +HFBlocks(5)
        #     sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = +1
        #     sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = +1
        #     sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = -1
        #     sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = -1
        # enddo
        iblock += 2
        for i_spwf in range(*self.hfblockrange[iblock]):
            symmetry[4*i_spwf:4*i_spwf + 4, 0] = ( 1,-1,-1, 1 ) # sx
            symmetry[4*i_spwf:4*i_spwf + 4, 1] = ( 1,-1, 1,-1 ) # sy
            symmetry[4*i_spwf:4*i_spwf + 4, 2] = ( 1, 1,-1,-1 ) # sz

        # do i=HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + 1,                      &
        # &       HFBlocks(1)+HFBlocks(3)+HFBlocks(5) + HFBLocks(7)
        #     sx(1,i) =  1 ; sy(1,i) = +1 ; sz(1,i) = -1
        #     sx(2,i) = -1 ; sy(2,i) = -1 ; sz(2,i) = -1
        #     sx(3,i) = -1 ; sy(3,i) = +1 ; sz(3,i) = +1
        #     sx(4,i) =  1 ; sy(4,i) = -1 ; sz(4,i) = +1
        # enddo
        iblock += 2
        for i_spwf in range(*self.hfblockrange[iblock]):
            symmetry[4*i_spwf:4*i_spwf + 4, 0] = ( 1,-1,-1, 1 ) # sx
            symmetry[4*i_spwf:4*i_spwf + 4, 1] = ( 1,-1, 1,-1 ) # sy
            symmetry[4*i_spwf:4*i_spwf + 4, 2] = (-1,-1, 1, 1 ) # sz

        super().__init__(name='HFPsi', mesh=mesh, data=data, symmetry=symmetry)

    def clone(self):
        """return an empty copy of self."""
        return HFPsi(
            init=self,
            n_neutrons=None, n_protons=None,
            n_neutron_wf=None, n_proton_wf=None,
            mesh=None,
        )


    def ilc(self, i4, i_wf):
        """Return the linear component index from the wave function component index `i4`
        (`0<=i4<4`) and the wave function index `i_wf`.

        In MOCCa HFPsi has shape `(mesh.linear_size, 4, self.n_total_wf)`. MOCCaPy Observables,
        however, are more comfortable with `(mesh.linear_size, 4 * self.n_total_wf)`. This
        method converts a MOCCa index `(i4,i_wf)` to a MOCCaPy index `iq = i4 + 4*i_wf`.

        Args:
            i4: the wave function component index of the part of the wave function, `0 <= i4 < 4`.
                This index refers to the real/imaginary spin-up (`i4` = 0, 1) and real/imag
                spin-down components (`i4` = 2, 3) of the wave function with index `i_wf`.
            i_wf: index of the wave function `0 <= i_wf < self.n_total_wf`.
        Returns:
            i4 + 4* i_wf
        """
        return 4*i_wf + i4

    def allocate_matrix_representation(self):
        """The matrix representation of an operator is a list of 8 square matrices,
        corresponding to the 8 symmetry blocks of the mean-field state.
        Its size is equal to the number of single particle wave functions in the symmetry
        blocks.

        Returns:
            a list of 8 square matrices, corresponding to the 8 symmetry blocks of the mean-field.
        """
        blocks = 8*[None]
        for ib in range(8):
            blockrange = self.hfblockrange[ib]
            n_spwfs = blockrange[1] - blockrange[0]
            blocks[ib] = np.empty((n_spwfs, n_spwfs), dtype=np.float64, order='F')

        return blocks

    def normalize(self):
        """Normalize the HFPsi state."""
        overlap = Overlap(self)
        norm = overlap.compute_diagonal_elements()
        np.sqrt(norm, out=norm)
        hfpsi_data3 = self.data.reshape((self.data.shape[0], 4, self.data.shape[1]//4), order='F')
        for i_spwf in range(norm.size):
            hfpsi_data3[:,:,i_spwf] /= norm[i_spwf]
