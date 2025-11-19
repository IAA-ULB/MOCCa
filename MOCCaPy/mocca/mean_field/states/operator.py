import numpy as np

from .hfpsi import HFPsi


class Operator:
    def __init__(self, hfpsi):
        self.ket = hfpsi
        self.bra = hfpsi
        self.matrix = hfpsi.allocate_matrix_representation()

    def compute_derivatives(self):
        pass

    def add_local_terms(self):
        pass

    def add_non_local_terms(self):
        pass

    def compute_matrix_representation(self):
        self.compute_derivatives()
        self.add_local_terms()
        self.add_non_local_terms()

        # Assuming no symmetry
        O_ket = self.O_ket
        bra = self.bra
        for ib in range(8):
            block = self.matrix[ib]
            for j in range(*self.bra.hfblockrange[ib]):
                for i in range(*self.ket.hfblockrange[ib]):
                    # TODO: use out= ? for performance
                    block[i-self.ket.hfblockrange[ib][0],j-self.bra.hfblockrange[ib][0]] = (
                        bra.data[:,   4*j] * O_ket.data[:,   4*i] +
                        bra.data[:, 1+4*j] * O_ket.data[:, 1+4*i] +
                        bra.data[:, 2+4*j] * O_ket.data[:, 2+4*i] +
                        bra.data[:, 3+4*j] * O_ket.data[:, 3+4*i]
                    ).sum(0)
        return self.matrix


class Overlap(Operator):
    def __init__(self, hfpsi):
        super().__init__(hfpsi)
        self.O_ket = hfpsi


class HamiltonianWoodsSaxon(Operator):
    def __init__(self, hfpsi, hbm=20.73553000, V0=-50.0, r0=1.25, a=0.5):
        super().__init__(hfpsi)
        self.hbm = hbm
        self.V0 = V0
        self.ainv = 1/a
        self.R = r0 * np.pow(hfpsi.n_neutrons + hfpsi.n_protons, 1/3)

    def compute_derivatives(self):
        if self.ket.mesh.dim >= 2:
            self.ket.differentiate('Laplacian')
        else:
            self.ket.differentiate('xx')

    def add_local_terms(self):
        self.O_ket = HFPsi.like(self.ket)

        def VWoodsSaxon1D(r):
            return self.V0 / (1. + np.exp(self.ainv * (r - self.R)))

        def VWoodsSaxon2D(x,y):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x*x + y*y) - self.R)))

        def VWoodsSaxon3D(x,y,z):
            return self.V0 / (1. + np.exp(self.ainv * (np.sqrt(x*x + y*y + z*z) - self.R)))

        dim = self.ket.mesh.dim
        if dim == 3:
            VWoodsSaxon = VWoodsSaxon3D
        elif dim == 2:
            VWoodsSaxon = VWoodsSaxon2D
        else: # dim == 1:
            VWoodsSaxon = VWoodsSaxon1D

        Vr = self.ket.mesh.apply(VWoodsSaxon)
        self.O_ket.data = Vr * self.O_ket.data
        axes = 'Laplacian' if self.ket.mesh.dim >= 2 else \
               'xx'
        self.O_ket.data += self.hbm * self.ket.derivatives[axes]


