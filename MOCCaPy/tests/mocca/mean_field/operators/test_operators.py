import numpy as np
import pytest

from mocca.mesh import LagrangeMesh
from mocca.mean_field.states import HFPsi
from mocca.mean_field.operators import *

def test_Operator():

    mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=False)
    nwn, nwp = 15, 15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = HFPsi(n_neutrons=20, n_protons=20,
                  n_proton_wf=nwp, n_neutron_wf=nwn,
                  mesh=mesh, init=None,
                  )
    hfpsi.data[:,:] = 1.
    hfpsi.hfblockrange = [(0,8),(8,8),
                          (8,15),(15,15),
                          (15,23),(23,23),
                          (23,30),(30,30)
                         ]
    Oij = Operator(hfpsi)
    with pytest.raises(AttributeError):
        Oij.compute_matrix_representation()

def test_Overlap():

    mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=False)
    nwn, nwp = 15, 15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = HFPsi(n_neutrons=20, n_protons=20,
                  n_proton_wf=nwp, n_neutron_wf=nwn,
                  mesh=mesh, init=None,
                  )
    hfpsi.data[:,:] = 1.
    hfpsi.hfblockrange = [(0,8),(8,8),
                          (8,15),(15,15),
                          (15,23),(23,23),
                          (23,30),(30,30)
                         ]
    overlap = Overlap(hfpsi)
    overlap.compute_matrix_representation()
    for block in overlap.matrix:
        assert np.all(block == 24.0)

def test_HamiltonianWoodsSaxon():

    mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=False)
    nwn, nwp = 15, 15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = HFPsi(n_neutrons=20, n_protons=20,
                  n_proton_wf=nwp, n_neutron_wf=nwn,
                  mesh=mesh, init=None,
                  )
    hfpsi.data[:,:] = 1.
    hfpsi.hfblockrange = [(0,8),(8,8),
                          (8,15),(15,15),
                          (15,23),(23,23),
                          (23,30),(30,30)
                         ]
    hws = hamiltonian.HamiltonianWoodsSaxon(hfpsi)
    hws.compute_matrix_representation()
    # for block in hws.matrix:
    #     assert np.all(block == 24.0)

def test_HamiltonianWoodsSaxon2():

    mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=False)
    nwn, nwp = 15, 15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = HFPsi(n_neutrons=20, n_protons=20,
                  n_proton_wf=nwp, n_neutron_wf=nwn,
                  mesh=mesh, init=None,
                  )
    hfpsi.data[:,:] = 1.
    hfpsi.hfblockrange = [(0,8),(8,8),
                          (8,15),(15,15),
                          (15,23),(23,23),
                          (23,30),(30,30)
                         ]
    hws = hamiltonian.HamiltonianWoodsSaxon2(hfpsi)
    hws.compute_matrix_representation()
    # for block in hws.matrix:
    #     assert np.all(block == 24.0)