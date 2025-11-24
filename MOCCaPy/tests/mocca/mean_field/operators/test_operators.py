import numpy as np
import pytest

from mocca.mesh import LagrangeMesh
from mocca.mean_field.states import HFPsi
from mocca.mean_field.operators import *
from MOCCaPy.mocca.mean_field.operators.hamiltonian import *

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

# def test_HamiltonianWoodsSaxon():
#
#     mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=False)
#     nwn, nwp = 15, 15
#     osc_freq = (0.2, 0.2, 0.2)
#     hfpsi = HFPsi(n_neutrons=20, n_protons=20,
#                   n_proton_wf=nwp, n_neutron_wf=nwn,
#                   mesh=mesh, init=None,
#                   )
#     hfpsi.data[:,:] = 1.
#     hfpsi.hfblockrange = [(0,8),(8,8),
#                           (8,15),(15,15),
#                           (15,23),(23,23),
#                           (23,30),(30,30)
#                          ]
#     hws = hamiltonian.HamiltonianWoodsSaxon(hfpsi)
#     hws.compute_matrix_representation()
#     # for block in hws.matrix:
#     #     assert np.all(block == 24.0)

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

def test_V_WoodsSaxon():
    for dim in range(1,4):
        mesh = LagrangeMesh(dim=dim, M=30, d=.8, reduced=True)
        V_WoodsSaxon = V_WoodsSaxon3D if mesh.dim == 3 else \
                       V_WoodsSaxon2D if mesh.dim == 2 else \
                       V_WoodsSaxon1D
        nwn, nwp = 15, 15
        osc_freq = (0.2, 0.2, 0.2)
        hfpsi = HFPsi(n_neutrons=20, n_protons=20,
                      n_proton_wf=nwp, n_neutron_wf=nwn,
                      mesh=mesh, init=None,
                      )
        ham = HamiltonianWoodsSaxon(hfpsi)
        Vr = mesh.apply(V_WoodsSaxon, ham.ws_parms)
        assert np.all(Vr > ham.ws_parms['V0'])
        assert np.all(Vr <= .0)
        pass


