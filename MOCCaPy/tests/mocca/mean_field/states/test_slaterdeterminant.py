import pytest
import numpy as np
from mocca.mesh.lagrange import LagrangeMesh
from mocca.mean_field import HFPsi, SlaterDeterminant
from mocca.mean_field.states.slaterdeterminant import _sd_occupancies

def test_SlaterDeterminant_ctor():

    mesh = LagrangeMesh(M=24, d=1., reduced=True)
    nwn,nwp = 15,15
    osc_freq = (0.2, 0.2, 0.2)
    mfs = SlaterDeterminant(HFPsi(
        n_neutrons=20,
        n_protons=20,
        n_proton_wf=nwp, n_neutron_wf=nwn, mesh=mesh,
        init='nilsson', osc_freq=osc_freq,
    ))
    mfs.validate()


def test_occupancies():
    n_neutrons = 20
    n_protons  = 20
    nwn = 15
    h_diag = np.array([
        -39.06333552, -17.64857292, -17.64859173, -17.64865844, -14.91986878,
        -17.64860771, -17.64860769, -29.07192623, -29.07192623, -29.07192623,
         -5.39339831,  -5.39327209,  -5.39344415,  -5.39366001,  -3.25459239,
        -39.06333552, -17.64857292, -17.64859173, -17.64865844, -14.91986865,
        -17.64860771, -17.64860768, -29.07192623, -29.07192623, -29.07192623,
         -5.39339387,  -5.39326771,  -5.39343963,  -5.39365522,  -3.25464135,
    ])
    rho,rho2,_,_ = _sd_occupancies(n_neutrons, n_protons, h_diag, nwn)
    print(f"\n{n_neutrons=}: {n_protons=}")
    for i in range(h_diag.size):
        print(f"{i}: {rho[i]} {rho2[i]}")

    n_neutrons = 19
    print(f"\n{n_neutrons=}: {n_protons=}")
    rho,rho2,_,_ = _sd_occupancies(n_neutrons, n_protons, h_diag, nwn)
    for i in range(h_diag.size):
        print(f"{i}: {rho[i]} {rho2[i]}")
