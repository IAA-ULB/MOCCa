import pytest
import numpy as np

from mocca.mesh import LagrangeMesh
import mocca.mean_field.states as states
import mocca.mean_field.operators as operators
from mocca.mean_field.solvers.dsp import DSP, gramm_schmidt
from mocca.util.timer import Timer

def test_dsp():
    mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=False)
    nwn, nwp = 15, 15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = states.HFPsi(
        n_neutrons=20, n_protons=20,
        n_proton_wf=nwp, n_neutron_wf=nwn,
        mesh=mesh, init=None,
    )
    hfpsi.data[:,:] = 1.
    hfpsi.hfblockrange = [(0,8),(8,8),
                          (8,15),(15,15),
                          (15,23),(23,23),
                          (23,30),(30,30)
                         ]
    h = operators.hamiltonian.HamiltonianWoodsSaxon2(hfpsi)
    dsp = DSP(h, alpha=.001)
    dsp.step()
    dsp.mu = 0.001
    dsp.step()
    dsp.step()

def test_gramm_schmidt():
    mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=True)
    nwn, nwp = 15, 15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = states.HFPsi(
        n_neutrons=20, n_protons=20,
        n_proton_wf=nwp, n_neutron_wf=nwn,
        mesh=mesh,
        init=None,
    )
    hfpsi.data[:, :] = np.random.random(hfpsi.data.shape)
    hfpsi.hfblockrange = [(0, 8), (8, 8),
                          (8, 15), (15, 15),
                          (15, 23), (23, 23),
                          (23, 30), (30, 30)
                          ]
    overlap = operators.Overlap(hfpsi)
    hfpsi.sp_energies = overlap.compute_diagonal_elements()
    gramm_schmidt(hfpsi, order=hfpsi.sp_energies, normalize=False)

    overlap_matrix = overlap.compute_matrix_representation()
    for ib in range(8):
        block_ib = overlap_matrix[ib]
        n = block_ib.shape[0]
        for i in range(n):
            for j in range(n):
                print(f"({i},{j}) -> {block_ib[i,j]}")
                if i != j:
                    assert block_ib[i,j] == pytest.approx(0.0)
    hfpsi.normalize()
    overlap_matrix = overlap.compute_matrix_representation()
    for ib in range(8):
        block_ib = overlap_matrix[ib]
        n = block_ib.shape[0]
        for i in range(n):
            for j in range(n):
                print(f"{ib=} ({i},{j}) -> {block_ib[i,j]}")
                if i != j:
                    assert block_ib[i,j] == pytest.approx(0.0)
                else:
                    assert block_ib[i,j] == pytest.approx(1.0)

def test_evolve():
    mesh = LagrangeMesh(M=30, d=.8, reduced=False)
    nwn, nwp = 15, 15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = states.HFPsi(
        n_neutrons=20, n_protons=20,
        n_proton_wf=nwp, n_neutron_wf=nwn,
        mesh=mesh,
        init='nilsson',osc_freq=osc_freq,
        orthogonalize=True, normalize=True
    )
    overlap = operators.Overlap(hfpsi)
    # overlap_matrix = overlap.compute_matrix_representation()
    hamiltonian = operators.hamiltonian.HamiltonianWoodsSaxon2(hfpsi)
    dsp = DSP(hamiltonian, alpha=.002)
    with Timer(name="DSP.evolve") as timer:
        dsp.evolve(nsteps=10)
    Timer.report()