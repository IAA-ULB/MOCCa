import pytest
import numpy as np

from mocca.mesh import LagrangeMesh
import mocca.mean_field.states as states
import mocca.mean_field.operators as operators
from mocca.mean_field.solvers.dsp import DSP
from mocca.mean_field.states import gramm_schmidt
from mocca.util.timer import Timer

# bad test, the spwfs are all teh same and cannot be orthnormalized.
# def test_dsp():
#     mesh = LagrangeMesh(dim=1, M=6, d=1, reduced=False)
#     nwn, nwp = 15, 15
#     osc_freq = (0.2, 0.2, 0.2)
#     hfpsi = states.HFPsi(
#         n_neutrons=20, n_protons=20,
#         n_proton_wf=nwp, n_neutron_wf=nwn,
#         mesh=mesh, init=None,
#     )
#     hfpsi.data[:,:] = 1.
#     hfpsi.hfblockrange = [(0,8),(8,8),
#                           (8,15),(15,15),
#                           (15,23),(23,23),
#                           (23,30),(30,30)
#                          ]
#     h = operators.hamiltonian.HamiltonianWoodsSaxon(hfpsi)
#     dsp = DSP(h, alpha=.001)
#     dsp.step()
#     dsp.mu = 0.001
#     dsp.step()
#     dsp.step()

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

def test_evolve(check=False, debug=False):
    mesh = LagrangeMesh(M=30, d=.8, reduced=True)
    nwn, nwp = 15, 15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = states.HFPsi(
        n_neutrons=20, n_protons=20,
        n_proton_wf=nwp, n_neutron_wf=nwn,
        mesh=mesh,
        init='nilsson',osc_freq=osc_freq,
        orthogonalize=True, normalize=True
    )
    hamiltonian = operators.hamiltonian.HamiltonianWoodsSaxon(hfpsi)
    dsp = DSP(hamiltonian, alpha=.002)
    # with Timer(name="DSP.evolve") as timer:
    dsp.evolve(nsteps=200, check=check, debug=debug)
    # Timer.report()
    """
iter = 200: h_ii = [
 -39.06333552 -17.64857292 -17.64859173 -17.64865844 -14.91986878 }
 -17.64860771 -17.64860769 -29.07192623 -29.07192623 -29.07192623  } neutronen (20)
  -5.39339831  -5.39327209  -5.39344415  -5.39366001  -3.25459239 }

 -39.06333552 -17.64857292 -17.64859173 -17.64865844 -14.91986865 }
 -17.64860771 -17.64860768 -29.07192623 -29.07192623 -29.07192623  } protonen (20)
  -5.39339387  -5.39326771  -5.39343963  -5.39365522  -3.25464135 }
]

hfblockrange[(0, 7), (7, 7), (7, 15), (15, 15), (15, 22), (22, 22), (22, 30), (30, 30)]

iter = 200: h_ii = [
 -39.06333552 -17.64857292 -17.64859173 -17.64865844 -14.91986878 
 -17.64860771 -17.64860769                                        } block 1 neutrons
                           -29.07192623 -29.07192623 -29.07192623  
  -5.39339831  -5.39327209  -5.39344415  -5.39366001  -3.25459239 } block 3 neutrons

 -39.06333552 -17.64857292 -17.64859173 -17.64865844 -14.91986865 
 -17.64860771 -17.64860768                                        } block 5 protons
                           -29.07192623 -29.07192623 -29.07192623 
  -5.39339387  -5.39326771  -5.39343963  -5.39365522  -3.25464135 } block 7 protons
]

iter = 200: h_ii = [
 -39.06333552 2
 -17.64865844 2
 -17.64860769 2
 -17.64860771 2
 -17.64859173 2
 -17.64857292 
 -14.91986878 
              } block 1 neutrons
 -29.07192623 2
 -29.07192623 2
 -29.07192623 2
  -5.39366001 2 alhoewel in block 1 nog 2 lagere niveau's zijn
  -5.39344415 2 alhoewel in block 1 nog 2 lagere niveau's zijn
  -5.39339831  
  -5.39327209 
  -3.25459239 
              } block 3 neutrons
 -39.06333552 2
 -17.64865844 2
 -17.64860771 2
 -17.64860768 2
 -17.64859173 2
 -17.64857292 
 -14.91986865 
              } block 5 protons
 -29.07192623 2
 -29.07192623 2
 -29.07192623 2
  -5.39365522 2 alhoewel in block 1 nog 2 lagere niveau's zijn
  -5.39343963 2 alhoewel in block 1 nog 2 lagere niveau's zijn
  -5.39339387  
  -5.39326771  
  -3.25464135 
              } block 7 protons
]

iter = 200: h_ii = [
 -39.06333552 2 block 1
 -29.07192623 2 block 3
 -29.07192623 2 block 3
 -29.07192623 2 block 3
 -17.64865844 2 block 1
 -17.64860769 2 block 1
 -17.64860771 2 block 1
 -17.64859173 2 block 1
 -17.64857292 2 block 1
 -14.91986878 2 block 1
  -5.39366001   block 3
  -5.39344415   block 3
  -5.39339831   block 3
  -5.39327209   block 3
  -3.25459239   block 3

 -39.06333552 2 block 5
 -29.07192623 2 block 7
 -29.07192623 2 block 7
 -29.07192623 2 block 7
 -17.64865844 2 block 5
 -17.64860771 2 block 5
 -17.64860768 2 block 5
 -17.64859173 2 block 5
 -17.64857292 2 block 5
 -14.91986865 2 block 5
  -5.39365522   block 7
  -5.39343963   block 7
  -5.39339387   block 7  
  -5.39326771   block 7
  -3.25464135   block 7
]
    """