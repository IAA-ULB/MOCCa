import pytest
import numpy as np

from mocca.mesh import LagrangeMesh
import mocca.mean_field.states as states
from mocca.mean_field.solvers.dsp import DSP

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
    h = states.HamiltonianWoodsSaxon2(hfpsi)
    dsp = DSP(h, alpha=.001)
    dsp.step()

