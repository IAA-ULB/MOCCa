import pytest
import numpy as np
from mocca.mesh.lagrange import LagrangeMesh
from mocca.mean_field import HFPsi, SlaterDeterminant

def test_slaterdeterminant_ctor():

    mesh = LagrangeMesh(M=24, d=1., reduced=True)
    nwn,nwp = 15,15
    osc_freq = (0.2, 0.2, 0.2)
    sd = SlaterDeterminant(HFPsi(
        n_neutrons=20,
        n_protons=20,
        n_proton_wf=nwp, n_neutron_wf=nwn, mesh=mesh,
        init='nilsson', osc_freq=osc_freq,
    ))
    pass