import numpy as np
from mocca.mesh.lagrange import LagrangeMesh
from mocca.mean_field.hfpsi import HFPsi

def test_hfpsi_ctor():
    # Input data taken from scripts/Ca40.f90/tant.data
    # Note that nwn, nwp differ from n_neutrons, n_protons
    mesh = LagrangeMesh(M=24, d=1., reduced=True)
    nwn,nwp = 15,15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = HFPsi(n_neutrons=20,
                  n_protons=20,
                  n_proton_wf=nwp, n_neutron_wf=nwn, mesh=mesh,
                  osc_freq=osc_freq
                  )