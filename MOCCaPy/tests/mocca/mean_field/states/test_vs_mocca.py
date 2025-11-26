import pytest
from pathlib import Path
import numpy as np
import h5py

from mocca.mean_field.states import HFPsi
from mocca.mesh import LagrangeMesh

def test_HFPsi():
    mesh = LagrangeMesh(M=30, d=.8, reduced=True)
    n_neutrons, n_protons = 20, 20
    nwn, nwp = 15, 15
    osc_freq = (0.2, 0.2, 0.2)
    hfpsi = HFPsi(
        n_neutrons=n_neutrons, n_protons=n_protons,
        n_proton_wf=nwp, n_neutron_wf=nwn,
        mesh=mesh,
        init='nilsson',osc_freq=osc_freq,
        # orthogonalize=True, normalize=True
    )
    wfs = HFPsi(
        n_neutrons=n_neutrons, n_protons=n_protons,
        n_proton_wf=nwp, n_neutron_wf=nwn,
        mesh=mesh,
        init='nilsson',osc_freq=osc_freq,
        # orthogonalize=True, normalize=True
    )

    nwt = nwn + nwp
    mocca_output = Path(__file__).parent / "mocca_output"
    x = np.empty(mesh.linear_size, dtype=np.float64)
    y = np.empty(mesh.linear_size, dtype=np.float64)
    z = np.empty(mesh.linear_size, dtype=np.float64)

    # Read wfs and mesh points
    with open(mocca_output/"wfs.txt", "r") as wfs_txt:
        lines = wfs_txt.readlines()
        for l,line in enumerate(lines):
            words = line.split()
            x[l] = float(words[0])
            y[l] = float(words[1])
            z[l] = float(words[2])
            for iw,w in enumerate(words[3:]):
                wfs.data[l,iw] = float(w)
                if iw % 4 ==0:
                    w0 = wfs.data[l,iw]
                else:
                    assert wfs.data[l,iw] == pytest.approx(w0)
                print(f"[{l}, {iw}] {hfpsi.data[l, iw]} ?= {wfs.data[l, iw]}")
                pass

    # compare mocca_output to mesh
    gridx = mesh.cast2linear(mesh.gridx)
    gridy = mesh.cast2linear(mesh.gridy)
    gridz = mesh.cast2linear(mesh.gridz)

    for ig in range(mesh.linear_size):
        assert gridx[ig] == pytest.approx(x[ig])
        assert gridy[ig] == pytest.approx(y[ig])
        assert gridz[ig] == pytest.approx(z[ig])
    # If arriving here, the grid points are in the expected order, so we can proceed with the

    for ig in range(mesh.linear_size):
        for iw in range(nwt):
            print(f"[{ig}, {iw}] {hfpsi.data[ig, iw]} ?= {wfs.data[ig, iw]}")
            # assert hfpsi.data[ig,iw] == pytest.approx(wfs.data[ig,iw], rel=1e-6, abs=1e-6)

    # Write to hdf5 file
    f5_path = mocca_output/"test_HFPsi.hdf5"
    with h5py.File(f5_path, "w") as f5:
        f5.create_dataset("x", data=x)
        f5.create_dataset("y", data=y)
        f5.create_dataset("z", data=z)
        f5.create_dataset("wfs", data=wfs.data)
        for name in f5:
            print(f"{f5_path.name} {f5.name} {name}")

def test_D_matrices():
    mocca_output = Path(__file__).parent / "mocca_output"

    N = 15
    D1p = np.empty(N*N, dtype=np.float64, order='F')
    D1m = np.empty(N*N, dtype=np.float64, order='F')
    D2p = np.empty(N*N, dtype=np.float64, order='F')
    D2m = np.empty(N*N, dtype=np.float64, order='F')

    f5_path = mocca_output/"test_D_matrices.hdf5"
    with h5py.File(f5_path, "w") as f5:
        for txt in ['lagx.txt', 'lagy.txt', 'lagz.txt' ]:
            with open(mocca_output/txt, "r") as D_txt:
                lines = D_txt.readlines()
                for l,line in enumerate(lines):
                    words = line.split()
                    D1p[l] = float(words[0])
                    D1m[l] = float(words[1])
                    D2p[l] = float(words[2])
                    D2m[l] = float(words[3])

            f5.create_dataset(f"D1p{txt[3]}", data=D1p.reshape((N, N),order='F'))
            f5.create_dataset(f"D1m{txt[3]}", data=D1m.reshape((N, N),order='F'))
            f5.create_dataset(f"D2p{txt[3]}", data=D2p.reshape((N, N),order='F'))
            f5.create_dataset(f"D2m{txt[3]}", data=D2m.reshape((N, N),order='F'))

        for name in f5:
            print(f"{f5_path.name} {f5.name} {name}")

        assert f5['D1px'] == f5['D1py']
