import pytest
from pathlib import Path
import numpy as np
import h5py

from mocca.mean_field.states import HFPsi
from mocca.mesh import LagrangeMesh

test_folders = [
    "spwf_0.2_0.2_0.2",
    "spwf_0.2_0.2_0.16",
]

@pytest.mark.parametrize("test_folder", test_folders)
def test_HFPsi(test_folder):
    mesh = LagrangeMesh(M=30, d=.8, reduced=True)
    n_neutrons, n_protons = 20, 20
    nwn, nwp = 15, 15
    osc_freq = [float(w) for w in test_folder.split('_')[1:]]
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
    mocca_output = Path(__file__).parent / test_folder
    x = np.empty(mesh.linear_size, dtype=np.float64)
    y = np.empty(mesh.linear_size, dtype=np.float64)
    z = np.empty(mesh.linear_size, dtype=np.float64)

    rel = 1e-5
    abs = 1e-4
    # Read wfs and mesh points
    txt = 'nilsson_wfs.txt'
    with open(mocca_output/txt, 'r') as wfs_txt:
        lines = wfs_txt.readlines()
        is_same = 0
        for l,line in enumerate(lines):
            words = line.split()
            x[l] = float(words[0])
            y[l] = float(words[1])
            z[l] = float(words[2])
            s =f'{l:4}'
            for iw,w in enumerate(words[3:]):
                if iw % 4 == 0:
                    s += ' '
                wfs.data[l,iw] = float(w)
                if hfpsi.data[l, iw] == pytest.approx(wfs.data[l, iw], rel=rel, abs=abs):
                    is_same += 1
                    s += '.'
                else:
                    s += 'F'

                # else:
                #     print(f"[{l}, {iw}] {hfpsi.data[l, iw]} ?= {wfs.data[l, iw]} {hfpsi.data[l, iw] == pytest.approx(wfs.data[l, iw], rel=1e-6, abs=abs)}")
                # assert hfpsi.data[l, iw] == pytest.approx(wfs.data[l, iw], rel=1e-6, abs=abs), \
                #     f"[{l}, {iw}] {hfpsi.data[l, iw]} ?= {wfs.data[l, iw]} {hfpsi.data[l, iw] == pytest.approx(wfs.data[l, iw], rel=1e-6, abs=1e-6)}"
            print(s)
        print(f"{is_same=}")
        pass

    same = []
    for iw in range(nwt):
        spwf = hfpsi.d3[:,:,iw]
        for jw in range(nwt):
            if jw in same:
                continue
            if np.all(spwf == pytest.approx(wfs.d3[:,:,jw], rel=rel, abs=abs)):
                same.append(jw)
                break
            elif np.all(spwf == pytest.approx(-wfs.d3[:, :, jw], rel=rel, abs=abs)):
                same.append(-jw)
                break
        else:
            same.append('?')

    print(same)
    if test_folder == "spwf_0.2_0.2_0.16":
        assert '?' not in same

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
    f5_path = mocca_output/(txt.replace('.txt','.hdf5'))
    with h5py.File(f5_path, "w") as f5:
        f5.create_dataset("x", data=x)
        f5.create_dataset("y", data=y)
        f5.create_dataset("z", data=z)
        f5.create_dataset("wfs", data=wfs.data)
        for name in f5:
            print(f"{f5_path.name} {f5.name} {name}")


@pytest.mark.parametrize("test_folder", test_folders)
def test_D_matrices(test_folder):
    mocca_output = Path(__file__).parent / test_folder

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

        # print(f5['D1px'][:])
        # print(f5['D1py'])
        assert np.all(f5['D1px'][:,:] == f5['D1py'][:,:])
        assert np.all(f5['D1px'][:,:] == f5['D1pz'][:,:])
        assert np.all(f5['D1mx'][:,:] == f5['D1my'][:,:])
        assert np.all(f5['D1mx'][:,:] == f5['D1mz'][:,:])

        mesh = LagrangeMesh(M=30, d=.8, reduced=True)
        DpEx = mesh.D[0,0]
        DmEx = mesh.DmE[0,0]
        print(DpEx[:4, :4])
        # print(DmEx[:4, :4])
        # print(f5['D1px'][:4,:4])
        print(f5['D1mx'][:4,:4])
        rel=1e-5
        abs=1e-5
        assert (DpEx[:, :] == pytest.approx(-f5['D1mx'][:, :],rel=rel,abs=abs))
        assert (DmEx[:, :] == pytest.approx(-f5['D1px'][:, :],rel=rel,abs=abs))
        pass

@pytest.mark.parametrize("test_folder", test_folders)
def test_derivatives(test_folder):
    mocca_output = Path(__file__).parent / test_folder
    hfpsi_h5 = mocca_output/"nilsson_wfs.hdf5"
    with h5py.File(hfpsi_h5, "r") as f5:
        wfs = f5['wfs']
        # shape = wfs.shape
        hfpsi_data = np.array(f5['wfs'][:,:])

    mesh = LagrangeMesh(M=30, d=.8, reduced=True)

    n_neutrons, n_protons = 20, 20
    nwn, nwp = 15, 15
    osc_freq = [float(w) for w in test_folder.split('_')[1:]]
    hfpsi = HFPsi(
        n_neutrons=n_neutrons, n_protons=n_protons,
        n_proton_wf=nwp, n_neutron_wf=nwn,
        mesh=mesh,
        init='nilsson', osc_freq=osc_freq,
        # orthogonalize=True, normalize=True
    )
    # store the MOCCa hfpsi in the MOCCaPy HFPsi
    hfpsi.data[:,:] = hfpsi_data


    f5_path = mocca_output/"test_derivatives.hdf5"
    with h5py.File(f5_path, "w") as f5:
        for axes in ['x','y','z', 'Laplacian']:
            dhfpsi_du = hfpsi.differentiate(axes)
            dwfs_du = np.empty_like(dhfpsi_du)
            txt = f'nilsson_nabla_{axes}_wfs.txt' if len(axes) == 1 else f'nilsson_laplacian_wfs.txt'
            with open(mocca_output / txt, 'r') as wfs_txt:
                lines = wfs_txt.readlines()
                for l, line in enumerate(lines):
                    words = line.split()
                    for iw,word in enumerate(words[3:]):
                        dwfs_du[l,iw] = float(word)

                f5.create_dataset(f'dwfs_d{axes}', data=dwfs_du)

            for ig in range(mesh.linear_size):
                print(f"{axes} {ig}")
                for iw in range(4*hfpsi.n_total_wf):
                    # assert dhfpsi_du[ig,iw] == dwfs_du[ig,iw]
                    if not dhfpsi_du[ig,iw] == pytest.approx(dwfs_du[ig,iw], rel=1e-5, abs=1e-3):
                        print(f"({ig},{iw}) {dhfpsi_du[ig,iw]} != {dwfs_du[ig,iw]}")

@pytest.mark.parametrize("test_folder", test_folders)
def test_create_hdf5(test_folder):
    npoints = 15**3
    nwt = 30
    nq = 4*nwt
    p_h5 = Path(__file__).parent / (test_folder + '.h5')

    if p_h5.exists():
        print(f"File {p_h5} already exists")
        return

    p_test_folder = Path(__file__).parent / test_folder
    with h5py.File(p_h5, "w") as f5:
        print(f"Creating {p_h5}")
        for p in p_test_folder.glob('*.txt'):
            shape = (npoints, nq) if ('wfs' in p.name) else \
                    (nwt * nwt, 4)
            data = np.empty(shape, dtype=np.float64, order='F')
            with open(p, 'r') as f:
                lines = f.readlines()
                for l, line in enumerate(lines):
                    words = line.split()
                    if 'wfs' in p.name:
                        for iq, word in enumerate(words[3:]):
                            data[l,iq] = float(word)
                    else:
                        for iq, word in enumerate(words):
                            data[l,iq] = float(word)
            print(f"Adding dataset {p.stem}")
            f5.create_dataset(p.stem, data=data)