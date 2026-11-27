# Helper for dumpfile.sh: checks the contents of the debug dumpfile.
# Usage: python3 dumpfile_check.py <dump.hdf5> <mf.hdf5>
import sys
import numpy as np
import h5py

dump_fn, wf_fn = sys.argv[1], sys.argv[2]

with h5py.File(dump_fn, 'r') as d, h5py.File(wf_fn, 'r') as w:
    # a) expected datasets present, with the expected shapes
    expected = {'wavefunctions/states': 3, 'wavefunctions/hpsi': 3,
                'wavefunctions/dpsi': 4, 'wavefunctions/ddpsi': 4,
                'fields/densities/D_I_I': 1, 'fields/potentials/F_I_I': 1}
    for name, ndim in expected.items():
        if name not in d:
            raise SystemExit(f'dumpfile: missing dataset {name}')
        if d[name].ndim != ndim:
            raise SystemExit(f'dumpfile: dataset {name} has rank {d[name].ndim}, expected {ndim}')
    nx = d.attrs['nx']; ny = d.attrs['ny']; nz = d.attrs['nz']
    mv = nx * ny * nz
    for name in expected:
        # densities/potentials are stored as flat 1d arrays of size
        # (mesh points) x (number of components)
        if d[name].size % mv != 0:
            raise SystemExit(f'dumpfile: dataset {name} of size {d[name].size} '
                             f'is not a multiple of the mesh size {mv}')

    # b) <psi|h|psi> from the stored hpsi matches the spenergies of the
    #    converged run for a few spwfs. The mesh volume element
    #    dv = dx^3 * 2^nsym is not stored in the file, so we determine the
    #    common normalization factor and check it against dx^3 * 2^nsym.
    states = d['wavefunctions/states'][:]
    hpsi = d['wavefunctions/hpsi'][:]
    nwt = states.shape[0]
    dx = d.attrs['dx']
    spenergies = w['wavefunctions/hfbasis/spenergies'][:]
    factors = []
    for wave in [0, nwt // 2, nwt - 1]:
        eh = np.sum(states[wave, :, :] * hpsi[wave, :, :])
        if abs(eh) < 1e-300:
            raise SystemExit(f'dumpfile: <psi{wave}|h|psi{wave}> vanishes')
        factors.append(spenergies[wave] / eh)
    factors = np.array(factors)
    # 1e-3 relative: on-the-fly recomputed derivatives differ slightly from
    # the stored ones used for the spenergies (benign reassociation drift).
    if np.max(np.abs(factors - factors[0])) > 1e-3 * abs(factors[0]):
        raise SystemExit(f'dumpfile: <psi|h|psi> normalisation inconsistent '
                         f'across spwfs: {factors}')
    dv = factors[0]
    nsym = round(np.log2(dv / dx**3))
    if abs(dv - dx**3 * 2**nsym) > 1e-5 * dv or not 0 <= nsym <= 3:
        raise SystemExit(f'dumpfile: <psi|h|psi> normalisation {dv} is not '
                         f'dx^3 * 2^nsym with 0 <= nsym <= 3 (dx={dx})')

print('Dumpfile check OK: all datasets present, <psi|h|psi> matches spenergies.')
