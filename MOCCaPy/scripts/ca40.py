import h5py
import sys
for p in sys.path: print(p)

import mocca
from mocca.edf.param  import Param
from mocca.mesh       import LagrangeMesh
from mocca.mean_field import SlaterDeterminant, HFPsi
from mocca.mean_field.operators import EdfHamiltonian
try:
    from mocca.solve_mfe  import MFESolver, HeavyBall, LinearMix
except ImportError:
    pass

mesh = LagrangeMesh(M=32, d=0.8)
#   dim=3 is default

# Initial mean-field state:
hfpsi = HFPsi(
    n_protons=20, n_neutrons=20,
    n_proton_wf=15, n_neutron_wf=15,
    mesh=mesh,
    init='nilsson', osc_freq=(0.2, 0.2, 0.2),
)
mfs = SlaterDeterminant(hfpsi)

# Energy density functional
param = Param("BSkG1")
bxl = param.create_EDF()
H = EdfHamiltonian(mfs, edf=bxl)

mfe_solver = MFESolver(
    wf0=mfs, edf=bxl,
    energy_tol=1e-9, moment_tol=1e-3,
    spwf_algo=HeavyBall(...),
    scf_algo=LinearMix(...)
)
# Strategy classes do separate responsibilities better.
# It is not the responsibility of the MFESolver to handle the arguments of the strategy classes.

wf, densities, potentials = mfe_solver.solve()

# Alternatively, these may be member functions of the mean-field state object (c.q. SlaterDeterminant or its base
# class), in which case you are just accessing a quantity from the wave function and 'calculate_' can be conveniently
# omitted:
E = wf.energy()
Espwf = wf.spwf_energy(potentials) # alternate way to calculate the energy; potentials are required
Qlm = wf.multipole_moments() # Qlm is a dictionary or list with values for all multipole moments

print(f"{E=}")
print(f"{Espwf=}")
print(f"{E-Espwf=}")

with h5py.File('myfile.hdf5','w') as hdf5file:
    mfe_solver.write_hdf5(hdf5file) # writes the necessary objects
    # other objects can be added
