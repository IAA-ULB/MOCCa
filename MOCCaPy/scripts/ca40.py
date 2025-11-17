import h5py


from mocca.param      import Param
from mocca.mesh       import LagrangeMesh
from mocca.mean_field import SlaterDeterminant
try:
    from mocca.solve_mfe  import MFESolver
except ImportError:
    pass

mesh = LagrangeMesh(M=32, d=0.8)
#   dim=3 is default

# Initial mean-field state:
wf0 = SlaterDeterminant(
    n_protons=20, n_neutrons=20,
    n_proton_wf=15, n_neutron_wf=15,
    mesh=mesh,
    init='nilsson', osc_freq=(0.2, 0.2, 0.2),
)

# Energy density functional
param = Param("BSkG1")
bxl = param.create_EDF()

mfe_solver = MFESolver(energy_tol=1e-9, moment_tol=1e-3, wf0=wf0, spwf_algo='heavy_ball', scf_algo='linear_mix', edf=bxl)
# note that we opted to pass strings for strategy and scf. Now we do not have to import the functions heavy_ball and
# linear_mix. They are only used internally, a registry for mapping strategy and scf strings to the corresponding
# functions (or objects, if they need state) may help to extend the possibilities without modifying the MFESolver class.
# Strategy classes would separate responsibilities better. It is not the responsibility of the MFESolver to handle the
# arguments of the strategy classes.

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
