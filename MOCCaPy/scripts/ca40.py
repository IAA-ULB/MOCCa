import h5py

from mocca.param     import Param
from mocca.meanfield import SlaterDeterminant
from mocca.mesh      import LagrangeMesh
# from mocca.observables import calculate_energy, calculate_spwf_energy, calculate_multipole_moments
from mocca.solve_mfe import MFESolver

mesh = LagrangeMesh(n=32, d=0.8, bc='antiperiodic')
# dim=3 is default
wf = SlaterDeterminant(Z=20, N=20, wf_init='Nilsson', mesh=mesh, nwp=40, nwn=40)
param = Param("BSkG1")
bxl = param.create_EDF()

mfe_solver = MFESolver(energy_tol=1e-9, moment_tol=1e-3, wf0=wf, spwf_algo='heavy_ball', scf_algo='linear_mix', edf=bxl)
# note that we opted to pass strings for strategy and scf. Now we do not have to import the functions heavy_ball and
# linear_mix. They are only used internally, a registry for mapping strategy and scf strings to the corresponding
# functions may help to extend the possibilities without modifying the MFESolver class.
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
