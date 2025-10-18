# A second attempt at MOCCaPy
from mocca.param      import Param
### I strongly recommend honouring PEP 8 (https://peps.python.org/pep-0008/)
from mocca            import mfstate # mfstate =  short for mean-field state
from mocca.densities  import densityvector
from mocca.potentials import potentialvector
### class names must start with capital (typically CamelCase)
### methods start with lowercase

from mocca.mesh   import Cartesian3D
### The name Cartesian3D does not cover its meaning: it misses that it is in fact a Lagrange mesh

from mocca.EDFs   import BXL
# from mocca.edf import BXL ### module names in lowercase
### overbodid: Param objegt genereert BXL automatisch

from mocca.evolve import heavy_ball
from mocca.scf    import linear_mix
from mocca.io     import write_mf_hdf5

# Define a Cartesian mesh with antiperiodic boundary conditions
mesh = Cartesian3D(n=(16,16,16), d=(0.8,0.8,0.8), bc='antiperiodic' )

# Define a mean-field state
# - that is a Slater determinant (as opposed to a BCS or Bogoliubov state)
# - with Z = 20 and N= 20 particles
# - bound to this specific mesh (ET: i would use 'discretized' rather than 'bound'.
# - initialized with an analytical Nilsson s.p. hamiltonian (alternative is random initialisation)
# - (by default) maximally symmetric
# - in a modelspace defined by 40 neutron and 40 proton spwfs                     NOTE: should modelspace be separated as a concept?
wf = mfstate(Z=20,N=20,iniwfs='Nilsson', mesh=mesh, type='Slater',nwn=40,nwp=40)

# Define an EDF ....
param = Param("BSkG1") # read BSkG1.param
# and create a BXL object:
bxl = param.create_EDF() ### the name bxl is arbitrary, but hints at it being a BXL objest

# Note, I really have no idea on how an "EDF" or a "parameterisation" should be stored.
#  - in external files in some home-made format like in MOCCa?
#  - "hardcoded" inside the Python routines?
#  - something else?
# Currently the .param files are read, storing the parameters as attributes

# set some convergence_criteria
cc = {
        'energy_precision' : 1e-9  # in units of MeV
        'moment_precision' : 1e-3  # dimensionless
        # ..... + others
        # those that are not specified explicitly should be set to sensible defaults
     }
# ... and then perform the evolution until convergence is reached; only at this point is an EDF actually needed
### OOPS evolve is not knowne... (as are none of the methods used belowe
wf, R, F = evolve(strategy=heavy_ball, scf=linear_mix, EDF, param, convergence_criteria = cc)
#               |> this chooses a strategy to evolve the single-particle wavefunctions
#                                      |> this chooses a strategy to evolve the single-particle potentials
#                                                      |> this chooses the EDF form to use
#                                                           |> this chooses the parameterisation
#
# Additional optional arguments to evolve could be
#   - R0: a density-vector to start from; if not provided, it should be one calculated from wf
#   - F0: a potential-vector to start from; if not provided, it should be calculated from wf

# After the evolution
# - wf is the final solution
# - R the vector of densities that corresponds to it; an instance of the densityvector class
# - F the vector of potentials that corresponds to it; an instance of the potentialvector class
#
# There is however, an unsolved problem here: the precise realisation of DensityVector and PotentialVector depend on the EDF selected.

# Calculate some observables and print them
# NOTE: I'm unsure if we should make the effort to save some CPU cycles and not recalculate things that are likely to have
# been calculated during the call to evolve
E = calculate_energy(wf,EDF,param)   # here the EDF and parameterisation are again needed;
                                   # you might conceivably want to calculate the energy wf with another model
Espwf = calculate_spwf_energy(wf,F,EDF,param) # alternate way to calculate the energy; potentials are required
Qlm = calculate_multipole_moments(wf) # Qlm is a dictionary or list with values for all multipole moments
# ..... + many others + .....

print *, 'Total energy (from functional)', E
print *, 'Total energy (from spwfs)', Espwf
print *, 'Difference', E-Espwf     # Should be very small when converged


# ... and write the mean-field state to file
write_mf_hdf5(wf,R,F0)
# NOTE
#  - this function should take care of meta-data related to wf
#  - I'm undecided whether observables should be stored
#  - the densityvector and potentialvector should be optional arguments
