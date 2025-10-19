# A second attempt at MOCCaPy
### opmerkingen van mij beginnen met `###`
# W: antwoorden van mij met # W:
from mocca.param      import Param
from mocca            import MfState # Mstate =  short for mean-field state
from mocca.densities  import DensityVector
from mocca.potentials import PotentialVector
### I strongly recommend honouring PEP 8 (https://peps.python.org/pep-0008/)
###     module names in lowercase
###     class names must start with capital (typically CamelCase)
###     methods start with lowercase
### names should be descriptive
### names of functions should be a verb(+noun) in it because they DO sth (to/with noun)
### names of objects should be nouns because they represent objects
### This helps the user understand what he is importing

# W: ok, I've not yet entirely internalized PEP8; adapted.

from mocca.mesh   import Cartesian3D
### The name Cartesian3D does not cover its meaning: it misses that it is in fact a Lagrange mesh
# W: understood; change depends on whether we implement ONLY Lagrange meshes or wish to also implement others
#             option a: from mocca.mesh import Lagrange    ; mesh = Lagrange('3D') or something along those lines
#             option b; from mocca.mesh import Cartesian3D and just do not support anything else.

#from mocca.EDFs   import BXL
# from mocca.edf import BXL ### module names in lowercase
### overbodig: Param instance genereert automatisch een BXL object op basis van func_file='BXL.func"
### dat is al geimplementeerd
# OK; commented

from mocca.evolve import heavy_ball ### = functie?
#                                       W: Ja. Nieuw naamvoorstel: from mocca.spwf_evolve import heavy_ball_evolve
from mocca.scf    import linear_mix ### = functie?
#                                       W: Ja. Nieuw naamvoorstel: from mocca.scf import linear_mix_potentials
from mocca.io     import write_mf_hdf5

# Define a Cartesian mesh with antiperiodic boundary conditions
mesh = Cartesian3D(n=(32,32,32), d=(0.8,0.8,0.8), bc='antiperiodic' )

# Define a mean-field state
# - that is a Slater determinant (as opposed to a BCS or Bogoliubov state)
# - with Z = 20 and N= 20 particles
# - bound to this specific mesh (ET: i would use 'discretized' rather than 'bound'.
# - initialized with an analytical Nilsson s.p. hamiltonian (alternative is random initialisation)
# - (by default) maximally symmetric
# - in a modelspace defined by 40 neutron and 40 proton spwfs                     NOTE: should modelspace be separated as a concept?
wf = MfState(Z=20, N=20, iniwfs='Nilsson', mesh=mesh, type='Slater', nwn=40, nwp=40)
### Het Nucleus concept is afgevoerd?
#   W: MfState vervangt het Nucleus concept. "Nucleus" was (in mijn hoofd) altijd een container voor de single-particle wavefunctions.
#      Expliciet voor MfState kiezen verduidelijkt wat dit object bevat: een boel single-particle wavefunctions + de noodzakelijke geassocieerde matrices (occupations etc...)
#      die samen een mean-field many-body golffunctie vormen. Naast duidelijkheid is het ook een voordeel dat MfState duidelijk toepasbaar is voor (a) Pasta en (b) toekomstige developments
#      richting beyond-mean-field technieken, waarin je vaak verschillende MfStates superposeert om een kern te beschrijven.
#
### mfstate is een functie? of een klasse waarvan hier de constructor wordt opgeroepen?  W: constructor zou ik denken.
### noch mfstate, noch wf worden hieronder gebruikt?                                     W: oeps, wf zou een argument voor evolve moeten zijn geweest!
### Mij lijkt het dat een (mean field) state een instance van een klasse is vermits state een noun is,
### Misschien moet je me even uitleggen wat het concept mean-field state juist inhoudt.
### ik vind de namen verwarrend. Is een mfstate een golf-functie?
### type='Slater' suggereert dat het misschien logischer zou zijn om iets te doen in de aard van
wf = SlaterDeterminant(Z=20, N=20, iniwfs='Nilsson', mesh=mesh, nwp=40, nwn=40)
# W: dit zou ook werken.
#    Drie verschillende opties moeten we dan wel voorzien
#    * SlaterDeterminant
#    * BCSState  (minst belangrijk van de drie)
#    * BogoliubovState
# Ik gaf een kleine uitleg van "State" hierboven, maar kan beter doen als dat nodig is.

# Define an EDF ....
param = Param("BSkG1") # read BSkG1.param file from mocca/parametrezations (a soft link to tantalus_full/parametrizations
# and create a BXL object:
bxl = param.create_EDF() ### the name bxl is arbitrary, but hints at it being a BXL objest which is selected by
                         ### func_file=BXL.func in BSkG1.param
                         ### param is accessible from bxl as bxl.param

# Note, I really have no idea on how an "EDF" or a "parameterisation" should be stored.
#  - in external files in some home-made format like in MOCCa?
#  - "hardcoded" inside the Python routines?
#  - something else?
### Currently the .param files are read, storing the parameters as attributes
### A different representation may be useful when the Fortran executable is phased out in the future.

# set some convergence_criteria
### misschien ook beter een klasse voor maken. de constructor is een handige plaats voor default values
# W: akkoord.
cc = {
        'energy_precision' : 1e-9  # in units of MeV
        'moment_precision' : 1e-3  # dimensionless
        # ..... + others
        # those that are not specified explicitly should be set to sensible defaults
     }
# ... and then perform the evolution until convergence is reached; only at this point is an EDF actually needed
### evolve is not known... (as are none of the methods used below)
### wat doet evolve precies? Het is een verb dat (in mijn gedachten) een tijdsevolutie berekend...
### names should be descriptive - al je huidige gebruikers weten ongetwijfeld wat daarmee bedoeld
### wordt, maar nieuwe gebruikers worden beter niet op het verkeerde been gezet. Een van de basisprincipes van Python
### is het principle of least surprise. Hoe beter een naam aangeeft wat er gebeurt, hoe minder verrassingen en hoe
### minder steil de leercurve is.
#
# W: ok, akkoord. "evolve" is inderdaad een anachronisme; wat ik bedoel is simpelweg "start de iteraties".
#                 "iterate" is daarentegen onduidelijk omdat iemand zonder voorkennis zou denken aan
#                 "doe 1 iteratie". "solve" is daartegen ook zo vaag. Eventueel "solve_mean_field_equations"?
#                 Andere voorstellen zeker welomen.
wf, densities, potentials = evolve(wf, strategy=heavy_ball, scf=linear_mix, edf=bxl, convergence_criteria = cc)  # W: dit voorstel is beter, ik voegde ook wf toe als initiële startpunt
### strategy en scf zijn misschien iets te vage namen.
#               |> this chooses a strategy to evolve the single-particle wavefunctions
#               |> this chooses a strategy to evolve the single-particle potentials
#               |> this chooses the EDF form to use
#               |> this chooses the parameterisation ### param is an attribute of bxi now
#
# Additional optional arguments to evolve could be
#   - R0: a density-vector to start from; if not provided, it should be one calculated from wf
#   - F0: a potential-vector to start from; if not provided, it should be calculated from wf
### moet wf hier ook geen argument zijn? W: Ja! simpelweg vergeten
### wf (van wf = mfstate(...) hierboven) wordt hier overschreven? maar werd nergens gebruikt

# After the evolution
# - wf is the final solution
# - "densities"  = the vector of densities that corresponds to it; an instance of the densityvector class
# - "potentials" =  the vector of potentials that corresponds to it; an instance of the potentialvector class
### wf vind ik acceptable als naam, R en F niet. om te ze zijn niet lowercase en niet descriptief
# W: Ok; veranderd naar densities, potentials
# W: Note that many iterative solution strategies require densities and potentials that are NOT linked to the current MfState; hence the need to separate these concepts.
#
# There is however, an unsolved problem here: the precise realisation of DensityVector and PotentialVector depend on the EDF selected.
### dat is het leuken van heet dynamische karakter van Python. Objecten kunnen veranderen. Zolang de EDF in kwestie weet
### hoe DensityVector/PotentialVector gemaakt moet worden is er geen probleem. De verantwoordelijkheid van EDFs wordt
### daarmee wel groter maar het is logisch dat als B afhangt van A, dat A dan ook weet hoe, en daar verantwoordelijk
### voor is.
### Wat bedoel je precies met een vector van densities/potentials hier? zijn dat scalaire functies gediscretiseerd op 
### het mesh? zijn dat  echte vectoren of gaat het eerder om een collectie van densities/potentials?
#
# W: 'vectoren' zijn collecties van dichtheden en potentialen; elke dichtheid is een functie op het mesh, mogelijk met indices.
#
# De allersimpelste densityvector die we gebruiken is bvb
#
#           R = ( \rho(r), \tau(r), J(r))
#
# waarbij
#   \rho(r) = de gewone nucleon dichtheid, met een neutron en proton component
#   \tau(r) = de kinetische dichtheid, met neutron en proton component
#   J(r)    = de spin-kinetische dichtheid, dat is een rank-2 tensor (i.e. 27 componenten!) met een neutron en proton component.
#
# In de FORTRAN code worden deze gedefinieerd als
#
#   \rho = array van dimensie (nx*ny*nz,4)
#   \tau = array van dimensie (nx*ny*nz,4)
#   J    = array van dimensie (nx*ny*nz,27,4)
#
# De laatste index gaat van 1 - 4, met de conventie
#   1 -> neutron
#   2 -> proton
#   3 -> isospin 0, ofte neutron + proton, som van component 1 & 2
#   4 -> isospin 1, ofte neutron - proton, verschil van component 1 & 2
#


# Calculate some observables and print them
# NOTE: I'm unsure if we should make the effort to save some CPU cycles and not recalculate things that are likely to have
# been calculated during the call to evolve
E = calculate_energy(wf,EDF,param)   # here the EDF and parameterisation are again needed;
                                   # you might conceivably want to calculate the energy wf with another model
E = calculate_energy(wf,bxl)   ### param is an attribute of bxl
Espwf = calculate_spwf_energy(wf, F, EDF,param) # alternate way to calculate the energy; potentials are required
Espwf = calculate_spwf_energy(wf, F, bxl) ### param is an attribute of bxl
Qlm = calculate_multipole_moments(wf) # Qlm is a dictionary or list with values for all multipole moments
# ..... + many others + .....
### die functies moeten nog geimporteerd worden.
print('Total energy (from functional)', E)
print('Total energy (from spwfs)', Espwf)
print('Difference', E-Espwf)     # Should be very small when converged

# ... and write the mean-field state to file
write_mf_hdf5(wf,R,F0)
### hdf5file openen en elk van de object apart wegschrijven
# W: dat is inderdaad duidelijker, maar leidt dat niet tot moeilijkheden voor het uitwisselen van files door de mindere standardisatie?
import h5py
hdf5file = h5py.File('myfile.hdf5','w')
wf.write_hdf5(hdf5file)
R.write_hdf5(hdf5file)
F0.write_hdf5(hdf5file)

# NOTE
#  - this function should take care of meta-data related to wf
#  - I'm undecided whether observables should be stored
#  - the densityvector and potentialvector should be optional arguments