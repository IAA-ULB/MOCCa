# MOCCa HDF5 files

## Overall structure
    
    Groups
     /                   : general information
     |
     | > /wavefunctions/ : single-particle wavefunctions and their properties 
     | | > /wavefunctions/hfbasis/  : spwf-related quantities in the Hartree-Fock basis 
     | | > /wavefunctions/compbasis/: spwf-related quantities in the computational basis 
     | | > /wavefunctions/canbasis/ : spwf-related quantities in the canonical basis
     | 
     | > /fields/        : functions on the mesh
     | | > /fields/densities/   : mean-field densities on the mesh
     | | > /fields/potentials/  : mean-field potentials on the mesh
     | 
     | > /multipoles/     : quantities related to multipole moments
     
## '/' group

 The root group contains all quantities required to specify the details of a MOCCa calculation, i.e. settings that should be known at the start of a calculation. These are all small quantities in terms of memory use and are all stored as HDF5 attributes.

 A complete list 

    - version             : string, describing the format version. Note: this is NOT the MOCCa version.
    - description         : string, string giving some comments. MOCCa writes an empty one for now. 
    - neutrons, protons   : double, average particle numbers for protons and neutrons. 
    - nx, ny, nz          : integers, numbers of mesh points. 
                            NOTE: these are the TOTAL number of mesh points in the whole of the simulation volume
    - dx, dy, dz          : double, mesh spacing in three directions; units of fm.
    - nwn, nwp, nwt       : integers,  number of neutron, proton and total wavefunctions
    - hfblocks            : integers, number of single-particle wavefunctions in each symmetry block.
    - param_name          : string, name of the Skyrme parameterization employed 
    - func_name           : string, name of the Skyrme functional employed 
    - diagsphamil         : integer, (1) the spwfs in MOCCa memory represent the Hartree-Fock basis 
                                     (0) the spwfs in MOCCa memory are not eigenstates of the single-particle hamiltonian
    - SYM_CODE            : character, encodes the precise symmetry options of the MOCCa run. 
    - mesh_type           : character, half-integer or integer
                                           |              |
                                           | > function values are stored in the middle of the discretisation cells
                                                          | > function values are stored on the 'left' of discretisation cells
                            for MOCCa calculations, this corresponds to calculations with anti-periodic and periodic boundary conditions.
    - BCtype              : character, anti-periodic or periodic.
    - calctype            : character, PASTA or NUCLEI.
    - pairingtype         : integer
                             0 -> Hartree-Fock
                             1 -> Hartree-Fock + BCS
                             2 -> Hartree-Fock-Bogoliubov
    - FermiEnergy         : double precision, rank 1 array of dimension 2 - first neutrons, then protons
    - omega               : double, rank 1 array of dimension 3; contains the cranking frequencies (\omega_x, \omega_y, \omega_z) in units of MeV \hbar^-1.

    THE ITEMS BELOW ARE NOT YET IMPLEMENTED
    ---------------------------------------
    - blocknumber         : integer, number of quasiparticle excitations were constructed. 
    - blocktype           : integer, specifies the type of blocking employed in the calculation 
                            0 -> no blocking 
                            1 equal filling blocking, asking for specific configuration
                            2 ordinary blocking, asking for lowest energy configuration
                            3 equal filling blocking, asking for specific configurations 
                            4 equal filling blocking, asking for lowest energy configuration 
                            5 ordinary blocking, spherical averaging; highly experimental feature.
                            6 UNUSED at the moment, reserved for future use.
                            7 ordinary blocking, index selection through overlap with a tagging state
    - blocklowest         : rank 1 string array of dimension (blocktype) ; exists only for blocktype = 2,4.
                            Possible entries: 'n+', 'n-', 'p+', 'p-', 'n0', 'p0'
    - blockindices        : rank 1 integer array of dimension (blocktype); exists only if blocktype = 1,3,5,7. 


## '/fields' group

 The fields group contains different quantities that live on the mesh - typically outputs of the calculation. 
 
 The mean-field densities and potentials as used by MOCCa are named according to the scheme proposed in [Ryssens et al., PRC 104 (2021)][https://link.aps.org/doi/10.1103/PhysRevC.104.044308]; these 'raw' quantities are all stored as *flat* arrays, i.e. arrays with rank 1! You 

 For example
    
    D_I_I   -> a (nx*ny*nz*2)   double precision array
    D_Nm_Nm -> a (nx*ny*nz*2)   double precision array
    C_I_N   -> a (nx*ny*nz*3*2) double precision array

 The last factor two in these quantities is an isospin index: first are neutrons, then protons.
 
 For user-friendliness, we also create a number of links to these objects that reflect the standard nomenclature; at the time of writing, these are
 
    rho -> D_I_I
    j   -> D_I_N
    s   -> C_I_N

 these links are stored within the '/fields' group itself. 
 
 The structure of the whole group is thus
 
    - rho
    - j
    - s
    - /densities/
      - /densities/D_I_I
      - /densities/D_Nm_Nm
      - ... etc ...
    - '/potentials/'
      - /densities/F_I_I
      - /densities/F_Nm_Nm
      - ... etc ...

## '/wavefunctions' group

 This group contains the information on the many-body state; first and foremost the actual values of the single-particle states. This group contains further subgroups that contain information 
 about these states in (a) the computational basis that is actually in memory, (b) the Hartree-Fock basis and (c) the canonical basis; note that all or some of these might coincide. Relevant 
 information is stored stored in subgroups /calcbasis/, /hfbasis/ and /canbasis/. Finally, there is information on the blocking configuration - needed to be able to restart calculations efficiently.

    /wavefunctions/
        - states: a double precision array of rank three containing the values of the single-particle wavefunctions on the mesh. 
              Dimensions: ( nx*ny*nz, 4, nwt ) 
                             |        |   |  wavefunction index
                             |        | spinor components
                             | spatial indices 
              By default, MOCCa stores the spwfs in the *computational* basis!

    -/wavefunctions/compbasis/ 
        - sphamil       : double precision rank 2 array of dimensions (nwn+nwp,nwn+nwp) containing the single-particle hamiltonian in the computational basis 
        - gaps          : double precision rank 2 array of dimensions (nwn+nwp,nwn+nwp) containing the full matrix of pairing gaps in the computational basis
        - rho_pairing   : double precision rank 2 array of dimensions (nwn+nwp,nwn+nwp) containing the density matrix \rho in the computational basis
        - kappa_pairing : double precision rank 2 array of dimensions (nwn+nwp,nwn+nwp) containing the anomalous density matrix \kappa in the computational basis
        - Bogo          : double precision rank 2 array of dimensions (nwn+nwp,nwn+nwp) containing the full Bogoliubov transformation in the computational basis

    -/wavefunctions/hfbasis/
        - spenergies    : double precision 1D array of dimension (nwn+nwp) that contains the diagonal matrix elements of the single-particle Hamiltonian in the Hartree-Fock basis 
        - HFtransfo     : double precision rank 2 array of dimensions (nwn+nwp,nwn+nwp) that contains the unitary transformation from the spwfs stored in states 
                          and those in the Hartree-Fock basis; this transformation might be trivial.
    -/wavefunctions/canbasis/
        - rho_can       : double precision 1D array of dimension (nwn+nwp) that contains the diagonal matrix elements of \rho in the canonical basis
        - cantransfo    : double precision rank 2 array of dimensions (nwn+nwp,nwn+nwp) that contains the unitary transformation from the spwfs stored in states 
                         and those in the canonical basis; this transformation might be trivial. 
        - configmatrix  : double precision rank 1 array of dimension (nwn+nwp) containing the quasiparticle occupation factors (see documentation elsewhere)

## '/multipoles' group

 This groups primary aim is to allow for efficient restarts of calculations with constraints on multipole moments. 

    /multipoles/ 
        - nmultipole     : number of multipole moments included in this file, i.e. the dimension of all arrays below.
        - ell            : integer rank 1 array of dimension (nmultipole); quantum numbers l of the multipole moments Q_lm 
        - m              : integer rank 1 array of dimension (nmultipole); quantum numbers m of the multipole moments Q_lm
        - impart         : integer rank 1 array of dimension (nmultipole); real (0) or imaginary (1) part of Q_lm
        - constrainttype : integer rank 1 array of dimension (nmultipole); type of constraint imposed on this multipole moment 
                           0 -> unconstrained 
                           1 -> Augmented Lagrangian constraint 
                           2 -> Projection on feasible set
        - value          : double rank 2 array of dimension (nmultipole,2), expectation value of the multipole moment - first neutrons, then protons. 
        - constraint     : double rank 1 array of dimension (nmultipole) targetted value of the multipole moment
        - deviation      : double rank 2 array of dimension (nmultipole,2), difference between current value and constraint
        - multiplier     : double rank 1 array of dimension (nmultipole), Lagrange multiplier of the constraint
        - intensity      : double rank 1 array of dimension (nmultipole), intensity parameter of the readjustment speed

