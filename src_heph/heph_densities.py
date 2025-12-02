#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# 
# HOW TO USE
# =========== 
#-------------------------------------------------------------------------------
# RULE NR 1: Preferably not by itself, as the generation of appropriate 
#            densities is better handled by the HEPH_functional module. 
#
#-------------------------------------------------------------------------------
#
# The workhorse of this module is the function
#   GenDensityExpression
#       
# That, when given a left-operator, right-operator for the density
# and which derivatives are needed of this density returns
#  a) and automatically generated name for the density
#  b) a Declaration expression, correctly declaring the density in FORTRAN
#  c) an Initialization expression, allocating and initialising the density
#  d) an Expression expression, calculating all of the components needed of the
#     density. 
#  e) A derivation expression, calculating all of the asked for derivatives of
#     this density. 
#
# All of this is of course based currently on templates in that function, but 
# those can be trivially adapted to correspond to different FORTRAN programs. 
#
# Indices that need to be traced over can be declared in the Contractions array
# as tuples of contracted indices. 
# examples:
#     tau   => Nabla, Nabla, [(0,1)]                 => Sum_mu tau_mumu
#
#     T_mu  => Nabla, Combine(Nabla, Sigma), [(0,1)] => Sum_mu T_mumunu
#   
#     F_mu  => Nabla, Combine(Nabla, Sigma), [(1,2)] => Sum_mu T_numumu
#
#
# Some tips:
# *) Always call initdensities first before doing anything, otherwise the 
#    attributes of the operators are not properly initialized. 
# *) When combining operators, it is MANDATORY to keep the possible indices
#    of the spin operators on the RIGHT of every operator, so that indices can
#    be grouped as
#           [ mu, nu     ...  | kappa ...]
#             derivative      |  spin    ]
#
#    So
#        Combine (Nabla, Sigma) is allowed
#        Combine (Sigma, Nabla) is NOT and will result in incorrect output.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script works from the observation that the calculation of any density
# in MOCCa/EV8/CR8/etc... can be achieved as
#
#    X(i,it) = \sum_{ij} \sum_{k=1,4} weight_{ji}                     (1)
#                        [O^{L} Psi_{i}]_k [O^{R} Psi_{j}]_k
#
# where the index k indicates the component in storage in the following storage
# scheme
#
#     Psi = ( Psi_1 + i Psi_2 )
#           ( Psi_3 + i Psi_4 )
#
# The sum can be over single-particle wavefunctions in different bases, and the
# weight can in principle be any matrix. This allows also for the creating of 
# expressions for pairing densities and mixed densities for projection codes.
#
# Note that there are no "extra sums" in (1), in the sense that every component
# of O^{L} Psi_{i} appears exactly once in the expression. This allows for the 
# work of this script, and the absence of "current of currents" is very
# important in this respect.
#
# Five operators are currently defined in this file, with which we can construct
# all of the densities necessary.
#
#               Identity, Nabla, Sigma, Current and TR
#
# All take a list of indices (even though Identity and Current don't need them) 
# and a 4-vector of components. As output, they permutate the components 
# according to their operator action and the index, with possible signs. 
# The action of sigma_y for example is coded as
#
#  Sigma(2, [1,2,3,4]) = [ -4, 3, 2, 1]
#
# representing  
#
#  sigma_y  ( Psi_1 + i Psi_2 ) = (   Psi_4 - i Psi_3)
#           ( Psi_3 + i Psi_4 )   ( - Psi_2 + i Psi_1)
#
# Operators can also be combined using the Combine function, which returns
# a function that combines both operators. 
#
# In addition, all of these operators have the following attributes
#
#   Operator.derorder = Order of the derivative in the operator. Used to 
#                       a) determine whether to calculate based on spwfs, their
#                          first-order derivatives or second-order derivatives 
#                          in the FORTRAN code.
#                       b) which of the indices in the (possible heavily 
#                          composite) operator are derivative indices, and not
#                          spin indices.
#   Operator.dimension= Order of the operator, scalar (0), vector (1) or tensor
#                       of rank (Operator.dimension)
#
#   Operator.time
#   Operator.signature_x
#   Operator.signature_y
#   Operator.signature_z
#   Operator.parity
#
#                       Symmetries of the operator. They are numbers for scalar
#                       operators and matrices for higher-order operators.
#                       From these the behaviour of the densities under all of 
#                       the D^TD_2h operators can be determined.
#
#-------------------------------------------------------------------------------

from string          import Template
from math            import log, factorial
from .               import heph_symmetries 
from src_heph.heph_symmetries  import *
import numpy         as np
import itertools

# Import some string templates
import src_heph.fortran_templates.Densities_templates as ta

# Horizontal line for printing
line  = 80*"-"
#-------------------------------------------------------------------------------
# Tab-character for the fortran code.
# 4 spaces for W.R., but I can imagine other people have different standards.
tab           = '    '

# Names of the arrays storing the single-particle states in Tantalus, from 
# which the densities need to be calculated.
ArrayNames=['DenPsi', 'DendPsi', 'DenddPsi', 'DendddPsi']

#-------------------------------------------------------------------------------
# Array containing all the different densities needed. Note that this contains 
# all of the densities that will get summed by Tantalus.
Densities_needed   = []
deriv_needed       = [] # what external derivatives are needed for each element
                        # in Densities needed
intermediate_status= [] # whether or not the density is an intermediate object.
                        # YES: the density does not feature in the expression
                        #      of the energy density and does not have to be
                        #      explicitly stored in the density vector
                        # NO:  the density gets stored explicitly
#-------------------------------------------------------------------------------
# Indices over which sums are supposed to go in both the FORTRAN code and the 
# naming scheme.
sumindices      = ['m', 'k', 'q', 'o', 'l']
crossindices    = ['x', 'y', 'z']
#-------------------------------------------------------------------------------
# Strings indicating the (external) laplacian and derivative of a density.  
lapstring = 'Lap'
derstring = 'Der'
#-------------------------------------------------------------------------------

def initdensities():
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Initialize the predefined operators with correct atributes, before they
    # ever get called.
    global ArrayNames
        
    Identity.dimension   = 0
    Identity.derorder    = 0  
    Identity.parity      = np.array([1])
    Identity.time        = np.array([1])
    Identity.signature_x = np.array([1])   
    Identity.signature_y = np.array([1])
    Identity.signature_z = np.array([1])
    Identity.name        ='I'
    
    Current.derorder     = 0
    Current.dimension    = 0
    Current.parity       = np.array([ 1])
    Current.time         = np.array([-1])
    Current.signature_x  = np.array([ 1])   
    Current.signature_y  = np.array([ 1])
    Current.signature_z  = np.array([ 1])
    Current.name         = 'C'
    
    Nabla.dimension      = 1
    Nabla.derorder       = 1
    Nabla.parity         = np.array([-1,-1,-1])
    Nabla.time           = np.array([ 1, 1, 1])
    Nabla.signature_x    = np.array([ 1,-1,-1])   
    Nabla.signature_y    = np.array([-1, 1,-1])
    Nabla.signature_z    = np.array([-1,-1, 1])   
    Nabla.name           = 'N'
    
    Sigma.dimension      = 1 
    Sigma.derorder       = 0
    Sigma.parity         = np.array([ 1, 1, 1])
    Sigma.time           = np.array([-1,-1,-1])
    Sigma.signature_x    = np.array([ 1,-1,-1])   
    Sigma.signature_y    = np.array([-1, 1,-1])
    Sigma.signature_z    = np.array([-1,-1, 1])
    Sigma.name           = 'S'
    
    TR.derorder     = 0
    TR.dimension    = 0
    TR.parity       = np.array([ 1])
    TR.time         = np.array([ 1])
    TR.signature_x  = np.array([-1])   
    TR.signature_y  = np.array([-1])
    TR.signature_z  = np.array([-1])
    TR.name         = 'T'
    
def ProcessDensities(fname, src, target, so, fam_active, density_spwf_summation, dry_run=False, verbose=True):
    """
    Master routine calling the other routines based on a list of densities.
    Also prints output.

    Input:
        fname                 : name of the file to write the output to
        src                   : source directory where the template file is
        target                : target directory where the output file should be
        so                    : a set of symmetry options
        fam_active            : a boolean indicating if we are generating a 
                                finite-amplitude linear response code (True)
                                or a static mean-field code (False).
        density_spwf_summation: if True, calculate derivatives of densities 
                                through summation over spwfs
                                If False, calculate them directly from the 
                                wavefunctions.
        dry_run               : if True, do not write any files (default: False)
        verbose               : if True, print detailed output (default: True)

    Output:
        Declaration           : string that declares the density in FORTRAN
        Memory                : string that gets the contribution to memory 
                                requirements for this density
    """
    from src_heph.heph_substitute import substitute

    global line
    Declaration      = ''
    Spwf_Declaration = ''
    Sph_declaration  = ''
    Initialisation   = ''

    # Actually complicated stuff
    Expression                       = '' # calculation statements for mean-field densities
    Expression_offdiag_symmetric     = '' # calculation statements for symmetric parts of more general densities
    Expression_offdiag_antisymmetric = '' # calculation statements for antisymmetric parts of more general densities
    Expression_offdiag_pp            = '' # calculation statements for perturbed pairing densities

    Derivation                       = ''
    Derivation_pair                  = ''
    Derivation_offdiag_symmetric     = ''
    Derivation_offdiag_antisymmetric = ''

    Isospincoupl                     = ''
    Isospincoupl_symmetric           = ''
    Isospincoupl_antisymmetric       = ''

    Expression_sph_sym     = ''
    Expression_sph_antisym = ''
    BCSExpression          = ''
    HFBExpression          = ''

    Zeroing          = ''
    Cleaning         = ''
    MPI_REDUCE       = ''
    Add              = ''
    Multiply         = ''
    Memory           = ''
    Write            = ''

    if(verbose):
        print (line)
        print (' Densities necessary for the functional                                    ')
        print (line)
        print ('      Name       Calc?     TR    Intermediate?  DIM with / out  Derivative combs.      ')
        print (line)

    for i in range(len(Densities_needed)):
        den = Densities_needed[i]
        
        owith    = OrderOfDen(den)
        owithout = OrderOfDen(den, contract=False)
        T        = TimeDen(den)
        D        = deriv_needed[i]
        I        = intermediate_status[i]
        if(verbose):
            print ('%15s %4s   %6d     %5s        %6d %6d     ' % (den,'y', T, I, owith, owithout), D )

    if(verbose):
        print (line)
        print (' Symmetries of the densities')
        print ('           DEN   LARG  RARG   T     P    RX    RY    RZ    SX    SY    SZ')
        print (line)

    for i in range(len(Densities_needed)):
      den = Densities_needed[i]
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # "Diagonal" summation of particle-hole densities in the canonical basis
      #         + "offdiagonal" summation of particle-particle densities
      (e,dec,spwf_dec,ini,der,isoi,mpii,zeroi,memi,cleani,addi,multi,writei)  = \
      GenDensityExpression(Densities_needed[i],deriv_needed[i],intermediate_status[i], \
                          'wave'     ,'wave',                                          \
                          'der_index', 'der_index',so, fam_active,                     \
                           density_spwf_summation, silent=not verbose)

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # "Off-diagonal" summation of densities in the HF-basis
      # (a) symmetric part
      if((not so.timelike) or TimeDen(Densities_needed[i]) == 1):
        # Symmetric parts of the perturbation of time-even mean-field densities are time-even
        # Symmetric parts of the perturbation of time-odd  mean-field densities are time-odd
        off_diag_tuple  = \
        GenDensityExpression(Densities_needed[i],deriv_needed[i],intermediate_status[i], \
                                'wave_j'     , 'wave_i',                                   \
                                'der_index_j', 'der_index_i',so,fam_active,density_spwf_summation,
                                complex_component=+1,weight='weight_sym', silent=True)
        e_off_sym   = off_diag_tuple[0]
        der_sym     = off_diag_tuple[4]
        iso_off_sym = off_diag_tuple[5]
      else:
        e_off_sym   = ''
        der_sym     = ''
        iso_off_sym = ''

      # (b) antisymmetric part
      if((not so.timelike) or TimeDen(Densities_needed[i])==-1):
        # Antisymmetric parts of the perturbation of time-even mean-field densities are time-odd
        # Antisymmetric parts of the perturbation of time-odd  mean-field densities are time-even
        off_diag_tuple  = \
        GenDensityExpression(Densities_needed[i],deriv_needed[i],intermediate_status[i], \
                            'wave_j'     , 'wave_i',                                   \
                            'der_index_j', 'der_index_i',so,fam_active,density_spwf_summation,
                             complex_component=-1,weight='weight_asym', silent=True)
        e_off_asym   = off_diag_tuple[0]
        der_asym     = off_diag_tuple[4]
        iso_off_asym = off_diag_tuple[5]
      else:
        e_off_asym   = ''
        der_asym     = ''
        iso_off_asym = ''

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Expressions required to evaluate the matrix elements of the single-particle
      # hamiltonian and pairing matrix Delta based on density-like expressions
      e_sph_sym     = ''
      e_sph_antisym = ''
      (t,t,left,right,t,t) = ParseOperators(Densities_needed[i], so.timelike)
      # TODO:  - [ ] write documentation!
      #        - [X] refactor this thing with loops
      #        - [ ] add if statements for time-reversal necessities?

      # Dirty hack: I do not want to deal (YET) with the symmetrization of pairing densities
      if('P' in den):
        symmetrization_needed = False 
      elif(left == right):
        symmetrization_needed = False 
      else:
        symmetrization_needed = True

      if(symmetrization_needed):
        # Explicit symmetrisation is required
        for symsign in [-1,+1]:
            if((not so.timelike) or TimeDen(Densities_needed[i]) == 1  or  'P' in den):
                # Symmetric parts of the perturbation of time-even mean-field densities are time-even
                # Symmetric parts of the perturbation of time-odd  mean-field densities are time-odd
                sph_tuple = \
                GenDensityExpression(Densities_needed[i],deriv_needed[i],intermediate_status[i], \
                                    'wave_j', 'wave_i',                                     \
                                    'wave_j', 'wave_i',so,fam_active,density_spwf_summation,
                                    complex_component=+1, weight='potential', symmetrize=symsign, silent=True)
                e_sph_sym += sph_tuple[0]
            if(not(so.timelike) or TimeDen(Densities_needed[i])==-1 or 'P' in den):
                # Antisymmetric parts of the perturbation of time-even mean-field densities are time-odd
                # Antisymmetric parts of the perturbation of time-odd  mean-field densities are time-even
                sph_tuple = \
                GenDensityExpression(Densities_needed[i],deriv_needed[i],intermediate_status[i], \
                                    'wave_j', 'wave_i',                                     \
                                    'wave_j', 'wave_i',so,fam_active,density_spwf_summation,
                                    complex_component=-1, weight='potential', symmetrize=symsign, silent=True)
                e_sph_antisym += sph_tuple[0]       
      else:
        # No explicit symmetrisation needed
        if((not so.timelike) or TimeDen(Densities_needed[i]) == 1):
            sph_tuple = \
            GenDensityExpression(Densities_needed[i],deriv_needed[i],intermediate_status[i], \
                                'wave_j', 'wave_i',                                     \
                                'wave_j', 'wave_i',so,fam_active,density_spwf_summation,
                                complex_component=+1, weight='potential', symmetrize=0, silent=True)
            e_sph_sym += sph_tuple[0]
        if((not so.timelike) or TimeDen(Densities_needed[i])==-1):
            sph_tuple = \
            GenDensityExpression(Densities_needed[i],deriv_needed[i],intermediate_status[i], \
                                'wave_j', 'wave_i',                                     \
                                'wave_j', 'wave_i',so,fam_active,density_spwf_summation,
                                complex_component=-1, weight='potential', symmetrize=0, silent=True)
            e_sph_antisym += sph_tuple[0] # we only need the calculation of this density
      if(verbose):
          print (' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')

      if( not intermediate_status[i]):
       Declaration     = Declaration    + '\n' + dec

      Spwf_Declaration = Spwf_Declaration + '\n' + spwf_dec
      if('P' in den):
        #-------------------------------------------------------------------------------------
        # Pairing densities 
        #-------------------------------------------------------------------------------------
        # The BCS expression is diagonal in 'wave'
        BCSExpression = BCSExpression + '\n' + e
        # But we also need the HFB expression 
        # So we recall the routine with different 'wave' indices
        # This summation is blockwise, hence the 'si+' in the indices 

        silent_call = False 
        if(not verbose):
            silent_call = True
        (e,dec,spwf_dec,ini,der,isoi,mpii,zeroi,memi,cleani,addi,multi, writei)  = \
                     GenDensityExpression(Densities_needed[i], deriv_needed[i],    \
                                          intermediate_status[i],                  \
                                         'si+wave2' , 'si+wave'  ,                 \
                                         'der_index', 'der_index', so, fam_active, \
                                          density_spwf_summation, silent=silent_call)
        HFBExpression = HFBExpression + '\n' + e

        Expression_offdiag_pp  = Expression_offdiag_pp + e_sph_sym
        # TODO: there seems to be no treatment for derivatives of pairing densities yes!
      else:
        #-------------------------------------------------------------------------------------
        # Particle-hole densities
        #-------------------------------------------------------------------------------------
        Expression                       = Expression                       + '\n' + e
        Expression_offdiag_symmetric     = Expression_offdiag_symmetric     + '\n' + e_off_sym
        Expression_offdiag_antisymmetric = Expression_offdiag_antisymmetric + '\n' + e_off_asym

        Derivation                       = Derivation                       + der
        Derivation_offdiag_symmetric     = Derivation_offdiag_symmetric     + der_sym
        Derivation_offdiag_antisymmetric = Derivation_offdiag_antisymmetric + der_asym
        Expression_sph_sym               = Expression_sph_sym               + '\n' + e_sph_sym
        Expression_sph_antisym           = Expression_sph_antisym           + '\n' + e_sph_antisym

      if( not intermediate_status[i]):
        Initialisation = Initialisation + '\n' + ini

        Isospincoupl                 = Isospincoupl                    + isoi
        Isospincoupl_symmetric       = Isospincoupl_symmetric          + iso_off_sym
        Isospincoupl_antisymmetric   = Isospincoupl_antisymmetric      + iso_off_asym

        MPI_REDUCE     = MPI_REDUCE     + '\n' + mpii
        Zeroing        = Zeroing        +        zeroi
        Cleaning       = Cleaning       + '\n' + cleani
        Add            = Add            + '\n' + addi
        Multiply       = Multiply       + '\n' + multi
        Memory         = Memory         + '\n' + memi
        Write          = Write          + '\n' + writei
    
    if(verbose):
        print (line)

    # Substitute into the densities.f90 file.
    dic={}
    dic['SPWF_DECLARATION'] = Spwf_Declaration
    dic['INITIALIZATION'  ] = Initialisation

    if(fam_active):
        dic['FAM'] = 1
    else:   
        dic['FAM'] = 0

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # These are the strings detailing the calculation of densities in all possible contexts
    #
    # a) particle-hole densities
    dic['EXPRESSION'        ]               = Expression                        # single summation for mean-field calculations
    dic['EXPRESSION_OFFDIAG_SYMMETRIC']     = Expression_offdiag_symmetric      # double summation for symmetric part of more general densities
    dic['EXPRESSION_OFFDIAG_ANTISYMMETRIC'] = Expression_offdiag_antisymmetric  # double summation for symmetric part of more general densities
    dic['EXPRESSION_DELTA_PP']              = Expression_offdiag_pp             # double summation for perturbed pairing densities

    dic['ISOSPINCOUPL'              ] = Isospincoupl
    dic['ISOSPINCOUPL_SYMMETRIC'    ] = Isospincoupl_symmetric
    dic['ISOSPINCOUPL_ANTISYMMETRIC'] = Isospincoupl_antisymmetric

    dic['EXPRESSION_SPH_SYM']     = Expression_sph_sym         # expression for the density-like calculation of the matrix elements of sph
    dic['EXPRESSION_SPH_ANTISYM'] = Expression_sph_antisym     # expression for the density-like calculation of the matrix elements of sph
    # b) particle-particle densities in the BCS case
    dic['BCSEXPRESSION'   ] = BCSExpression            # single summation for mean-field calculations
    # c) particle-particle densities in the BCS case
    dic['HFBEXPRESSION'   ] = HFBExpression            # single summation for mean-field calculations
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    if(density_spwf_summation):
        dic['DERIVATION'              ] = ""
        dic['DERIVATION_SUM_SPWF_PH'  ] = Derivation
        dic['DERIVATION_SUM_SPWF_BCS' ] = Derivation_pair
        dic['DERIVATION_SUM_SPWF_HFB' ] = Derivation_pair
        dic['DERIVATION_OFFDIAG_SYMMETRIC']     = '' # Put in empty strings for the FAM expressions
        dic['DERIVATION_OFFDIAG_ANTISYMMETRIC'] = ''

    else:
        dic['DERIVATION'              ]         = Derivation + Derivation_pair
        dic['DERIVATION_OFFDIAG_SYMMETRIC']     = Derivation_offdiag_symmetric
        dic['DERIVATION_OFFDIAG_ANTISYMMETRIC'] = Derivation_offdiag_antisymmetric
        dic['DERIVATION_SUM_SPWF_PH'  ] = ""
        dic['DERIVATION_SUM_SPWF_BCS' ] = ""
        dic['DERIVATION_SUM_SPWF_HFB' ] = ""

    dic['ZEROING'         ] = Zeroing
    dic['CLEANING'        ] = Cleaning
    dic['MPIDEN']           = MPI_REDUCE
    dic['ADD'             ] = Add
    dic['MULTIPLY'        ] = Multiply
    dic['WRITEDENSITIES_HDF5'] = Write

    # Symmetry options
    if(so.timelike):
      dic['TR']  = ''
      dic['NTR'] = '!'
    else:
      dic['TR']  = '!'
      dic['NTR'] = ''

    axes = ['X', 'Y', 'Z']
    for k in range(3):
     if(so.ReduceAxes[k] == 1):
      dic['REDU%s'%axes[k]] = ' '
      dic['FULL%s'%axes[k]] = '!'
     else:
      dic['REDU%s'%axes[k]] = '!'
      dic['FULL%s'%axes[k]] = ' '

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Here we figure out the symmetries of the ordinary density rho = D_I_I
    # as that one is often necessary outside of the densities module
    #
    # TODO: refactor this to use the AxisReflection function
    axes = ['X', 'Y', 'Z']
    for k in range(3):
     if(so.ReduceAxes[k] == 1):
      dic['S%s_RHO'%axes[k]] ='+1'
     else:
      dic['S%s_RHO'%axes[k]] =' 0'

    for k in range(3):
     if(so.ReduceAxes[k] == 1):
       if(k == 0 or k == 1):
         dic['S%s_RHO_ANTISYM'%axes[k]] ='-1'
       else:
         dic['S%s_RHO_ANTISYM'%axes[k]] ='+1'
     else:
      dic['S%s_RHO_ANTISYM'%axes[k]] =' 0'

    # And the same for D_I_S
    for j in range(3):
      P   = AxisReflection(Identity, Sigma,[()],[(j)],so,False)   
      for k in range(3):  
        if(so.ReduceAxes[k] == 1):
         dic['S%s_S%s'%(axes[k],axes[j])] ='%s'%P[k]
        else:
         dic['S%s_S%s'%(axes[k],axes[j])] =' 0'
         
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Ugly manual checking if parity is part of the generator set and 
    #  signalling this to the FORTRAN code for use in the calculation of the 
    #  two-body center of mass correction.
    symdic  = populatesymmetries()
    dic['PBROKEN'] = ' '
    for sym in so.generators:
      if (sym == symdic['P']):
        dic['PBROKEN'] = '!'

    if(not dry_run):
        substitute(src+fname, target+fname, dic)           

    # Return the declaration of all densities for use in vectors.f90
    return Declaration, Memory

def ParseOperators(density, timelike, findindices=False):    
    """    
     Parse the operators that are used to construct a density from its name.
     
     density  :  name of the density, for example D_N_N
     timelike :  logical indicating if there is an antilinear, antihermitian 
                 symmetry that is conserved. If so, the pairing densities get
                 an EXTRA T on the left, on top of the one they already have.
     findices :  TO BE DOCUMENTED
    """
    left  = ''
    right = '' 

    der   = density.count(derstring)
    lap   = density.count(lapstring)

    split = density.replace(derstring + '_', '').replace(lapstring + '_', '')

    split = split.split('_')

    left = split[1]
    right= split[2]
    #---------------------------------------------------------------------------
    # Don't put the couplings in the definition of left- and right-operators
    for l in sumindices + crossindices:
        left  =  left.replace(l, '')
        right = right.replace(l, '')
       
    if('C' in density):      
            right = 'C' + right
    
    if('P' in density):
        # Add a timereversal operator on the left for pairing densities
        left  = left + 'T'
        # If however, an antilinear, antihermitian symmetry is conserved, we 
        # add an EXTRA time-reversal operator, ending in a net overall sign.
        if(timelike): 
          left  = left + 'T'
    #---------------------------------------------------------------------------
    # Find the coupling over the sumindices
    coupling  = []
    foundsums = []
    indices   = []
    for l in sumindices:
        c   = ()
        ind = 0
        for i in range(1,len(density)):      
            if(density[i] == l):
                c = c+ (ind,)
            if(density[i-1] in ['N', 'S']):
                ind = ind + 1 
        if(len(c) > 0) :
            coupling.append(c)
            indices.append(l)
        if(len(c) == 3):            
            foundsums.append(l)
    #---------------------------------------------------------------------------
    # Find the coupling over the crossindices, but only if not contracted with
    # another index
    cross  = []
    for l in crossindices:
        c   = ()
        ind = 0
        for i in range(len(density)):      
            if(density[i] == l):
                if(i == len(density) - 1):
                    c = c+ (ind,)
                elif (density[i+1] not in foundsums):
                    c = c+ (ind,)
            if(density[i] in crossindices+sumindices):
                ind = ind + 1 
        if(len(c) > 0) :
            cross.append(c)    
            indices.append(l)
    if(findindices):
        return(der, lap, left, right, coupling, cross, indices)
    else:
        return(der, lap, left, right, coupling, cross)

def ReconstructDensity(der, lap, left, right):
  """
    Reconstruct the string for a density starting from its decomposition.
    
    Example:
        der = 1, lap = 0, left = 'I', right = 'CNS'

        should be translated to
        
        Der_C_I_NS
        
    Attention: this function silently assumes that there is no index coupling 
               going on, i.e. there are no contractions or vector products
               involved.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Input:
      der  : number of external (single) derivatives
      lap  : number of external laplacians
      left : left operators
      right: right operators

      These are identical to the output of ParseOperators.

    Output:
      density : string representing the density
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  """

  # Derivatives  
  density = der * 'Der_' + lap * 'Lap_'
  # Dealing with a current?
  if('C' in right or 'C' in left):
    density = density + 'C'
  else:
    density = density + 'D'
  
  # Dealing with a pairing density?
  if('T' in right or 'T' in left):
    density = density + 'P'
    
  #Adding in the operators
  density = density + '_' + left.replace('T', '').replace('C', '')
  density = density + '_' + right.replace('T', '').replace('C', '')
  
  return density

def GenDensityExpression(denin,derivative_combinations,intermediate, 
                         leftwave     , rightwave     ,
                         left_der_wave, right_der_wave, so,
                         fam_active, density_spwf_summation,  
                         complex_component=0, weight = 'weight',
                         silent=False, symmetrize=0):
    """
      Generate all the necessary strings to plug into FORTRAN source code 
      template Densities.f90.

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Input : 
      * denin                  : Expression for the density
      * Derivative_combinations: declares the different external derivatives of 
                                 the density that need to be calculated 
      * intermediate           : is this density one that needs to be explicitly
                                 stored (False) or one that is just an intermediate
                                 object for the calculation of other densities (True).
                                 If the latter, a lot of the output string should be empty
      * leftwave, rightwave    : Strings indicating to the summation what the
                                 left and right spwf is.
      * left/right_der_wave    : Strings indicating to the summation what the
                                 left and right spwf is when a derivative is involved
      * so                     : a set of symmetry options
      * density_spwf_summation : if True, calculate derivatives of densities 
                                 through summation over spwfs
      * fam_active            : a boolean indicating if we are generating a
                                    finite-amplitude linear response code (True)
                                    or a static mean-field code (False).
      * complex_component      : determines how the contribution of a given 
                                 pair of (leftwf, rightwf) gets added to the
                                 total density 
                                 
                                 (-1) => this is the imaginary part of a complex number
                                 (0)  => the density is real
                                 (+1) => this is the real part of a complex number
      * weight                 : string for the weight associated with the sum
      * silent                 : If True  => don't print the symmetry output 
                                 If False => print symmetry output for the 
                                             reflection symmetries of the 
                                             generated density
      * symmetrize             : TODO document
      - - - - - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      Output :
      
        Expression     :  string that calculates the density
        Declaration    :  string that declares the density in the FORTRAN code
        spwf_dec       :  string that declares the temporary density for use
                          in the densit function
        Initialisation :  string that (if necessary) allocates the density
        Derivation     :  string that handles all the derivatives that need to 
                          be calculated of the density
        Isospincoupl   :  string that handles the calculation of isospin 
                          densities
        Zeroing        :  string that sets the density to zero
        Memory         :  string that gets the contribution to memory requirements
                          for this density
        Cleaning       :  string deallocating the density
        Write          :  string writing the density to the HDF5 file
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      
    """

    assert abs(complex_component) <= 1

    #---------------------------------------------------------------------------
    # Initialisation 
    Expression    = ''
    Declaration   = ''
    spwf_dec      = ''
    Initialisation= ''
    Derivation    = ''
    Zeroing       = ''
    Isospincoupl  = ''
    Cleaning      = ''
    Add           = ''
    Multiply      = ''
    MPI_reduce    = ''
    Memory        = ''
    Write         = ''
    #---------------------------------------------------------------------------
    # Parse the structure from the name
    density = denin
    (x, y, left, right, coupling, cross) = ParseOperators(density, so.timelike)
    #---------------------------------------------------------------------------
    # If symmetrize is non-zero, change left <-> right and the couplings
    # accordingly
    if(symmetrize == -1 ):
        switch_den = density.split('_')
        R = switch_den[2]
        L = switch_den[1]
        # We put all spin indices on the r.h.s.
        # This allows us to
        # 1. Be more efficient with "on the fly" derivatives
        #    (They can be grouped more)
        # 2. make it easier to figure out the symmetries of the intermediate
        #    functions for the "on the fly" derivatives
        if('S' in R):
          for l in sumindices:
            if('S' + l in R):
               R = R.replace('S'+l, '')
               L = L + 'S' + l
          if('S' in R):
            R = R.replace('S', '')
            L = L + 'S'
        switch_den = switch_den[0]+'_'+R+'_'+ L
        (x,y,left,right,coupling,cross) \
                                 = ParseOperators(switch_den,so.timelike)

    if(complex_component == -1):
      right = 'C' + right
      # Note how this takes into account automatically the -i in the antisymmetric part
      # of C-objects!
    
    # Construct the left/right operators
    operatordic = {}
    operatordic['I'] = Identity
    operatordic['N'] = Nabla
    operatordic['S'] = Sigma
    operatordic['C'] = Current
    operatordic['T'] = TR
    
    LeftOperator = Identity
    for i in range(len(left)):
        # this needs to be done in reverse order
        l = left[len(left) - i -1 ] 
        LeftOperator  = Combine(operatordic[l], LeftOperator)
        
    RightOperator = Identity
    for i in range(len(right)):
        r = right[(len(right)) -i -1] 
        RightOperator = Combine(operatordic[r], RightOperator)
              
    #---------------------------------------------------------------------------
    # Declaration and initialisation, also for the derivatives.
    dic= {}
    dic['NAME']    = density

    dic['LEFTWF']    = ArrayNames[ LeftOperator.derorder]
    dic['RIGHTWF']   = ArrayNames[RightOperator.derorder]
    if(LeftOperator.derorder == 0) :
        dic['LEFTWAVE']  = leftwave
    else:
        dic['LEFTWAVE']  = left_der_wave

    if(RightOperator.derorder == 0) :
        dic['RIGHTWAVE']  = rightwave
    else:
        dic['RIGHTWAVE']  = right_der_wave

    totalind= ''
    dim     = ''
    
    #---------------------------------------------------------------------------
    # Find the correct dimensions
    ndim = LeftOperator.dimension + RightOperator.dimension 
    # Subtract two dimensions for every summation, and only one for every 
    # vector product    
    ndim = ndim - 2*len(coupling) - len(cross)

    # But we have double-counted possible contractions of vector products,
    # couplings of length 3
    threedim=0
    threecoupl = []
    for c in coupling:
        if(len(c) == 3):
            threedim = threedim + 1
            threecoupl.append(c)
    ndim = ndim - threedim
    
    totalind=''
    dim     =''
    sizecount = 1
    for r in range(ndim):
        totalind   = totalind + ',:' 
        dim        = dim      + ',3'
        sizecount  = sizecount * 3
    dic['TOTALIND']= totalind
    dic['DIM']     = dim
    if(len(dim)>1):
        dic['DIM_spwf'] = '(' + dim[1:] + ')' # remove leading comma and add brackets
    else:
        dic['DIM_spwf'] = ''
    # SIZE of the density to pass onto the MPI_ALLREDUCE call
    # the factor two reflects isospin
    dic['TRANS_SIZE'] = '%d*mv'%(sizecount*2) 
    
    if('P' not in density):
      dic['ISOSIZE'] = 4 # Full complement of neutron, proton, isoscalar, isovector
    else:
      dic['ISOSIZE'] = 2 # Only neutron, proton components for pairing densities

    #---------------------------------------------------------------------------
    # Get the declaration of the density and its derivatives right
    if( not intermediate ):
        if(fam_active):
            Declaration    = ta.Dec_complex.substitute(dic) # Densities are complex in the FAM code
        else:
            Declaration    = ta.Dec_real.substitute(dic)         # ... but real in the mean-field code

        MPI_reduce     = ta.mpi.substitute(dic)
        Initialisation = ta.Ini.substitute(dic)
        Zeroing        = ta.Zero_template.substitute(dic)
        Cleaning       = ta.Clean_template.substitute(dic)
        Add            = ta.Add_template.substitute(dic)
        Multiply       = ta.Multiply_template.substitute(dic)
        Write          = ta.write_template.substitute(dic)
        # memory requirement for this density
        Memory         = ta.Memory.substitute(dic)
    spwf_dec    = ta.Dec_spwf.substitute(dic)

    for c in derivative_combinations:
        l = c[0]
        d = c[1]
        
        if( l == 0 and d == 0):
            continue
        
        dic['NAME']    = l*'Lap_' + d*'Der_' + density
        if(d>0):
            dic['TOTALIND']= totalind + ',:'
        else:
            dic['TOTALIND']= totalind
        # Instead of adding dimensions to the FORTRAN array, we add simply
        # number of derivatives, since they are a completely symmetric tensor
        derind = Number_symmetric(3,d)
        if(derind != 1):
            dic['DIM']     = ',' + str(int(derind))  + dim  
        else:
            dic['DIM']     = dim

        if(fam_active): 
            Declaration    = Declaration    + '\n' + ta.Dec_complex.substitute(dic)
        else:   
            Declaration    = Declaration    + '\n' + ta.Dec_real.substitute(dic)    
        if(density_spwf_summation):
          # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          # Dealing with the MPI ALLREDUCE call for the derivatives
          sizecount = 1
          for r in range(ndim):
            sizecount  = sizecount * 3
          sizecount = sizecount * Number_symmetric(3,d) 
          # SIZE of the density to pass onto the MPI_ALLREDUCE call
          # the factor two reflects isospin
          dic['TRANS_SIZE'] = '%d*mv'%(sizecount*2)
          MPI_reduce     = MPI_reduce     + '\n' + ta.mpi.substitute(dic)
          Memory         = Memory         + '\n' + ta.Memory.substitute(dic)
          # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

        Initialisation = Initialisation + '\n' + ta.Ini.substitute(dic)
        Zeroing        = Zeroing         + ta.Zero_template.substitute(dic)
        Cleaning       = Cleaning + '\n' + ta.Clean_template.substitute(dic)
        Add            = Add      + '\n' + ta.Add_template.substitute(dic)
        Multiply       = Multiply + '\n' + ta.Multiply_template.substitute(dic)
        dic['NAME']    = density

    if((not so.timelike) or (leftwave != rightwave and 'P' not in density) or TimeDen(density) == +1):
        #---------------------------------------------------------------------------
        # Generate the expression to calculate the density
        # Start from standard wave-functions, [ Psi_1, Psi_2, Psi_3, Psi_4 ]^T
        start = np.zeros((4, 1))
        start[0, 0] = 1
        start[1, 0] = 2
        start[2, 0] = 3
        start[3, 0] = 4

        #---------------------------------------------------------------------------
        # TODO: define a dedicated function for this next particular piece of 
        #       code. It is highly complex, yet identical to the one in the 
        #       src_heph/heph_fields module. It should be put in one spot
        #       and abstracted.
        #---------------------------------------------------------------------------

        # Construct an iterator with all possible combinations of uncontracted 
        # indices. Note that the ordering is [scalar indices, vector_indices]
        # Note that it is not important in which order they are, 
        # since we loop over all of them, only that they are consistently applied. 
        args = itertools.product(range(3), repeat=ndim)

        Expression = Expression + ta.Den_line.substitute(dic)
        if not intermediate:
            Expression = Expression + ta.Den_comment.substitute(dic)
        else:
            Expression = Expression + ta.Den_intermediate_comment.substitute(dic)

        for arg in args:
            # We have the uncontracted indices. Now construct the combinations of
            # indices, including contracted ones, that correspond to this. 
            uncontracted = []
            if len(coupling) + len(cross) == 0:
                # Nothing to do if no couplings needed
                uncontracted = [arg]
            else:
                #-------------------------------------------------------------------
                # These are all of the combinations needed for the summation indices
                cont = itertools.product(range(3), repeat=len(coupling))

                # All of the possibilities for the vector products
                crossind = []
                for i in range(len(cross)):
                    crossind = crossind + (Rot_ind(arg[len(coupling) + i]))

                # Combine all possibilities
                if len(crossind) != 0:
                    # all the combinations , including scalar and vector products
                    fullcont = []
                    for s in cont:
                        for r in crossind:
                            fullcont.append(s + r)
                elif (threedim != 0 and len(crossind) == 0):
                    fullcont = []
                    for s in cont:
                        # Modify the number of possibilities for couplings between vector
                        # products and scalar products. 
                        threeopt = itertools.product(range(2), repeat=threedim)
                        for r in threeopt:
                            fullcont.append(s + r)
                elif (threedim != 0 and len(crossind) != 0):
                    fullcont = []
                    for s in cont:
                        for r in crossind:
                            threeopt = itertools.product(range(2), repeat=threedim)
                            for t in threeopt:
                                fullcont.append(s + r + t)
                else:
                    # No vector indices, and no scalar-vector
                    fullcont = cont

                #-------------------------------------------------------------------
                # Note that now fullcont contains all of the terms needed for the
                # particular argument of the left-hand side. 
                #
                #  The ordering of the indices is:
                #    
                #  ( mu, nu, ...., xsi , mx, nx, ....,zx  ,  ex, ey, ....., ez )  
                #   < scalar indices >  < vector indices >  < contracted vectors)
                #
                #  corresponding to things of the form
                #
                #  coupling, (0,1)        cross (0,1)         coupling (0,1,2)
                #
                #  meaning 
                #
                #  the value of the   | the values of the  |  whether it is the 
                #  indices in the     | vector indices     |  first term or the 
                #  summation          |                    |  second in the vector
                #                     |                    |  product
                #-------------------------------------------------------------------
                for c in fullcont:
                    p = ()
                    ii = 0
                    for i in range(LeftOperator.dimension + RightOperator.dimension):
                        found = False
                        for combination in coupling:
                            if i in combination:
                                if len(combination) == 2:
                                    p = p + (c[coupling.index(combination)],)
                                    found = True
                                elif len(combination) == 3:
                                    found = True
                                    if i == combination[0]:
                                        p = p + (c[coupling.index(combination)],)
                                    elif i == combination[1]:
                                        rot = Rot_ind(c[coupling.index(combination)])
                                        o = threecoupl.index(combination)
                                        p = p + (rot[c[-1 - o]][0],)
                                    elif i == combination[2]:
                                        rot = Rot_ind(c[coupling.index(combination)])
                                        o = threecoupl.index(combination)
                                        p = p + (rot[c[-1 - o]][1],)
                        for combination in cross:
                            if i == combination[0]:
                                p = p + (c[cross.index(combination) + len(coupling)],)
                                found = True
                            if i == combination[1]:
                                p = p + (c[cross.index(combination) + len(coupling) + 1],)
                                found = True
                        if not found:
                            p = p + (arg[ii],)
                            ii = ii + 1
                    uncontracted.append(p)
            IND = ''
            for mu in arg:
                IND = IND + ',' + str(int(abs(mu) + 1))  # Python indexes 0:N-1

            dic['IND'] = IND
            if len(dic['IND']) > 1:
                dic['IND_nocomma'] = '(' + IND[1:] + ')'  # identical, but without leading comma and in brackets
            else:
                dic['IND_nocomma'] = ""

            Expression = Expression + ta.Den_1_spwf.substitute(dic)

            #-----------------------------------------------------------------------
            # Now loop over the uncontracted indices
            for true_arg in uncontracted:
                # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                # Check if the indices for contraction are not superfluous
                # 
                # The ugly tuple(np.abs( construction is simply because abs doesn't 
                # accept tuples as arguments, for whatever reasons.
                larg = tuple(np.abs(true_arg[:LeftOperator.dimension]))
                rarg = tuple(np.abs(true_arg[LeftOperator.dimension:]))

                #-------------------------------------------------------------------
                # print some output on the densities
                (px, py, pz) = AxisReflection(LeftOperator, RightOperator, larg, rarg, so, 'P' in denin)

                direction = ['x', 'y', 'z']
                pl = ''
                for i in larg:
                    pl += direction[i]
                pr = ''
                for i in rarg:
                    pr += direction[i]

                mu = larg
                if len(larg) == 0:
                    mu = (0)
                nu = rarg
                if len(rarg) == 0:
                    nu = (0)
                T = LeftOperator.time[mu] * RightOperator.time[nu]
                par = LeftOperator.parity[mu] * RightOperator.parity[nu]

                try:
                    if len(T) > 0:
                        T = T[0]
                        par = par[0]
                except:
                    pass

                if not silent:
                    print(r'%15s %4s %4s  $%+d$ & $%+d$ & $%+d$ & $%+d$ & $%+d$ & $%+d$ & $%+d$ & $%+d$ \\'
                        % (denin, pl, pr, T, par, par * int(px), par * int(py), par * int(pz), int(px), int(py), int(pz)))
                #-------------------------------------------------------------------

                # Get the index of the reduced storage scheme for all of the 
                # derivative indices
                larg_stor = Storage_Mapping(larg[:LeftOperator.derorder])
                rarg_stor = Storage_Mapping(rarg[:RightOperator.derorder])

                leftind = LeftOperator(larg, start)
                rightind = RightOperator(rarg, start)

                lcolumns = leftind.shape[1]
                rcolumns = rightind.shape[1]

                LIND = ''
                RIND = ''

                # only nablas are symmetric
                if len(larg[:LeftOperator.derorder]) > 0:
                    LIND = LIND + ',' + str(int(larg_stor + 1))
                if len(rarg[:RightOperator.derorder]) > 0:
                    RIND = RIND + ',' + str(int(rarg_stor + 1))

                dic['LIND'] = LIND
                dic['RIND'] = RIND

                for i in range(4):
                    SIGN = np.sign(leftind[i]) * np.sign(rightind[i])
                    # We add a minus sign for the pairing densities
                    if 'P' in denin:
                        SIGN = -SIGN

                    for l in true_arg:
                        if l == 0:
                            SIGN = SIGN
                        else:
                            SIGN = SIGN * np.sign(l)
                    if SIGN > 0:
                        dic['SIGN'] = '+'
                    else:
                        dic['SIGN'] = '-'
                    dic['RCOMP'] = str(int(abs(rightind[i, 0])))
                    dic['LCOMP'] = str(int(abs(leftind[i, 0])))
                    Expression = Expression + \
                        '& \n               &' + \
                        ta.Den_diag.substitute(dic)
            Expression = Expression + '\n'

            # Final summation
            if not intermediate:
                # Only sum for storage if the object is not intermediate
                if weight != 'potential':
                    dic['WEIGHT'] = weight
                    if complex_component == +1:
                        Expression = Expression + ta.Den_sum_realpart.substitute(dic) + '\n\n'
                    elif complex_component == -1:
                        Expression = Expression + ta.Den_sum_imagpart.substitute(dic) + '\n\n'
                    elif complex_component == +0:
                        Expression = Expression + ta.Den_sum_real.substitute(dic) + '\n\n'
                else:
                    if symmetrize == -1 and ('C' in denin):  # dirty hack!
                        mult = '(-0.5d0) *'
                    elif symmetrize == +1 or symmetrize == -1:
                        mult = '(+0.5d0) *'
                    else:
                        mult = ''
                    dic['WEIGHT'] = mult + 'F%' + density.replace('D', 'F').replace('C', 'G') + '(i' + IND + ',it)'

                    if('P' not in density):
                        # Calculation for the single-particle hamiltonian
                        if complex_component == +1:
                            Expression = Expression + ta.Sph_sum_realpart.substitute(dic) + '\n\n'
                        elif complex_component == -1:
                            Expression = Expression + ta.Sph_sum_imagpart.substitute(dic) + '\n\n'
                        elif complex_component == +0:
                            Expression = Expression + ta.Sph_sum_real.substitute(dic) + '\n\n'
                    else:
                        # Calculation for the pairing matrix
                        if complex_component == +1:
                            Expression = Expression + ta.Delta_sum_realpart.substitute(dic) + '\n\n'
                        elif complex_component == -1:
                            Expression = Expression + ta.Delta_sum_imagpart.substitute(dic) + '\n\n'
                        elif complex_component == +0:
                            Expression = Expression + ta.Delta_sum_real.substitute(dic) + '\n\n'

            # And add a line for the isospin coupling
            if 'P' not in density:
                Isospincoupl = Isospincoupl + ta.Den_iso_comment.substitute(dic)
                Isospincoupl = Isospincoupl + ta.iso_normal.substitute(dic)

            #-----------------------------------------------------------------------
            # Generate expressions for the calculation of derivatives of densities
            if len(derivative_combinations) > 1:
                Derivation = Derivation + ta.Den_line.substitute(dic)
                Derivation = Derivation + ta.Den_comment_deriv.substitute(dic)

            for c in derivative_combinations:
                if c == (0, 0):
                    continue
                Derivation = Derivation + ta.Den_comment_deriv_b % (c[0], c[1])

                deriv_args = list(itertools.product(range(3), repeat=c[1]))

                # Arguments of the derivative operators
                # name of the derivative
                if c[0] == 0:
                    dic['NAME'] = (c[1] - 1) * 'Der_' + density
                else:
                    dic['NAME'] = (c[0] - 1) * 'Lap_' + (c[1]) * 'Der_' + density

                if density_spwf_summation:
                    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
                    # Add the calculation of derivatives through summation of
                    # the extra densities.
                    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    if c[1] != 0:  # There is a gradient here
                        # Get only the independent derivative operations
                        deriv_args = sorted(list(set([tuple(sorted(da)) for da in deriv_args])))
                        for darg in deriv_args:
                            directions = ['X', 'Y', 'Z']
                            dic['DIR'] = directions[darg[0]]

                            if c[0] > 0:
                                print("HEPHAESTOS cannot yet combine laplacians and gradients with DENSUM=1.")
                            else:
                                if c[1] == 1:
                                    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                    # There is no laplacian and but a single gradient.
                                    lleft = 'N' + left.replace('I', '')
                                    rright = 'N' + right.replace('I', '')

                                    # nabla D^L,R = D^nabla L, R + D^L, nabla R 
                                    #   Den       = LEFTDEN      + RIGHTDEN
                                    dic['LEFTDEN'] = ReconstructDensity(c[1] - 1, c[0], lleft, right)
                                    dic['RIGHTDEN'] = ReconstructDensity(c[1] - 1, c[0], left, rright)

                                    # DERIND       , DERLIND, DERRIND
                                    # Note that the first uses reduced storage mapping due to 
                                    # symmetries of multiple derivative operators.
                                    # LEFTIND and RIGHTIND do not use this yet!
                                    dic['DERIND'] = ',' + str(int(Storage_Mapping(darg) + 1)) + IND

                                    dargstring = ''
                                    for dargind in darg:
                                        dargstring = dargstring + str(dargind + 1) + ','
                                    dargstring = dargstring[:-1]

                                    # The nabla operator is added to the start of left, i.e.
                                    # it can be added to the left of the indices string.
                                    dic['DERLIND'] = dargstring + IND
                                    # The nabla operator is added to the start of right; we 
                                    # have to find out which indices are "left" and "right"
                                    # and insert the dargstring in between
                                    ldim = LeftOperator.dimension
                                    dic['DERRIND'] = IND[:ldim] + dargstring + IND[ldim:]

                                    Derivation = Derivation + ta.Der_sum.substitute(dic)
                                    # Add a line for the isospin coupling while we are here
                                    if 'P' not in density:
                                        Isospincoupl = Isospincoupl + ta.iso_der.substitute(dic)
                                elif c[1] == 2:
                                    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                                    # There is no laplacian but a double gradient.
                                    lleft = 'NN' + left.replace('I', '')
                                    rright = 'NN' + right.replace('I', '')

                                    # NN D^L,R = D^NN L, R + D^L, NN R + 2 D^N L, N R 
                                    #   Den       = LEFTDEN      + RIGHTDEN
                                    dic['LEFTDEN'] = ReconstructDensity(0, 0, lleft, right)
                                    dic['RIGHTDEN'] = ReconstructDensity(0, 0, left, rright)

                                    lleft = 'N' + left.replace('I', '')
                                    rright = 'N' + right.replace('I', '')
                                    dic['CENTRALDEN'] = ReconstructDensity(0, 0, lleft, rright)

                                    # DERIND       , DERLIND, DERRIND
                                    # Note that the first uses reduced storage mapping due to 
                                    # symmetries of multiple derivative operators.
                                    # LEFTIND and RIGHTIND do not use this yet!
                                    dic['DERIND'] = ',' + str(int(Storage_Mapping(darg) + 1)) + IND

                                    # On the left, the new Nabla's precede all other arguments
                                    dic['DERLIND'] = str(darg[0] + 1) + ',' + str(darg[1] + 1) + IND
                                    # On the right, the new Nabla's precede only the indices of the right operator
                                    ldim = LeftOperator.dimension
                                    if ldim > 0:
                                        dic['DERRIND'] = IND[0:ldim] + ',' \
                                            + str(darg[0] + 1) + ',' \
                                            + str(darg[1] + 1) + \
                                            IND[ldim:]
                                    else:
                                        dic['DERRIND'] = str(darg[0] + 1) + ',' \
                                            + str(darg[1] + 1) + \
                                            IND[ldim:]

                                    # In the center, one nabla is on the left
                                    dic['DERCIND'] = str(darg[0] + 1) \
                                        + IND[0:ldim] + ',' \
                                        + str(darg[1] + 1) + IND[ldim:]
                                    Derivation = Derivation + ta.Der_der_sum.substitute(dic)

                                    # Add a line for the isospin coupling while we are here
                                    if 'P' not in density:
                                        Isospincoupl = Isospincoupl + ta.iso_der.substitute(dic)
                                else:
                                    print('Hephaestos cannot yet handle more than 2 gradients with DENSYM!')
                                    sys.exit(1)

                    else:  # there is only a Laplacian here
                        lleft = 'NN' + left.replace('I', '')
                        rright = 'NN' + right.replace('I', '')

                        # Delta D^L,R = D^Delta L, R + D^L, Delta R + 2 D^Nabla L, Nabla R 
                        #   Den       = LEFTDEN      + RIGHTDEN
                        dic['LEFTDEN'] = ReconstructDensity(c[1] - 1, 0, lleft, right)
                        dic['RIGHTDEN'] = ReconstructDensity(c[1] - 1, 0, left, rright)

                        lleft = 'N' + left.replace('I', '')
                        rright = 'N' + right.replace('I', '')
                        dic['CENTRALDEN'] = ReconstructDensity(c[1] - 1, 0, lleft, rright)

                        dic['IND'] = IND
                        Derivation = Derivation + ta.Lap_sum_a.substitute(dic)
                        for k in range(3):
                            dic['DERLIND'] = str(k + 1) + ',' + str(k + 1) + IND
                            dic['DERRIND'] = str(k + 1) + ',' + str(k + 1) + IND
                            dic['DERCIND'] = str(k + 1) + ',' + str(k + 1) + IND
                            Derivation = Derivation + ta.Lap_sum_b.substitute(dic)
                        Derivation = Derivation[:-2]

                        # Add a line for the isospin coupling while we are here
                        if 'P' not in density:
                            Isospincoupl = Isospincoupl + ta.iso_lap.substitute(dic)
                else:
                    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    # Add the calculation of the derivatives of the original  
                    # density through calls to derivative routines.
                    #
                    # Note that larg and rarg need not be redefined here. They 
                    # take the value of the last combination of uncontracted 
                    # indices. This is sufficient, because necessarily all of the 
                    # combinations of larg and rarg need to exhibit the same 
                    # symmetries. 
                    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    if c[1] != 0:
                        # Get only the independent derivative operations
                        deriv_args = sorted(list(set([tuple(sorted(da)) for da in deriv_args])))
                        for darg in deriv_args:
                            directions = ['X', 'Y', 'Z']
                            dic['DIR'] = directions[darg[0]]

                            # Note that the symmetries put into a certain call to the 
                            # derivatives are determined by 
                            # a) left- and right-operator
                            # b) indices (arguments) of these
                            # c) but also all the arguments of previously applied
                            #    derivatives.
                            # 
                            # This makes this particular bit of code rather complicated.

                            # Decide if we need to use a gradient or a laplacian routine
                            if c[0] > 0:
                                # There is a Laplacian involved, and we first calculate
                                # all derivatives, and then only afterwards laplacians.

                                (px, py, pz) = AxisReflection(LeftOperator, RightOperator, larg, rarg, so, 'P' in denin, darg)
                                dic['PX'] = str(px)
                                dic['PY'] = str(py)
                                dic['PZ'] = str(pz)

                                if len(darg) > 0:
                                    dic['IND'] = ',' + str(int(Storage_Mapping(darg) + 1)) + IND
                                else:
                                    dic['IND'] = IND
                                Derivation = Derivation + ta.Lap.substitute(dic)

                                # Add a line for the isospin coupling while we are here
                                if 'P' not in density:
                                    Isospincoupl = Isospincoupl + ta.iso_lap.substitute(dic)

                            else:
                                # There is no laplacian, so we only calculate partial
                                # derivatives

                                (px, py, pz) = AxisReflection(LeftOperator, RightOperator, larg, rarg, so, 'P' in denin, darg[1:])
                                dic['PX'] = str(px)
                                dic['PY'] = str(py)
                                dic['PZ'] = str(pz)

                                syms = (px, py, pz)
                                dic['PS'] = syms[darg[0]]

                                dic['DERIND'] = ',' + str(int(Storage_Mapping(darg) + 1)) + IND
                                if len(darg) > 1:
                                    dic['IND'] = ',' + str(int(Storage_Mapping(darg[1:]) + 1)) + IND
                                else:
                                    dic['IND'] = IND
                                Derivation = Derivation + ta.Der_indep.substitute(dic)

                                # Add a line for the isospin coupling while we are here
                                if 'P' not in density:
                                    dic['DARG'] = str(darg)
                                    dic['DC'] = str(c)
                                    Isospincoupl = Isospincoupl + ta.iso_der.substitute(dic)
                    else:
                        # Pure Laplacian operators
                        (px, py, pz) = AxisReflection(LeftOperator, RightOperator, larg, rarg, so, 'P' in denin)
                        dic['PX'] = str(px)  # + 'd0'
                        dic['PY'] = str(py)  # + 'd0'
                        dic['PZ'] = str(pz)  # + 'd0'

                        dic['IND'] = IND
                        Derivation = Derivation + ta.Lap.substitute(dic)

                        # Add a line for the isospin coupling while we are here
                        if 'P' not in density:
                            Isospincoupl = Isospincoupl + ta.iso_lap.substitute(dic)

            if len(derivative_combinations) > 1:
                Derivation = Derivation + '\n'
                Isospincoupl = Isospincoupl + '\n'

            dic['NAME'] = density

        Expression = Expression + ta.Den_line.substitute(dic)
        
    return (Expression, Declaration, spwf_dec, Initialisation, Derivation, Isospincoupl,\
                                   MPI_reduce, Zeroing, Memory, Cleaning, Add, Multiply,\
                                   Write)
    
def GenVecProd(density, coupling):
    #---------------------------------------------------------------------------
    # Generate the appropriate expressions for the calculation of a vector
    # product. 
    #---------------------------------------------------------------------------

    import src_heph.fortran_templates.GenVecProd as tb
    
    index_encountered = 0
    name              = ''
    ini               = ''
    
    dic = {}
    
    #---------------------------------------------------------------------------
    # Constructing the correct name: inserting xi, xj etc at the correct spot.
    for i in range(len(density)):
        l = density[i]
        if( density[i:i+3] == 'Der'):
            for c in coupling: 
                if index_encountered in c:
                    name = name + l + 'x' #'%d'%(coupling.index(c)+1)
                else :
                    name = name + l
            index_encountered=index_encountered+1
        elif( (l.isupper() and l != 'I' and l != 'D' and l!= 'C')):
            for c in coupling:
                if index_encountered in c:
                    name = name + l + 'x' #'%d'%(coupling.index(c)+1)
            index_encountered=index_encountered+1
        else:   
            name = name + l
        
    dic['NAME'] = name
    
    #---------------------------------------------------------------------------
    # Constructing the correct set of indices for allocation and declaration.
    order      =  OrderOfDen(density)

    allocind       = ''
    allocindhdf5   = ''
    for i in range(order - len(coupling)):
        allocind     = allocind     + ',3'
        allocindhdf5 = allocindhdf5 + '*3'
    
    dic['ALLOCIND']     = allocind
    dic['ALLOCINDHDF5'] = allocindhdf5
    dic['DECLIND']  = allocind.replace('3', ':')
    
    decl = tb.decl_template.substitute(dic)
    ini  = tb.Ini.substitute(dic)
    
    #---------------------------------------------------------------------------
    # Constructing the calculation segments.
    # Generate all of the possible values for the remaining indices.
    # Note that in this generation, we take
    #  [ m n, ..., z, x1, x2 ]
    #    -----------  -------
    #    uncoupled    coupled  indices
    redorder = order - 2*len(coupling) 
    args     = itertools.product(range(3), repeat=order - len(coupling))
    
    calc = ''
    for arg in args:
        uncontracted = []
    
        if(len(coupling) == 0):
            uncontracted = [arg]
        else:
           comb = []
           for k in range(len(coupling)):
                (i,j) = Rot_ind(arg[redorder + k])
                # Add a minus sign to indicate which one of the combinations
                # is reversed.
                comb = comb + [[(i+1,j+1), (-j-1,i+1)]] 
  
           combinations = itertools.product(*comb)
           for c in combinations:
                p = ()
                ii = 0
                for i in range(order):
                    found = False                    
                    for k in range(len(coupling)):
                            if(i in coupling[k]): 
                                    # Look at this beast of an expression :)
                                    p = p + (c[k][coupling[k].index(i)],)  
                                    found = True
                    if(not found): 
                            p = p + (arg[ii]+1,)
                            ii = ii +1
                uncontracted.append(p)
                
                
        dic['ARG'] = ''
        dic['DEN'] = density
        for i in range(len(arg)):
            dic['ARG'] = dic['ARG'] + ',%s'%(arg[i]+1)
        calc = calc + tb.calc_a.substitute(dic)
        calc = calc + 2 * tab +'&'
        
        count = 0
        for true_arg in uncontracted:
            dic['IND']  = ''
            s = 1
            for i in range(len(true_arg)):
                dic['IND']  = dic['IND']+ ',%s'%(abs(true_arg[i]))
                s           = s * true_arg[i]
            if(s > 0 ): 
                dic['SIGN'] = '+'
            else:
                dic['SIGN'] = '-'
            calc = calc + tb.calc_b.substitute(dic)
            count = count + 1
            if(count%3 == 0):
                calc = calc + '& \n' + 2 * tab +'&'
        calc = calc  + '\n'       
            
    return (decl, ini, calc)

#===============================================================================
# Operator routines.
#===============================================================================

def Identity(mu,indices):
    
    return indices

def Nabla(mu,indices):
    # - - - - - - - - - - - - - - - - - - - -
    # Derivative in direction 
    # 
    # Takes every column  
    #   ( 1 )     ( 0  1  0 )
    #   ( 2 ) = > ( 0  2  0 )
    #   ( 3 ) = > ( 0  3  0 )
    #   ( 4 )     ( 0  4  0 )
    dim     = indices.shape
    rows    = dim[0]
    columns = dim[1]
    
    out = np.zeros((rows, columns))
    
    for i in range(columns):
        out[:,i] = indices[:,i]

    return out
    
def Sigma(mu,indices):
    # - - - - - - - - - - - - - - - - - - - -
    # Sigma in direction 
    
    dim     = indices.shape
    rows    = dim[0]
    columns = dim[1]
    
    out = np.zeros((rows, columns))
    try:
      if(mu[0] == 0):
        for i in range(columns):
            out[0,i] =   indices[2,i]
            out[1,i] =   indices[3,i]
            out[2,i] =   indices[0,i]
            out[3,i] =   indices[1,i] 
        return out
      elif(mu[0] == 1):
        for i in range(columns):
            out[0,i] =   indices[3,i]
            out[1,i] = - indices[2,i]
            out[2,i] = - indices[1,i]
            out[3,i] =   indices[0,i] 
        return out
      else:
        for i in range(columns):
            out[0,i] =   indices[0,i]
            out[1,i] =   indices[1,i]
            out[2,i] = - indices[2,i]
            out[3,i] = - indices[3,i] 
        return out

    except IndexError:
        print ('Index Error in sigma!')
        exit(1)

def Current(mu,indices):
    # Operates on indices to get a current C, instead of a density D
    out = np.zeros_like(indices)

    out[0,:] =   indices[1,:] 
    out[1,:] = - indices[0,:]
    out[2,:] =   indices[3,:]
    out[3,:] = - indices[2,:]
    
    return out

def TR(mu,indices):
    # Operates on indices to a time-reverse of the spwf
    # Implemented as
    # i sigma_y K 
    # where K is the complex conjugation
    out = np.zeros_like(indices)

    out[0,:] =   indices[2,:] 
    out[1,:] = - indices[3,:]
    out[2,:] = - indices[0,:]
    out[3,:] =   indices[1,:]
    
    return out

def Combine( L , R ):
    #---------------------------------------------------------------------------
    # Combines two symmetry-operators to define a new one. Note that the order
    # of the operators is important for non-Abelian systems.
    #---------------------------------------------------------------------------

    LR = lambda mu,x : L(mu[0:L.dimension],R(mu[L.dimension:],x))
    LR.dimension = L.dimension + R.dimension
    LR.derorder  = L.derorder  + R.derorder
    
    LR.parity = BroadCastSymmetries(L.parity, R.parity)
    LR.time = BroadCastSymmetries(L.time, R.time)
    LR.signature_x = BroadCastSymmetries(L.signature_x, R.signature_x)
    LR.signature_y = BroadCastSymmetries(L.signature_y, R.signature_y)
    LR.signature_z = BroadCastSymmetries(L.signature_z, R.signature_z)
    
    LR.name =(L.name + R.name).replace('I', '')
    
    return LR

def BroadCastSymmetries(L, R):
    #-----------------------------------------
    # Broadcasts symmetries of two operators.
    #-----------------------------------------
    S       = np.outer(L, R)
    order   = int(log(np.size(S),3))
    
    t= ()
    for r in range(order):
       t = t + (3,) 
    if(order != 0):
        S = S.reshape(t)
    return S
    
def AxisReflection(LO, RO, larg, rarg, so, pairing, nabla_arg = []):
    #---------------------------------------------------------------------------
    # Calculates the sign under axis reflection for the symmetries of an EV8/CR8
    # calculation.
    #
    #  LO        = left-operator
    #  RO        = right-operator
    #  larg      = indices of the left operator
    #  rarg      = indices of the right operator
    #  nabla_arg = indices of the nabla_operators acting possibly on the density
    #  so        = a set of symmetry options
    #  pairing   = whether or not we are dealing with a pairing density, i.e. 
    #              a density that gets "summed over kappa".
    #---------------------------------------------------------------------------

    # Make sure that the empty tuple get recognised as simply indicating a number
    if(len(larg) == 0):
        mu = (0)
    else :
        mu = larg
    if(len(rarg) == 0):
        nu = (0)
    else:
        nu = rarg

    # Startindices
    start  = np.zeros((4,1))
    start[0,0] = 1 
    start[1,0] = 2  
    start[2,0] = 3 
    start[3,0] = 4

    st_xp = np.zeros((4,1)) ; st_yp = np.zeros((4,1)) ;  st_zp = np.zeros((4,1))
    st_xm = np.zeros((4,1)) ; st_ym = np.zeros((4,1)) ;  st_zm = np.zeros((4,1))

    # Then figure out the transformation of symmetry
    # Currently, hardcoded CR8-like symmetries 
#    sxp = [+1,-1,-1,+1] ; sxm = [-1,+1,+1,-1] 
#    syp = [+1,-1,+1,-1] ; sym = [+1,-1,+1,-1]
#    szp = [+1,+1,-1,-1] ; szm = [-1,-1,+1,+1]

    if(so.ReduceAxes[0] == 1):
      sxp = multiply_quantum_numbers(so.syms[0].permutation, so.combs[0],1)
      sxm = multiply_quantum_numbers(so.syms[0].permutation, so.combs[0],3)
    else:
      sxp = [0,0,0,0]
      sxm = [0,0,0,0]

    if(so.ReduceAxes[1] == 1):
      syp = multiply_quantum_numbers(so.syms[1].permutation, so.combs[1],1)
      sym = multiply_quantum_numbers(so.syms[1].permutation, so.combs[1],3)
    else:
      syp = [0,0,0,0]
      sym = [0,0,0,0]

    if(so.ReduceAxes[2] == 1):
      szp = multiply_quantum_numbers(so.syms[2].permutation, so.combs[2],1)
      szm = multiply_quantum_numbers(so.syms[2].permutation, so.combs[2],3)
    else:
      szp = [0,0,0,0]
      szm = [0,0,0,0]

    for i in range(4):
        st_xp[i,0] = sxp[i] * start[i,0] ; st_xm[i,0] = sxm[i] * start[i,0] 
        st_yp[i,0] = syp[i] * start[i,0] ; st_ym[i,0] = sym[i] * start[i,0] 
        st_zp[i,0] = szp[i] * start[i,0] ; st_zm[i,0] = szm[i] * start[i,0] 

    # We take the right spwf to transform as a positive signature spwf 
    lind_x = LO(larg,st_xp)  ; lind_y = LO(larg,st_yp); lind_z = LO(larg,st_zp)
    if(pairing and not so.timelike):    
        rind_x = RO(rarg,st_xm); rind_y=RO(rarg,st_ym); rind_z = RO(rarg,st_zm)
    else:
        rind_x = RO(rarg,st_xp); rind_y=RO(rarg,st_yp); rind_z = RO(rarg,st_zp)

    lind   = LO(larg, start) ;  rind   = RO(rarg, start)

    # Get the ratios for the symmetry signs
    px = lind_x[:,0] * rind_x[:,0]/(lind[:,0] * rind[:,0])
    py = lind_y[:,0] * rind_y[:,0]/(lind[:,0] * rind[:,0])
    pz = lind_z[:,0] * rind_z[:,0]/(lind[:,0] * rind[:,0])

    #---------------------------------------------------------------------------
    # Taking into account the signs of the derivatives
    for i in range(LO.derorder):
        if mu[i] == 0:  
            px = -1 * px
        if mu[i] == 1:  
            py = -1 * py
        if mu[i] == 2:  
            pz = -1 * pz

    for i in range(RO.derorder):
        if nu[i] == 0:  
            px = -1 * px
        if nu[i] == 1:  
            py = -1 * py
        if nu[i] == 2:  
            pz = -1 * pz

    # Sanity check, are all ratios equal? 
   # for i in range(4):
   #     if(px[i] != px[0]):
   #         print ("Error determining symmetries, X!")
   #         exit()            
   #     if(py[i] != py[0]):
   #         print ("Error determining symmetries, Y!")
   #         exit()            
   #     if(pz[i] != pz[0]):
   #         print ("Error determining symmetries, Z!")
   #         exit()            

    px = px[0]
    py = py[0]
    pz = pz[0]

    #---------------------------------------------------------------------------
    # Symmetries of the left- and rightoperator
    #pxl = LO.signature_z[mu] * LO.parity[mu] * LO.signature_y[mu] * LO.time[mu]
    #pxr = RO.signature_z[nu] * RO.parity[nu] * RO.signature_y[nu] * RO.time[nu]
    #px  = pxl * pxr

    #pyl = LO.signature_y[mu] * LO.parity[mu] * LO.time[mu]
    #pyr = RO.signature_y[nu] * RO.parity[nu] * RO.time[nu]
    #py  = pyl * pyr
    
    #pzl = LO.parity[mu]  * LO.signature_z [mu]
    #pzr = RO.parity[nu]  * RO.signature_z[nu] 
    #pz  = pzl * pzr

    #---------------------------------------------------------------------------    
    #Taking into account extra external derivatives 
    for nablaind in nabla_arg:
        # Note that python indexes starting from 0, so (x,y,z) = (0,1,2)
        if(nablaind == 0):
            px = - px
        elif(nablaind == 1):
            py = - py
        elif(nablaind == 2):
            pz = - pz
    if(px > 0):
        px = '+1'
    elif(px < 0):
        px = '-1'
    else:
        px = '0'
        
    if(py > 0):
        py = '+1'
    elif(py < 0):
        py = '-1'
    else:
        py = '0'
    
    if(pz > 0):
        pz = '+1'
    elif(pz < 0):
        pz = '-1'
    else:
        pz = '0'

    return(px,py,pz)

#===============================================================================
# Auxiliary routines.  
#===============================================================================
def TimeDen(density):
    """
      Obtain the behavior under time-reversal of the density.

      Input:
        density: a string representing a local density
      Output:
        T      : +1 if the density is time-even
                 -1 is the density is time-odd
    """
    # I don't pass in the symmetry option, as an extra T cannot affect the 
    # end result
    (x, y, left, right, coupling, cross) = ParseOperators(density, False)

    # Construct the left/right operators
    operatordic = {}
    operatordic['I'] = Identity
    operatordic['N'] = Nabla
    operatordic['S'] = Sigma
    operatordic['C'] = Current
    operatordic['T'] = TR
    
    LeftOperator = Identity
    for i in range(len(left)):
        # this needs to be done in reverse order
        l = left[len(left) - i -1 ] 
        LeftOperator  = Combine(operatordic[l], LeftOperator)
        
    RightOperator = Identity
    for i in range(len(right)):
        r = right[(len(right)) -i -1] 
        RightOperator = Combine(operatordic[r], RightOperator)

    # The .time property for both left and right operators is an array in 
    # general, but all of the indices should have the same value, so we 
    # simply average
    T = np.average(LeftOperator.time) * np.average(RightOperator.time)

    return T

def OrderOfDen(density, contract=True):
    """
     Returns the order (= number of indices) of the density represented 
     by a string. Simply checks the number of capital letters N and S in 
     in the name. 
     Either: a) disregard contractions           contract = False
          or b) take into account contractions   contract = True
    
     Note that it is safe to use this on a field too.
    """

    # Add a dimension for every derivative and 
    # don't count the capital D,C,F,G, P 
    test  = density.replace('_', '')
    order = test.count(derstring) - 1 
    test  = test.replace(derstring, '')
    test  = test.replace(lapstring, '')
    
    for i in range(len(test)):
        letter = test[i]
        if(letter.isupper() and letter != 'I' and letter != 'P'):
            # Every capital letter that is not I adds an index
            order = order + 1
        if( (letter in sumindices + crossindices) and contract):
            # But if the a non-capital letter follows, the index is contracted.
            order = order - 1
    for l in crossindices:  
        if(density.count(l) == 2):
            order = order + 1
        
    return order
    
    
def Rot_ind(k):
    """
     C_k = A_i B_j + A_j B_i
     Input is k, output is the corresponding pair
         [ (i,j) ] if s == 0
         [ (j,i) ] if s == 1
     where (i,j) carries the plus sign 
    
     Note that the minus sign is only relative and always
     carried by the non-zero index
    """
    if(k == 0):
        i = 1
        j = 2
        return [(i,j), (-j,i)]
    elif(k == 1):
        i = 2
        j = 0
        return [(i,j), (j,-i)]
    elif(k == 2):
        i = 0
        j = 1
        return [(i,j), (-j,i)]
        
def Storage_Mapping(indices):
    """
     Map the indices (i,j,k,...) of a totally symmetric tensor unto indices
     that are used for efficient storage
    
     Note that negative numbers are treated as positive, in order to not
     upset the vector products.
    """

    # Sort the indices into lexicographical order
    s_indices = sorted(np.abs(indices)) 
    # np.abs, because normal abs doesn't accept tuples
    
    # Find the number of total elements that are possible
    k = len(indices)
    n = Number_symmetric(3,k)
    
    # Find the index in a lexicographical sorting scheme
    # We calculate how many indices are before the current one
    index = 0
    for i in range(k):
        for j in range(s_indices[i]):
            index = index + Number_symmetric(3-j-1,k-i-1)
    
    return(index)
    
def Number_symmetric(n,k):
    """
     Returns the number of independent elements in a totally symmetric tensor
     of order k over dimension n.
    """    
    return(factorial(n + k - 1)/(factorial(k) * factorial(n-1)))
    
def Multiplicity(indices):
    """
     Get the total number of independent combinations that can be obtained by
     permutation the indices.
    """
    a = set(list(itertools.permutations(indices)))
    m = len(a)
    return(m)
    
def Isospinindices(iso):
  """
    Translate a string indicating a specific component of a density (0 or 1 for 
    isoscalar and isovector densities, p/n for proton/neutron densities) into 
    an array index. In the FORTRAN code, densities are stored as
    
      D_X_(SPATIAL INDEX, CARTESIAN COMPONENTS, 1:4)
    
    where 1 : neutron density
          2 : proton density
          3 : isoscalar density
          4 : isovector density
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    Input:
      iso: string indicating an "isospin component of a density"
    Output:
      index : 1-4 indicating the index in the FORTRAN array
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      
  """

  assert len(iso) == 1
  
  index = 0
  
  if(iso == '0'):
    index = 3
  elif(iso == '1'):
    index = 4
  elif(iso == 'n'):
    index = 1
  elif(iso == 'p'):
    index = 2  
  return index


