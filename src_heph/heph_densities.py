#--------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#--------------------------------------------------------------------   
# 
# HOW TO USE
# =========== 
#
# The workhorse of this module is the function
#   GenDensityExpression(LeftOperator, RightOperator, Der, Lap)
#       
# That, when given a left-operator, right-operator for the density
# and which derivatives are needed of this density returns
#  a) and automatically generated name for the density
#  b) a Declaration expression, correctly declaring the density for a FORTRAN code
#  c) an Initialization expression, allocating and initialising the density
#  d) an Expression expression, calculating all of the components needed of the
#     density. 
#  e) A derivation expression, calculating all of the asked for derivatives of
#     this density. 
#
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
#
# Python module to treat densities.f90 file. Contains routines at the moment 
# that will likely migrate to more general modules as they are needed.
#
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
# work of this script, and the concious choice of not including currents of 
# currents is very important here.
#
# Four operators are currently defined in this file, with which we can construct
# all of the densities necessary.
#               Identity, Nabla, Sigma and Current
# All take a list of indices (even though Identity and Current don't need them) 
# and a 4-vector of components. As output, they permutate the components 
# according to their operator action and the index, with possible signs. 
# The action of sigma_y for example is coded as
#
#   Sigma(2, [1,2,3,4]) = [ -4, 3, 2, 1]
# representing  
#
#  sigma_y  ( Psi_1 + i Psi_2 ) = ( -Psi_4 + i Psi_3)
#           ( Psi_3 + i Psi_4 )   (  Psi_2 - i Psi_1)
#
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
#                          composite operator) are derivative indices, and not
#                          spin indices.
#   Operator.dimension= Order of the operator, scalar (0), vector (1) or tensor
#                       of rank( Operator.dimension)
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
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# 
# Important stuff to still do  (C = code, NC = noncode)
#
# In the near future
# ====================
#
#    C Support for detection of symmetric combinations 
#      (and subsequently not letting the FORTRAN code calculate them.)
#    C Support for packed storage schemes in FORTRAN: don't let FORTRAN store
#      all of the symmetric combinations
#      a) for the spwfs storage
#      b) for the density storage
#    C Support for laplacians inside the density 
#    C Currently for EV8-like symmetries, due to no decision yet about 
#      bookkeeping in Hephaestos itself. 
#    C Automatic continuation when time-reversal is conserved
#    C Add pairing densities (simply add another flag)
#    
#
#   NC Figure out how the densities should be named in the code. Historical
#      names have the advantage of lisibility, but cannot be extended to N2/3LO.
#      Then again things as D(dd) for tau are easily mistyped...
#   NC To be debated: do I want to create the concept of a densityvector again?
#
# In the distant future
# ======================
# 
#   C support for finite-range interactions: add another different position 
#     index. This is more FORTRAN work however.
#   
#-------------------------------------------------------------------------------

from string         import Template
from math           import log
import numpy        as np
import itertools
from heph_functional import Densities_needed, Deriv_den_needed, Lapla_den_needed, Densities_coupling

#-------------------------------------------------------------------------------
# Definition of lists needed by the preprocessing. Need to be global so that 
# Hephaestos can decide based on symmetries and parameterization asked which 
# ones to initialise.
#
#-------------------------------------------------------------------------------
densities     =[]
dimensions    =[]
ArrayNames    =[]
matsizes      =[]
leftoperators =[]
rightoperators=[]
currents      =[]
deriv_needed  =[]
lapla_needed  =[]
contractions  =[]
rotcoupl      =[]

#-------------------------------------------------------------------------------
# Tab-character for the fortran code.
# 4 spaces for W.R., but I can imagine other people have different standards.
tab       ='    '

def Rot_ind(k):
    # C_k = A_i B_j + A_j B_i
    # Input is k, output is the corresponding pair
    #     [ (i,j) ] if s == 0
    #     [ (j,i) ] if s == 1
    # where (i,j) carries the plus sign 
    
    if(k == 0):
        i = 1
        j = 2
    elif(k == 1):
        i = 2
        j = 0
    elif(k == 2):
        i = 0
        j = 1
        
    return (i,j)
    
def GenVecProd(density, coupling):
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Generate the appropriate expressions for the calculation of a vector
    # product. 
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    decl_template = Template(   tab + 'real(KIND=dp), allocatable :: $NAME(:$DECLIND,:)')
    ini_template  = Template(   tab + 'if(.not.allocated($NAME)) then \n'      \
                            + 2*tab + 'allocate($NAME(mv$ALLOCIND,2)) \n'      \
                            +   tab + 'endif  \n'                              \
                            +   tab + '$NAME = 0.0_dp  \n')
    calc_temp     = Template(   tab + '$NAME(:$ARG,:) = & \n')
    calc_a_temp   = Template( '$SIGN $DEN(:$IND,:) ')
    
    index_encountered = 0
    name              = ''
    ini               = ''
    
    dic = {}
    
    #---------------------------------------------------------------------------
    # Constructing the correct name: inserting xi, xj etc at the correct spot.
    for i in range(len(density)):
        l = density[i]
        if( density[i:i+3] == 'der'):
            for c in coupling: 
                if index_encountered in c:
                    name = name + l + 'x%d'%(coupling.index(c)+1)
                else :
                    name = name + l
            index_encountered=index_encountered+1
        elif( (l.isupper() and l != 'I' and l != 'D' and l!= 'C')):
            for c in coupling:
                if index_encountered in c:
                    name = name + l + 'x%d'%(coupling.index(c)+1)
            index_encountered=index_encountered+1
        else:   
            name = name + l
        
    dic['NAME'] = name
    
    
    #---------------------------------------------------------------------------
    # Constructing the correct set of indices for allocation and declaration.
    order      =  OrderOfDen(density)

    allocind   = ''
    for i in range(order - len(coupling)):
        allocind = allocind + ',3'
    
    dic['ALLOCIND'] = allocind
    dic['DECLIND']  = allocind.replace('3', ':')
    
    decl = decl_template.substitute(dic)
    ini  =  ini_template.substitute(dic)
    
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
                            p = p + (arg[ii],)
                            ii = ii +1
                uncontracted.append(p)
                
                
        dic['ARG'] = ''
        dic['DEN'] = density
        for i in range(len(arg)):
            dic['ARG'] = dic['ARG'] + ',%s'%(arg[i]+1)
        calc = calc + calc_temp.substitute(dic)
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
            calc = calc + calc_a_temp.substitute(dic)
            count = count + 1
            if(count%3 == 0):
                calc = calc + '& \n' + 2 * tab +'&'
        calc = calc  + '\n'       
            
    return (decl, ini, calc)
    
def OrderOfDen(density):
    #---------------------------------------------------------------------------
    # Returns the order (= number of indices) of the density represented 
    # by a string.Simply checks the number of capital letters N and S in the name
    # that are not contracted.
    
    # Add a dimension for every derivative and don't count the capital D or C
    test  = density.replace('_', '')
    order = test.count("der") - 1 
    test  = test.replace('der', '')
    test  = test.replace('Lap', '')
    
    # Add a dimension for every cross-product x1,x2, ...., x9
    for i in range(9):
        order = order + test.count("x%d"%i) 
        test  = test.replace('x%d'%i, '')
    
    for i in range(len(test)):
        letter = test[i]
        if(letter.isupper() and letter != 'I'):
            # Every capital letter that is not I adds an index
            order = order + 1
        if(not letter.isupper()):
            # But if the a non-capital letter follows, the index is contracted.
            order = order - 1
    return order

def initdensities():
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Initialize the predefined operators with correct atributes, before they
    # ever get called.
    global densities, dimensions, leftoperators, rightoperators
    global currents, deriv_needed, lapla_needed, ArrayNames, contractions
    global rotcoupl
   
        
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
    
    ArrayNames=['HFPsi', 'HFdPsi', 'HFddPsi']

    N   = Nabla
    I   = Identity
    NS  = Combine(Nabla, Sigma)
    CN  = Combine(Current,Nabla)
    CNS = Combine(Current,  NS)
    NNS = Combine(Nabla,    NS)
    NN  = Combine(Nabla, Nabla)

    deriv_needed = Deriv_den_needed
    lapla_needed = Lapla_den_needed
        

    for den in Densities_needed: 

        (x,y, left, right) =  ParseOperators(den)

        O = Identity
        for l in left:
                if(l == 'N'):
                        O = Combine(N, O)
                elif(l=='C'):
                        O = Combine(C, O)
                elif(l=='S'):
                        O = Combine(S, O)

        leftoperators.append(O)

        O = Identity
        for l in right:
                if(l == 'N'):
                        O = Combine(N, O)
                elif(l=='C'):
                        O = Combine(C, O)
                elif(l=='S'):
                        O = Combine(S, O)

        rightoperators.append(O)
        rotcoupl.append([])

def ProcessDensities(fname, src, target):
    #===========================================================================
    # Writing the code to compute the densities into densities.f90 file.
    #
    #===========================================================================

    global densities, dimensions, MATSIZE, leftoperators, rightoperators
    global currents, deriv_needed, lapla_needed, contractions, rotcoupl

    Expression     = ''
    Declaration    = ''
    Initialisation = ''
    Derivation     = ''
    VecProd        = ''
    DeclVecProd    = ''
    IniVecProd     = ''
    print ' - - - - - - - - - - - - - - - - - - - - - - - - - '
    print ' Generated densities '
    print ' - - - - - - - - - - - - - - - - - - - - - - - - - '
    print '         Name    RDim   CDim    Der    Contracted  '
    print ' - - - - - - - - - - - - - - - - - - - - - - - - - '
    for i in range(len(leftoperators)): 
        realorder= leftoperators[i].dimension + rightoperators[i].dimension
        order    = realorder - 2 * len(Densities_coupling[i])
        derorder = leftoperators[i].derorder + rightoperators[i].derorder
        
        # Getting all of the expression for all of the densities.
        (N,E,D,I,Der) = GenDensityExpression( leftoperators[i],rightoperators[i],deriv_needed[i],lapla_needed[i], Densities_coupling[i])

        print ' %12s %6d %6d %6d %6d %6d '%(N, realorder, order, derorder, deriv_needed[i], lapla_needed[i]), Densities_coupling[i]

        Expression     = Expression     + '\n'  + E
        Declaration    = Declaration    + '\n'  + D
        Initialisation = Initialisation + '\n'  + I
        if(len(Der)>0) :
             Derivation     = Derivation     + '\n'  + Der   
        if(len(rotcoupl[i])>0):
            print 'NAME', N
            (decl, ini, calc) = GenVecProd('der_' + N,rotcoupl[i])
            DeclVecProd = DeclVecProd + decl + '\n'
            IniVecProd  = IniVecProd  + ini  + '\n'
            VecProd     = VecProd     + calc + '\n'
    print ' - - - - - - - - - - - - - - - - - - - - - - - - - '

    
    # Substitute into the densities.f90 file.        
    dic={}
    dic['DECLARATION'   ] = Declaration
    dic['INITIALIZATION'] = Initialisation
    dic['EXPRESSION'    ] = Expression
    dic['DERIVATION'    ] = Derivation 
    dic['VECPROD'       ] = VecProd 
    dic['DECL_VECPROD'  ] = DeclVecProd
    dic['INIVECPROD'    ] = IniVecProd    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))  
    
def GenDensityExpression(LeftOperator, RightOperator, Der, Lap, Contract=[]):
    
    Expression    = ''
    Declaration   = ''
    Initialisation= ''
    Derivation    = ''
    #---------------------------------------------------------------------------
    # Density calculation template to fill in
    # Could be defined globally, but is nice to have here for quick reference
    Den_template_1 = Template(2*tab+'$NAME(i$IND,it) = $NAME(i$IND,it) + $WEIGHT * (')
    Den_template_2 = Template(tab+'$SIGN $LEFTWF(i,1,1,$LCOMP$LIND,wave) * $RIGHTWF(i,1,1,$RCOMP$RIND,wave)')
    
    Ini_template   = Template(   tab+'if(.not.allocated($NAME)) then     \n' + \
                               2*tab+'allocate($NAME(mv$DIM,2)) \n' + \
                                 tab+'endif \n'                              + \
                                 tab+'$NAME = 0.0d0')
    Dec_template   = Template(   tab + 'real*8,allocatable :: $NAME(:$TOTALIND,:)')
    Der_template   = Template( 2*tab +'call Derive_grad($NAME(:$IND,it),$PX,$PY,$PZ,der_$NAME(:,1$IND,it), &\n') 
    Der_template_b = Template( 2*tab + ' &  $DERSPACE der_$NAME(:,2$IND,it), &\n')
    Der_template_c = Template( 2*tab + ' &  $DERSPACE der_$NAME(:,3$IND,it))  \n')
    Lap_template   = Template( 2*tab +'call Derive_lap ($NAME(:$IND,it), $PX,$PY,$PZ, lap_$NAME(:$IND,it)) \n')
    
    #---------------------------------------------------------------------------
    # Some templates for comments to put into the densities file
    Den_comment          = Template(2*tab+'! Calculation of density $NAME \n')
    
    Den_comment_deriv    = Template(2*tab+'! Derivation of density $NAME  \n')
    Den_comment_deriv_b  = Template(2*tab+'! first order derivatives: $DER\n')
    Den_comment_deriv_c  = Template(2*tab+'! laplacians             : $LAP\n')
    
    Den_line             = Template(2*tab+'! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  \n')
   
    #---------------------------------------------------------------------------
    # Contraction indices to cycle through for the contractions
    cont_ind = ['m', 'n', 'k', 'r', 'q']
   
    #---------------------------------------------------------------------------
    # Find the correct name in our grand naming scheme
    autoname = 'D' 
    if( 'C' in LeftOperator.name  or  'C' in RightOperator.name):
        autoname = 'C'
    if( 'C' in LeftOperator.name  and 'C' in RightOperator.name):
        print 'HEPHAESTOS cannot manage currents of currents.'
        exit()
    
    left =LeftOperator.name.replace('C', '')
    right=RightOperator.name.replace('C', '')
    lname=''
    rname=''
    for l in range(len(left)):
        lname = lname + left[l]
        for c in Contract:
            if(l in c):
                lname = lname + cont_ind[Contract.index(c)]
    for r in range(len(right)):
        rname = rname + right[r]
        for c in Contract:
            if((r+len(left)) in c):
               rname = rname + cont_ind[Contract.index(c)]
            
    Name = autoname + '_' + lname + '_' + rname
   
    #---------------------------------------------------------------------------
    # Declaration and initialisation, also for the derivatives.
    dic= {}
    dic['NAME']    = Name
    
    dic['LEFTWF']  = ArrayNames[ LeftOperator.derorder]
    dic['RIGHTWF'] = ArrayNames[RightOperator.derorder]
    dic['WEIGHT']  = 'occupations(wave)' # For now only simply the occupations
    
    totalind= ''
    dim     = ''
    
    #---------------------------------------------------------------------------
    # Find the correct dimensions
    ndim = LeftOperator.dimension + RightOperator.dimension - 2*len(Contract)
    totalind=''
    dim     =''
    for r in range(ndim):
        totalind = totalind + ',:' 
        dim      = dim      + ',3'
    
    dic['TOTALIND']= totalind
    dic['DIM']     = dim
    
    Declaration    = Dec_template.substitute(dic)
    Initialisation = Ini_template.substitute(dic)
    
    for l in range(Lap):
        dic['NAME']    = (l+1)*'Lap_' + Name
        Declaration    = Declaration    + '\n' + Dec_template.substitute(dic)
        Initialisation = Initialisation + '\n' + Ini_template.substitute(dic)
        dic['NAME']    = Name
        
        
    for l in range(Der):
        dic['NAME']    = (l+1)*'Der_' + Name
        dic['TOTALIND']= dic['TOTALIND'] + ',:' 
        dic['DIM']     = dic['DIM']      + ',3' 
        Declaration    = Declaration    + '\n' + Dec_template.substitute(dic)
        Initialisation = Initialisation + '\n' + Ini_template.substitute(dic)
        dic['NAME']    = Name
    
    #---------------------------------------------------------------------------
    # Generate the expression to calculate the density
    # Start from standard wave-functions, [ Psi_1, Psi_2, Psi_3, Psi_4 ]^T
    start  = np.zeros((4,1))
    start[0,0] = 1 
    start[1,0] = 2  
    start[2,0] = 3 
    start[3,0] = 4     
    
    # Construct an iterator with all possible combinations of uncontracted indices
    args = itertools.product(range(3), repeat=ndim)
    
    Expression = Expression +  Den_line.substitute(dic)
    Expression = Expression +  Den_comment.substitute(dic)
    
    dic['DER'] = Der
    dic['LAP'] = Lap
    
    if(Der >=1 or Lap >= 1):
        Derivation = Derivation +  Den_line.substitute(dic)
        Derivation = Derivation +  Den_comment_deriv.substitute(dic)
        Derivation = Derivation +  Den_comment_deriv_b.substitute(dic)
        Derivation = Derivation +  Den_comment_deriv_c.substitute(dic)
    for arg in args:
        # We have the uncontracted indices. Now construct the combinations of
        # indices, including contracted ones, that correspond to this. 
        uncontracted = []
        if(len(Contract) == 0):
            uncontracted = [arg]
        else:
            cont = itertools.product(range(3), repeat=len(Contract))
            for c in cont:
                p  = ()   
                ii = 0
                for i in range(LeftOperator.dimension + RightOperator.dimension):
                        found = False                    
                        for combination in Contract:
                                if(i in combination): 
                                        p = p + (c[Contract.index(combination)],)
                                        found = True
                        if(not found): 
                                p = p + (arg[ii],)
                                ii = ii +1
                uncontracted.append(p)


        IND = ''
        for mu in arg: 
            IND = IND + ',' + str(mu+1) # Python indexes 0:N-1
        
        dic['IND'] = IND
        
        Expression = Expression +  Den_template_1.substitute(dic)
        
        # Now loop over the uncontracted indices
        for true_arg in uncontracted: 

            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            # Check if the indices for contraction are not superfluous
            #       
            larg = true_arg[:LeftOperator.dimension]
            rarg = true_arg[LeftOperator.dimension:]
            
            leftind  = LeftOperator(larg, start)
            rightind = RightOperator(rarg, start)
            
            lcolumns =  leftind.shape[1] 
            rcolumns = rightind.shape[1]
            
            # Note that the derivative operators have their indices on the left
            LIND = ''
            RIND = ''
            for lder in range(LeftOperator.derorder):
                LIND = LIND +  ',' + str(larg[lder]+1)
            for rder in range(RightOperator.derorder):
                RIND = RIND +  ',' + str(rarg[rder]+1)      
            dic['LIND']   = LIND
            dic['RIND']   = RIND
            
            for i in range(4):
                SIGN          = np.sign(leftind[i])*np.sign(rightind[i])
                if(SIGN > 0) :
                    dic['SIGN']   = '+'
                else :
                    dic['SIGN']   = '-'
                dic['RCOMP'] = str(int(abs(rightind[i,0])))
                dic['LCOMP'] = str(int(abs( leftind[i,0])))
                Expression = Expression +  \
                            '& \n               &' +  \
                            Den_template_2.substitute(dic)
                            
        # Don't forget the closing bracket
        Expression = Expression +  ')\n'
        
        #----------------------------------------------------------------------- 
        # Add the derivatives of the original density
        for l in range(Lap):
        
            #Preparing symmetries for derivatives
            (px,py,pz)   = AxisReflection(LeftOperator, RightOperator,larg,rarg)
            dic['PX']    = str(px)
            dic['PY']    = str(py)
            dic['PZ']    = str(pz)
        
            dic['NAME'] = l*'Lap_' + Name
            Derivation  = Derivation + Lap_template.substitute(dic)
            # Note that this is simple, since Laplacians don't change the 
            # symmetry properties of functions.
        
        for l in range(Der):
            dic['NAME'] = l*'der_' + Name
            
            # All components
            deriv_args = itertools.product(range(3), repeat=l)
            for darg in deriv_args: 
                dic['IND']  = ''    
                for i in darg:
                    dic['IND']  = dic['IND'] + ',%d'%(i+1)
                dic['IND']  = dic['IND'] +  IND    
                #Preparing symmetries for derivatives
                (px,py,pz)   = AxisReflection(LeftOperator, RightOperator,larg,rarg,darg)
                dic['PX']    = str(px)
                dic['PY']    = str(py)
                dic['PZ']    = str(pz)
            
                
                Derivation   = Derivation  + Der_template.substitute(dic)
                dic['DERSPACE'] = (18+len(Name)+ len(dic['IND']) +10) * ' '
                Derivation   = Derivation  + Der_template_b.substitute(dic)
                Derivation   = Derivation  + Der_template_c.substitute(dic)
            Derivation = Derivation + '\n'
        dic['NAME'] = Name
    Expression = Expression + Den_line.substitute(dic) 
    return (Name,Expression, Declaration, Initialisation, Derivation)

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
    
    if(mu[0] == 0):
        for i in range(columns):
            out[0,i] =   indices[2,i]
            out[1,i] =   indices[3,i]
            out[2,i] =   indices[0,i]
            out[3,i] =   indices[1,i] 
        return out
    elif(mu[0] == 1):
        for i in range(columns):
            out[0,i] = - indices[3,i]
            out[1,i] =   indices[2,i]
            out[2,i] =   indices[1,i]
            out[3,i] = - indices[0,i] 
        return out
    else:
        for i in range(columns):
            out[0,i] =   indices[0,i]
            out[1,i] =   indices[1,i]
            out[2,i] = - indices[2,i]
            out[3,i] = - indices[3,i] 
        return out

def Current(mu,indices):
    # Operates on indices to get a current C, instead of a density D
    out = np.zeros_like(indices)

    out[0,:] =   indices[1,:] 
    out[1,:] = - indices[0,:]
    out[2,:] =   indices[3,:]
    out[3,:] = - indices[2,:]
    
    return out

def Combine( L , R ):
    #---------------------------------------------------------------------------
    # Combines two symmetry-operators to define a new one. Note that the order
    # of the operators is important for non-Abelian systems.
    #
    # Note that this currently only works for operators of the same order.
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
    S = S.reshape(t)
    
    return S
    
def AxisReflection(LeftOperator, RightOperator, larg, rarg, nabla_arg = []):
    #---------------------------------------------------------------------------
    # Calculates the sign under axis reflection for the symmetries of an EV8/CR8
    # calculation.
    #
    # larg      = indices of the left operator
    # rarg      = indices of the right operator
    # nabla_arg = indices of the nabla_operators acting possibly on the density
    #
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

    # Symmetries of the left- and rightoperator
    pxl = LeftOperator.signature_x [mu] * LeftOperator.parity [mu]
    pxr = RightOperator.signature_x[nu] * RightOperator.parity[nu]
    px  = pxl * pxr

    pyl = LeftOperator.signature_y [mu] * LeftOperator.parity [mu] 
    pyr = RightOperator.signature_y[nu] * RightOperator.parity[nu] 
    py  = pyl * pyr
    
    pzl = LeftOperator.parity [mu]  * LeftOperator.signature_z [mu]
    pzr = RightOperator.parity[nu]  * RightOperator.signature_z[nu] 
    pz  = pzl * pzr
    
    #Taking into account extra nabla's
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
    else:
        px = '-1'
            
    if(py > 0):
        py = '+1'
    else:
        py = '-1'
    
    if(pz > 0):
        pz = '+1'
    else:
        pz = '-1'
    
    return(px,py,pz)

def ParseOperators(density):    
        #
        # Parse the operators that are used to construct a density from its name.
        #
        #
        left  = ''
        right = '' 

        der   = density.count('der_')
        lap   = density.count('Lap_')

        split = density.replace('der_', '').replace('Lap_', '').split('_')
        
        left = split[1]
        right= split[2]

        if(density[0] == 'C'):      
                right = 'C' + right

        return(der, lap, left, right)
    
