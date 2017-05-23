from operatortospwf import *
from string         import Template
import numpy        as np
from math           import log
import itertools

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

#-------------------------------------------------------------------------------
# Tab-character for the fortran code.
# 4 spaces for W.R., but I can imagine other people have different standards.
tab       ='    '


def initdensities():
    global densities, dimensions, leftoperators, rightoperators
    global currents, deriv_needed, lapla_needed, ArrayNames

    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Initialize the predefined operators with correct atributes, before they
    # ever get called.
        
    Identity.dimension   =0
    Identity.derorder    =0  
    Identity.parity      =np.array([1])
    Identity.time        =np.array([1])
    Identity.signature_x =np.array([1])   
    Identity.signature_y =np.array([1])
    Identity.signature_z =np.array([1])
    
    Current.derorder     =0
    Current.dimension    =0
    Current.parity      =np.array([1])
    Current.time        =np.array([1])
    Current.signature_x =np.array([1])   
    Current.signature_y =np.array([1])
    Current.signature_z =np.array([1])
    
    Nabla.dimension      =1
    Nabla.derorder       =1
    Nabla.parity      =np.array([-1,-1,-1])
    Nabla.time        =np.array([ 1, 1, 1])
    Nabla.signature_x =np.array([ 1,-1,-1])   
    Nabla.signature_y =np.array([-1, 1,-1])
    Nabla.signature_z =np.array([-1,-1, 1])   
    
    Sigma.dimension      =1 
    Sigma.derorder       =0
    Sigma.parity      =np.array([ 1, 1, 1])
    Sigma.time        =np.array([-1,-1,-1])
    Sigma.signature_x =np.array([ 1,-1,-1])   
    Sigma.signature_y =np.array([-1, 1,-1])
    Sigma.signature_z =np.array([-1,-1, 1])
   
    
    ArrayNames=['HFPsi', 'HFdPsi', 'HFddPsi']

    NS  = Combine(Nabla, Sigma)
    CN  = Combine(Current,Nabla)
    NNS = Combine(     NS,Nabla)

    densities     = ['rho', 'vecs', 'tau'      , 'TN2LO', 'vecj']
    leftoperators = [Identity, Identity, Nabla , Nabla  , Identity]
    rightoperators= [Identity, Sigma, Nabla    , NS     , CN]

def processdensities(fname, src, target):

    #===========================================================================
    # Writing the code to compute the densities
    #
    #===========================================================================

    global densities, dimensions, MATSIZE, leftoperators, rightoperators
    global currents, deriv_needed, lapla_needed

    Expression     = ''
    Declaration    = ''
    Initialisation = ''
    Derivation     = ''

    #-----------------------------------------------------------------------
    # Fixing correct calculation
    for i in range(len(densities)): 
        
        print i, densities[i]
        
        # Generating all of the strings
        (E,D,I,Der) = GenDensityExpression( leftoperators[i],rightoperators[i],densities[i],0,0)

        # Appending
        Expression     = Expression     + '\n'  + E
        Declaration    = Declaration    + '\n'  + D
        Initialisation = Initialisation + '\n'  + I
        Derivation     = Derivation     + '\n'  + Der   
        
    dic={}
    dic['DECLARATION']    = Declaration
    dic['INITIALIZATION'] = Initialisation
    dic['EXPRESSION']     = Expression
    dic['DERIVATION']     = Derivation    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))  
    
def GenDensityExpression(LeftOperator, RightOperator, Name, Der, Lap):
    #---------------------------------------------------------------------------
    #
    
    Expression    = ''
    Declaration   = ''
    Initialisation= ''
    Derivation    = ''
    #---------------------------------------------------------------------------
    # Density calculation template to fill in
    # Could be defined globally, but is nice to have here for quick reference
    Den_template_1 = Template(tab+'$NAME(i$IND,it) = $NAME(i$IND,it) + $WEIGHT * (')
    Den_template_2 = Template(tab+'$SIGN $LEFTWF(i,1,1,$LCOMP$LIND,wave) * $RIGHTWF(i,1,1,$RCOMP$RIND,wave)')
    
    Ini_template   = Template( tab+'if(.not.allocated($NAME)) then     \n' + \
                               tab+tab+'allocate($NAME(mv$DIM,2)) \n' + \
                               tab+'endif \n'                              + \
                               tab+'$NAME = 0.0d0')
    Dec_template   = Template( tab + 'real*8,allocatable :: $NAME(:$TOTALIND,:)')
    Der_template   = Template( tab + tab +'call Derive_grad($NAME(:$IND,it), $PX,$PY,$PZ, der_$NAME(:$IND,1,it),der_$NAME(:$IND,2,it),der_$NAME(:$IND,3,it))')
    Lap_template   = Template( tab + tab +'call Derive_lap ($NAME(:$IND,it), $PX,$PY,$PZ, lap_$NAME(:$IND,it))')
    #---------------------------------------------------------------------------
    # Declaration and initialisation
    
    dic= {}
    dic['NAME']    = Name
    
    dic['LEFTWF']  = ArrayNames[ LeftOperator.derorder]
    dic['RIGHTWF'] = ArrayNames[RightOperator.derorder]
    dic['WEIGHT']  = 'occupations(wave)' # For now only simply the occupations
    
    totalind= ''
    dim     = ''
    
    for ls in range(LeftOperator.dimension):
        totalind = totalind + ',:' 
        dim      = dim      + ',3'
    for rs in range(RightOperator.dimension):
        totalind = totalind + ',:' 
        dim      = dim      + ',3'
    
    dic['TOTALIND']= totalind
    dic['DIM']     = dim
    
    Declaration    = Dec_template.substitute(dic)
    Initialisation = Ini_template.substitute(dic)
    
    #---------------------------------------------------------------------------
    # Generate the expression to calculate the density

    #Start from standard wave-functions, [ Psi_1, Psi_2, Psi_3, Psi_4 ]^T
    start  = np.zeros((4,1))
    start[0,0] = 1 
    start[1,0] = 2  
    start[2,0] = 3 
    start[3,0] = 4     
    # Construct an iterator with all possible combinations of indices
    argcomb_right = itertools.combinations_with_replacement(range(3), RightOperator.dimension) 
    # Now construct the calculation expression for every value of indices
    for rarg in argcomb_right:
        # Need to restart the iterator for every rarg
        argcomb_left  = itertools.combinations_with_replacement(range(3),  LeftOperator.dimension) 
        for larg in argcomb_left:
            leftind  = LeftOperator(larg, start)
            rightind = RightOperator(rarg, start)
            
            IND = ''
            for mu in larg: 
                IND = IND + ',' + str(mu+1) # Python indexes 0:N-1
            for nu in rarg: 
                IND = IND + ',' + str(nu+1)
            
            dic['IND']     = IND
            Expression = Expression +  Den_template_1.substitute(dic)
            
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
            
    return (Expression, Declaration, Initialisation, Derivation)

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
    
