from operatortospwf import *
from string         import Template
import numpy        as np
from math           import factorial

densities     =[]
dimensions    =[]
ArrayNames    =[]
matsizes      =[]
leftoperators =[]
rightoperators=[]
currents      =[]
deriv_needed  =[]
lapla_needed  =[]
# Tab-character. 4 spaces for W.R., but I can imagine other people have
# different standards.
tab       ='    '

def initdensities():
    global densities, dimensions, MATSIZE, leftoperators, rightoperators
    global currents, deriv_needed, lapla_needed, ArrayNames
    
#    densities     = [['rho','vecs'], ['tau','TN2LO'], ['vecj', 'Jmunu'], ['tau_c','TN2LO_C'] ]
#    leftoperators = [Identity, Nabla, Identity,  Nabla]
#    rightoperators= [Identity, Nabla, Nabla   ,  Nabla]
#    currents      = [0,0,1,1,1]
#    deriv_needed  = [1,0,0,0,0]
#    lapla_needed  = [1,0,0,0,0]

    densities     = [['tau_c','TN2LO_C'] ]
    leftoperators = [ Nabla]
    rightoperators= [ Nabla]
    currents      = [1]
    deriv_needed  = [0]
    lapla_needed  = [0]
    #---------------------------------------------------------------------------
    # This defines the maximum of derivatives of the single-particle 
    # wave-functions that will be needed in the code. 
    # For reference:
    # LO   densities:  0
    # NLO  densities:  1
    # N2LO densities:  2
    # N3LO densities:  3
    #---------------------------------------------------------------------------
    MAXDERIV = 2
    #---------------------------------------------------------------------------
    # Arraynames
    #
    # Defines the names of the arrays in the fortran code that contain the 
    # relevant derivative of the single-particle wave-functions.
    #
    # Note that the script assumes 'upper-diagonal' storage in the sense that 
    #     HFdPsi(:,:,:,1,2,wave) contains           \partial_xy \Psi
    # and HFdPsi(:,:,:,1,2,wave) is not initialized
    #
    ArrayNames=['HFPsi', 'HFdPsi', 'HFddPsi']

    #---------------------------------------------------------------------------
    # The l'th order of derivatives adds (2*l+1) possible combinations
    # Note that the total number of possible non-redundant entries of a symmetric
    # tensor of rank n over a space of dimension d is
    # 
    # ( n + d - 1)!
    # -------------
    # n! (d-1)!
    #---------------------------------------------------------------------------    
    for l in range(MAXDERIV+1):
        matsizes.append(factorial(l + 3 -1)/(factorial(l)*2))

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
        
        # Generating all of the strings
        (E,D, I, Der) = GenDensityExpression( leftoperators[i],rightoperators[i],densities[i], currents[i], deriv_needed[i], lapla_needed[i])

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

def GenDensityExpression(LeftOperator, RightOperator, Names, Current, Der, Lap):
    #---------------------------------------------------------------------------
    # This function generates three strings, needed for the calculation of 
    # densities, starting from (already computed) single-particle wave-functions
    # and their derivatives. 
    #
    # Expression    = expression of the density as a function of the spwfs
    # Declaration   = FORTRAN declaration of the dimension of the density
    # Initialisation= Setting the density to 0
    # Derivation    = calculation of derivatives of this density
    #
    # Current determines whether (Current = 1) or not (Current = 0) we are 
    # dealing with a normal density or with a current. 
    #
    # Der keeps track if we need derivatives of the density in the program or
    # not. (Divergences and curls can later be calculated too.)
    # 
    # Der
    # 0    No gradient needed
    # 1    Gradient needed
    #
    # Lap keeps track if we need laplacians of this density or not.
    #
    # Lap
    # 0    No laplacian needed
    # 1    Laplacian of density needed
    # 2    Laplacian of laplacian of density needed
    #---------------------------------------------------------------------------
    global tab
    
    Expression     = ''
    Initialisation = ''
    Declaration    = ''
    Derivation     = ''
    #---------------------------------------------------------------------------
    #Start from standard wave-functions, [ Psi_1, Psi_2, Psi_3, Psi_4 ]^T
    
    leftind  = np.zeros((4,sum(matsizes)))
    rightind = np.zeros((4,sum(matsizes)))    

    leftind[0,0] = 1 ; rightind[0,0] = 1
    leftind[1,0] = 2 ; rightind[1,0] = 2 
    leftind[2,0] = 3 ; rightind[2,0] = 3
    leftind[3,0] = 4 ; rightind[3,0] = 4
    
    #---------------------------------------------------------------------------
    # If Current == True then we modify the indices to reflect the imaginary 
    # spinor product
    if(Current == 1):
        New_RightOperator = Combine( MakeCurrent, RightOperator)
        
    if(Current != 0 and Current != 1):
        print 'Hephaestos can not yet deal with currents of currents.' 
        exit

    
    #---------------------------------------------------------------------------
    # Act with the derivative operators on the left and right
    leftind   = LeftOperator(leftind)
    rightind  = RightOperator(rightind)
    
    #---------------------------------------------------------------------------
    # Determine if the result is scalar, vector or tensor in derivatives.
    leftorder  = determine_order(leftind) 
    rightorder = determine_order(rightind)

    #---------------------------------------------------------------------------
    # Density calculation template to fill in
    # Could be defined globally, but is nice to have here for quick reference
    Den_template_1 = Template(tab+'$NAME(i,j,k$IND,it) = $NAME(i,j,k$IND,it) + $WEIGHT * (')
    Den_template_2 = Template(tab+'$SIGN $LEFTWF(i,j,k,$LIND,wave) * $RIGHTWF(i,j,k,$RIND,wave)')
    
    Ini_template   = Template( tab+'if(.not.allocated($NAME)) then     \n' + \
                               tab+tab+'allocate($NAME(nx,ny,nz$DIM,2)) \n' + \
                               tab+'endif \n'                              + \
                               tab+'$NAME = 0.0d0')
    Dec_template   = Template( tab + 'real(KIND=dp),allocatable :: $NAME(:,:,:$TOTALIND,:)')
    Der_template   = Template( tab + tab +'call Derive_grad($NAME(:,:,:$IND,it), $PX,$PY,$PZ, der_$NAME(:,:,:$IND,1,it),der_$NAME(:,:,:$IND,2,it),der_$NAME(:,:,:$IND,3,it))')
    Lap_template   = Template( tab + tab +'call Derive_lap ($NAME(:,:,:$IND,it), $PX,$PY,$PZ, lap_$NAME(:,:,:$IND,it))')
    
    
    dic={}
    dic['LEFTWF']  = ArrayNames[leftorder]
    dic['RIGHTWF'] = ArrayNames[rightorder]
    dic['WEIGHT']    = 'occupations(wave)' # For now only simply the occupations

    #---------------------------------------------------------------------------
    #First the spin-scalar
    dic['NAME']   = Names[0]

    size_left = matsizes[leftorder]        
    size_right= matsizes[rightorder]     
    
    dic['TOTALIND']= ''
    dic['DIM']     = ''
    for ls in range(leftorder):
              dic['TOTALIND']= dic['TOTALIND'] + ',:' 
              dic['DIM']     = dic['DIM'] + ',3' 
    for rs in range(rightorder):
              dic['TOTALIND']= dic['TOTALIND'] + ',:'
              dic['DIM']     = dic['DIM'] + ',3' 
    
    Declaration    = Dec_template.substitute(dic)
    Initialisation = Ini_template.substitute(dic)
    
    #---------------------------------------------------------------------------   
    # Add the derivatives of the original density
    if(Lap == 1):
        dic['NAME']    = 'Lap_' + Names[0]
        
        Declaration    = Declaration    + '\n' + Dec_template.substitute(dic)
        Initialisation = Initialisation + '\n' + Ini_template.substitute(dic)
        dic['NAME']    = Names[0]
    if(Der == 1):
        dic['NAME']    = 'Der_' + Names[0]
        dic['TOTALIND']= dic['TOTALIND'] + ',:' 
        dic['DIM']     = dic['DIM'] + ',3' 
    
        Declaration    = Declaration    + '\n' + Dec_template.substitute(dic)
        Initialisation = Initialisation + '\n' + Ini_template.substitute(dic)
        dic['NAME']    = Names[0]
    #---------------------------------------------------------------------------
    # Now add the calculation statement
    for ls in range(size_left):
        for rs in range(ls, size_right):
            
            #-------------------------------------------------------------------
            # First declaration statement
            IND           = ''
            if(leftorder != 0):
                IND = IND + ',' + str(ls+1)
            if(rightorder != 0):
                IND = IND + ',' + str(rs+1)
            
            dic['IND']     = IND
                        
            Expression = Expression +  Den_template_1.substitute(dic)

            leftcolumn = sum(matsizes[0:leftorder])  + ls
            rightcolumn= sum(matsizes[0:rightorder]) + rs - ls
            
            # Loop over the wave-function components            
            for i in range(4):
                LIND=str(int(leftind[i,leftcolumn]))
                
                for l in range(leftorder):
                    LIND = LIND +  ',' + str(ls+1)  
                dic['LIND']   = LIND
                
                RIND=str(int(rightind[i,rightcolumn]))
                for r in range(rightorder):
                    RIND = RIND + ',' + str(rs+1) 
                dic['RIND']   = RIND
                
                SIGN          = np.sign(leftind[i,leftcolumn])*np.sign(rightind[i,rightcolumn])
                if(SIGN > 0) :
                    dic['SIGN']   = '+'
                else :
                    dic['SIGN']   = '-'
                Expression = Expression +  \
                            '& \n               &' +  \
                            Den_template_2.substitute(dic)
                            
            # Don't forget the closing bracket
            Expression = Expression + ') \n'
            #-------------------------------------------------------------------  
            # Add the derivatives of the original density
            if(Lap == 1):

                (px,py,pz)     = AxisReflection(LeftOperator, RightOperator,ls,rs)
                
                dic['PX']      = str(px)
                dic['PY']      = str(py)
                dic['PZ']      = str(pz)
                
                dic['NAME']    = Names[0]
                Derivation     = Derivation + '\n' + Lap_template.substitute(dic)
        
            if(Der == 1):
            
                (px,py,pz)     = AxisReflection(LeftOperator, RightOperator,ls,rs)
                
                dic['PX']      = px
                dic['PY']      = py
                dic['PZ']      = pz
                
                dic['NAME']    = Names[0]
                Derivation     = Derivation + '\n' + Der_template.substitute(dic)
                    
    #---------------------------------------------------------------------------      
    # Now the spin-vector
    Expression    = Expression + '\n'
    Derivation    = Derivation + '\n'
    dic['NAME']   = Names[1]

    size_left = matsizes[leftorder]        
    size_right= matsizes[rightorder]  
    
    old_right_ind = rightind
    
    dic['TOTALIND']= ',:'
    dic['DIM']     = ',3'
    for ls in range(leftorder):
              dic['TOTALIND']= dic['TOTALIND'] + ',:' 
              dic['DIM']     = dic['DIM'] + ',3' 
    for rs in range(rightorder):
              dic['TOTALIND']= dic['TOTALIND'] + ',:'
              dic['DIM']     = dic['DIM'] + ',3' 
              
    Declaration    = Declaration    + '\n' + Dec_template.substitute(dic)
    Initialisation = Initialisation + '\n' + Ini_template.substitute(dic)
    #---------------------------------------------------------------------------   
    # Add the derivatives to the spin density
    if(Lap == 1):
        dic['NAME']    = 'Lap_' + Names[1]
                    
        Declaration    = Declaration    + '\n' + Dec_template.substitute(dic)
        Initialisation = Initialisation + '\n' + Ini_template.substitute(dic)
        dic['NAME']    = Names[1]
    if(Der == 1):
        dic['NAME']    = 'Der_' + Names[1]
        dic['TOTALIND']= dic['TOTALIND'] + ',:' 
        dic['DIM']     = dic['DIM'] + ',3' 
                
        Declaration    = Declaration    + '\n' + Dec_template.substitute(dic)
        Initialisation = Initialisation + '\n' + Ini_template.substitute(dic)
        dic['NAME']    = Names[1]
    
    #---------------------------------------------------------------------------
    # Now add the corresponding statements
    
    Sigma = [Sigma_x, Sigma_y, Sigma_z]
    for k in range(3): 
        for ls in range(size_left):
            for rs in range(ls, size_right):
                # Spatial indices
                IND           = ''
                if(leftorder != 0):
                    IND = IND + ',' + str(ls+1)
                if(rightorder != 0):
                    IND = IND + ',' + str(rs+1)
                    
                # Add the spin index
                IND  = IND + ',' + str(k+1)
                
                dic['IND']    = IND
                
                Expression = Expression + Den_template_1.substitute(dic)

                leftcolumn = sum(matsizes[0:leftorder])  + ls
                rightcolumn= sum(matsizes[0:rightorder]) + rs - ls
                
                #Apply sigma on the right
                rightind = Sigma[k](old_right_ind)
                
                # Loop over the wave-function components            
                for i in range(4):
                    LIND=str(abs(int(leftind[i,leftcolumn])))
                    
                    for l in range(leftorder):
                        LIND = LIND +  ',' + str(ls+1)  
                    dic['LIND']   = LIND
                    
                    RIND=str(abs(int(rightind[i,rightcolumn])))
                    for r in range(rightorder):
                        RIND = RIND + ',' + str(rs+1) 
                    dic['RIND']   = RIND
                    
                    SIGN          = np.sign(leftind[i,leftcolumn])*np.sign(rightind[i,rightcolumn])
                    if(SIGN > 0) :
                        dic['SIGN']   = '+'
                    else :
                        dic['SIGN']   = '-'
                    Expression = Expression +  \
                                '& \n               &' +  \
                                Den_template_2.substitute(dic)
                # Don't forget the closing bracket
                Expression = Expression + ') \n'
                
                
                #-------------------------------------------------------------------  
                # Add the derivatives of the original density
                if(Lap == 1):

                    New_right = Combine(RightOperator, Sigma[k])
                    
                    (px,py,pz)     = AxisReflection(LeftOperator,New_right,ls,rs)
                    
                    dic['PX']      = str(px)
                    dic['PY']      = str(py)
                    dic['PZ']      = str(pz)
                    
                    Derivation     = Derivation + '\n' + Lap_template.substitute(dic)
            
                if(Der == 1):
                
                    New_right = Combine(RightOperator, Sigma[k])
                    (px,py,pz)     = AxisReflection(LeftOperator, New_right,ls,rs)
                    
                    dic['PX']      = px
                    dic['PY']      = py
                    dic['PZ']      = pz
                    
                    Derivation     = Derivation + '\n' + Der_template.substitute(dic)
             
#            
    return (Expression, Declaration, Initialisation, Derivation)

def determine_order(indices):
    #---------------------------------------------------------------------------
    # Calculate if a given set of indices represents a scalar, vector or tensor
    # Currently idiotically implemented, can be generalized.
    for i in range(len(indices[0,:])):
        if(indices[0,i] != 0) :
            if(i == 0):
                return 0
            if(i == 1):
                return 1
            if(i == 4):
                return 2
#-------------------------------------------------------------------------------
# Definition of operators for use in the densities.
#
# Note that python allows the functions to carry attributes. In our case this
# is used to store the behaviour under the symmetry operators. 
#-------------------------------------------------------------------------------
def Identity(indices):
    # Operates on the indices of the wave-function
    # Identity operator
    
    Identity.parity      =np.array([1])
    Identity.time        =np.array([1])
    Identity.signature_x =np.array([1])   
    Identity.signature_y =np.array([1])
    Identity.signature_z =np.array([1])
    
    return indices
    
def Nabla(indices):

    Nabla.parity      =np.array([-1,-1,-1])
    Nabla.time        =np.array([ 1, 1, 1])
    Nabla.signature_x =np.array([ 1,-1,-1])   
    Nabla.signature_y =np.array([-1, 1,-1])
    Nabla.signature_z =np.array([-1,-1, 1])   

    # Operates on the indices of the wave-function
    out = np.zeros_like(indices)
    
    # First order derivative of an ordinary wave-function
    out[:,1] = indices[:,0]
    out[:,2] = indices[:,0]
    out[:,3] = indices[:,0]
    
    # First order derivative of a first derivative
    out[:,4] = indices[:,1] 
    out[:,5] = indices[:,2] # dxy is calculated as dx(dy)
    out[:,7] = indices[:,3] # dxz is calculated as dx(dz)
    out[:,6] = indices[:,2] # dyy is calculated as dy(dy)
    out[:,8] = indices[:,3] # dyz is calculated as dy(dz)
    out[:,9] = indices[:,3] # dzz is calculated as dz(dz)
    
    return out

def NablaNabla(indices):
    #                            xx xy xz yy yz zz
    NablaNabla.parity      =np.array([+1,+1,+1,+1,+1,+1])
    NablaNabla.time        =np.array([ 1, 1, 1, 1, 1, 1])
    NablaNabla.signature_x =np.array([ 1,-1,-1, 1, 1, 1])   
    NablaNabla.signature_y =np.array([ 1,-1, 1, 1,-1, 1])   
    NablaNabla.signature_z =np.array([ 1, 1,-1, 1,-1, 1])      

    # Operates on the indices of the wave-function
    out = np.zeros_like(indices)
    
    # First order derivative of an ordinary wave-function
    out[:,1] = indices[:,0]
    out[:,2] = indices[:,0]
    out[:,3] = indices[:,0]
    
    # First order derivative of a first derivative
    out[:,4] = indices[:,1] 
    out[:,5] = indices[:,2] # dxy is calculated as dx(dy)
    out[:,7] = indices[:,3] # dxz is calculated as dx(dz)
    out[:,6] = indices[:,2] # dyy is calculated as dy(dy)
    out[:,8] = indices[:,3] # dyz is calculated as dy(dz)
    out[:,9] = indices[:,3] # dzz is calculated as dz(dz)
    
    return out


def Sigma_x(indices):
    # Operates on the indices of the wave-function
    # \sigma_x operator
    
    Sigma_x.parity      =np.array([ 1])
    Sigma_x.time        =np.array([-1])
    Sigma_x.signature_x =np.array([ 1])   
    Sigma_x.signature_y =np.array([-1])
    Sigma_x.signature_z =np.array([-1])
    
    out = np.zeros_like(indices)
    out[0,:] =   indices[2,:]
    out[1,:] =   indices[3,:]
    out[2,:] =   indices[0,:]
    out[3,:] =   indices[1,:] 
    return out

def Sigma_y(indices):
    # Operates on the indices of the wave-function
    # \sigma_y operator
    out = np.zeros_like(indices)
    
    Sigma_y.parity      =np.array([ 1])
    Sigma_y.time        =np.array([-1])
    Sigma_y.signature_x =np.array([-1])   
    Sigma_y.signature_y =np.array([ 1])
    Sigma_y.signature_z =np.array([-1])
    
    out[0,:] = - indices[3,:]
    out[1,:] =   indices[2,:]
    out[2,:] =   indices[1,:]
    out[3,:] = - indices[0,:] 
    return out
    
def Sigma_z(indices):
    # Operates on the indices of the wave-function
    # \sigma_z operator
    out = np.zeros_like(indices)
    
    Sigma_z.parity      =np.array([ 1])
    Sigma_z.time        =np.array([-1])
    Sigma_z.signature_x =np.array([-1])   
    Sigma_z.signature_y =np.array([-1])
    Sigma_z.signature_z =np.array([ 1])
    
    out[0,:] =   indices[0,:]
    out[1,:] =   indices[1,:]
    out[2,:] = - indices[2,:]
    out[3,:] = - indices[3,:] 
    return out

def MakeCurrent(indices):
    # Operates on indices to get a current C, instead of a density D
    out = np.zeros_like(indices)

    MakeCurrent.parity      =np.array([ 1])
    MakeCurrent.time        =np.array([-1])
    MakeCurrent.signature_x =np.array([ 1])   
    MakeCurrent.signature_y =np.array([ 1])
    MakeCurrent.signature_z =np.array([ 1])

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
    
    LR = lambda x:L(R(x))
    
    LR.parity      = R.parity      * L.parity
    LR.signature_x = R.signature_x * L.signature_x
    LR.signature_y = R.signature_y * L.signature_y
    LR.signature_z = R.signature_z * L.signature_z
    
    return LR
    
#-------------------------------------------------------------------------------
# Define behaviour under reflection of an axis under symmetry operators.
# Currently works for EV8/CR8-like calculations.
#-------------------------------------------------------------------------------

def AxisReflection(LeftOperator, RightOperator, leftorder, rightorder):
    #---------------------------------------------------------------------------
    #
    # Calculates the sign under axis reflection for the symmetries of an EV8/CR8
    # calculation.
    #
    # leftorder = index of the left operator
    #
    # rightorder= index of the right operator
    #---------------------------------------------------------------------------
    pxl = LeftOperator.signature_x[leftorder]  * LeftOperator.parity[leftorder]
    pxr = RightOperator.signature_x[rightorder] * RightOperator.parity[rightorder]
    px  = pxl * pxr
    
    if(px > 0):
        px = '+1'
    else:
        px = '-1'
    
    pyl = LeftOperator.signature_y[leftorder]  * LeftOperator.parity[leftorder]
    pyr = RightOperator.signature_y[rightorder] * RightOperator.parity[rightorder]
    py  = pyl * pyr
    
    if(py > 0):
        py = '+1'
    else:
        py = '-1'
    
    pzl = LeftOperator.signature_z[leftorder]  * LeftOperator.parity[leftorder]
    pzr = RightOperator.signature_z[rightorder] * RightOperator.parity[rightorder]
    pz  = pzl * pzr
    
    if(pz > 0):
        pz = '+1'
    else:
        pz = '-1'
    
    return(px,py,pz)
