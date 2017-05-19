from operatortospwf import *
from string         import Template
import numpy        as np
from math           import factorial

densities =[]
dimensions=[]
ArrayNames=[]
matsizes  =[]

def initdensities():
    global densities, dimensions, MATSIZE, ArrayNames
    
    densities = ['rho']
    order     = [  0  ]
    
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

    
  
#   for den in densities:

    #-----------------------------------------------------------------------
    # Fixing correct calculation
    (Expression, Declaration, Initialisation) = GenDensityExpression(Identity, Identity, ['rho', 'vecs'], 0)
    
    (E,D,I) = GenDensityExpression(Nabla, Nabla, ['tau', 'TN2LO'], 0)
    
    Expression     = Expression     + '\n '  + E
    Declaration    = Declaration    + '\n '  + D
    Initialisation = Initialisation + '\n '  + I
        
    dic={}
    dic['DECLARATION']    = Declaration
    dic['INITIALIZATION'] = Initialisation
    dic['EXPRESSION']     = Expression    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))  

def GenDensityExpression(LeftOperator, RightOperator, Names, Current):
    #---------------------------------------------------------------------------
    # This function generates three strings, needed for the calculation of 
    # densities, starting from (already computed) single-particle wave-functions
    # and their derivatives. 
    #
    # Expression    = expression of the density as a function of the spwfs
    # Declaration   = FORTRAN declaration of the dimension of the density
    # Initialisation= Setting the density to 0
    
    Expression     = ''
    Initialisation = ''
    Declaration    = ''
    #---------------------------------------------------------------------------
    #Start from standard wave-functions, [ Psi_1, Psi_2, Psi_3, Psi_4 ]^T
    
    leftind  = np.zeros((4,sum(matsizes)))
    rightind = np.zeros((4,sum(matsizes)))    

    leftind[0,0] = 1 ; rightind[0,0] = 1
    leftind[1,0] = 2 ; rightind[1,0] = 2 
    leftind[2,0] = 3 ; rightind[2,0] = 3
    leftind[3,0] = 4 ; rightind[3,0] = 4
    
    #---------------------------------------------------------------------------
    # Act with the derivative operators on the left and right
    leftind   = LeftOperator(leftind)
    rightind  = RightOperator(rightind)
    
    #---------------------------------------------------------------------------
    # If Current == True then we modify the indices to reflect the imaginary 
    # spinor product
    if(Current == 1):
        leftind   = MakeCurrent(leftind)
    if(Current != 0 and Current != 1):
        print 'Hephaestos can not yet deal with currents of currents.' 
        exit
    #---------------------------------------------------------------------------
    # Determine if the result is scalar, vector or tensor in derivatives.
    leftorder  = determine_order(leftind) 
    rightorder = determine_order(rightind)

    #---------------------------------------------------------------------------
    # Density calculation template to fill in
    # Could be defined globally, but is nice to have here for quick reference
    Den_template_1 = Template('$NAME(i,j,k$IND,it) = $NAME(i,j,k$IND,it) + $WEIGHT * (')
    Den_template_2 = Template('$SIGN $LEFTWF(i,j,k,$LIND,wave) * $RIGHTWF(i,j,k,$RIND,wave)')
    
    Ini_template   = Template( '    if(.not.allocated($NAME)) then     \n' + \
                               '       allocate($NAME(nx,ny,nz$DIM,2)) \n' + \
                               '    endif \n'                              + \
                               '    $NAME = 0.0d0')
    Dec_template   = Template('real(KIND=dp),allocatable :: $NAME(:,:,:$TOTALIND,:)')

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
    # Now add the corresponding statements
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
            
    #---------------------------------------------------------------------------      
    # Now the spin-vector
    Expression    = Expression + '\n'
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
    # Now add the corresponding statements
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
                rightind = Sigma(k+1,old_right_ind)
                
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
            
    return (Expression, Declaration, Initialisation)

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
#-------------------------------------------------------------------------------
def Identity(indices):
    # Operates on the indices of the wave-function
    # Identity operator
    return indices
    
def Nabla(indices):
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

def Sigma(Dir,indices):
    # Operates on the indices of the wave-function
    # \sigma_{Dir} operator
    out = np.zeros_like(indices)
    if(Dir == 1 or Dir == 'x'):
        out[0,:] =   indices[2,:]
        out[1,:] =   indices[3,:]
        out[2,:] =   indices[0,:]
        out[3,:] =   indices[1,:] 
        return out
    elif(Dir == 2 or Dir == 'y'):
        out[0,:] = - indices[3,:]
        out[1,:] =   indices[2,:]
        out[2,:] =   indices[1,:]
        out[3,:] = - indices[0,:] 
        return out
    elif(Dir == 3 or Dir == 'z'):
        out[0,:] =   indices[0,:]
        out[1,:] =   indices[1,:]
        out[2,:] = - indices[2,:]
        out[3,:] = - indices[3,:] 
        return out

def MakeCurrent(indices):
    # Operates on indices to get a current C, instead of a density D
    out = np.zeros_like(indices)

    out[0,:] =   indices[1,:] 
    out[1,:] = - indices[0,:]
    out[2,:] =   indices[3,:]
    out[3,:] = - indices[2,:]
    
    return out
    

