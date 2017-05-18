#===============================================================================
# Module that allows the user to automatically translate the action of an
# operator on a wave-function:
#
#  \hat{O} | \Psi \rangle or \langle  \Psi | \hat{O}
#  (O acting on the right)   (O acting on the left)
#
# into Fortran code, based on the storage scheme
#
#        ( Re Psi (s = +1) )
#  Psi = ( Im Psi (s = +1) )
#        ( Re Psi (s = -1) )
#        ( Im Psi (s = -1) )
# 
# Note that there is no reference to signature of the states, as in cr8. 
#
# This will (hopefully) be particularly useful in automatically writing 
# statements of the form:
#
#    Density(x,y,z) = Density(x,y,z) + &
# &                 weight * Psi_L(r,sigma) \hat{O}_R  \hat{O}_L Psi_R(r,sigma)
#
# where further indices can be easily accomodated, O_R is not necessarily equal
# to O_L, and they both act in their respective directions. 
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
#
#===============================================================================
# Examples 
#
#
#
#===============================================================================
from string import *



#def GenDensityExpression( LeftWF, leftn, RightWF, rightn, LeftOperator, RightOperator):
#    Expression = ''

#    #Start from standard wave-functions
#    leftind  = [1,2,3,4]
#    rightind = [1,2,3,4]

#    # Find the correct action of right-operator
#    rightind = RightOperator(rightind)
#    leftind  = LeftOperator (leftind)

#    # Determine the signs to put in the code
#    signs=['','','','']
#    for i in range(4):
#        if(leftind[i] * rightind[i] < 0) :
#            signs[i] = '-'
#        else :
#            signs[i] = '+'

#    # Fill in the expression
#    for i in range(4):
#        Expression = Expression   +  \
#                    '& \n               &' +  \
#                    '%s %s(i,1,1,%d,%s)*%s(i,1,1,%d,%s)'%(signs[i], \
#                                  LeftWF, abs( leftind[i]),  leftn, \
#                                 RightWF, abs(rightind[i]), rightn  )
#    return Expression
       
def Identity(indices):
    # Operates on the indices of the wave-function
    # Identity operator
    return [ indices[0], indices[1],  indices[2],  indices[3] ]
    
def Sigma(Dir,indices):
    # Operates on the indices of the wave-function
    # \sigma_{Dir} operator
    if(Dir == 1 or Dir == 'x'):
        return [ indices[2], indices[3],  indices[0],  indices[1]]
    elif(Dir == 2 or Dir == 'y'):
        return [-indices[3], indices[2],  indices[1], -indices[0]]
    elif(Dir == 3 or Dir == 'z'):
        return [ indices[0], indices[1], -indices[2], -indices[3]]
    
