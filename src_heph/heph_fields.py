#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# TODO
#   
#   Q Is the isospin coupling of the fields ok?
#   A NO, see Sadoudi.
#-------------------------------------------------------------------------------

from string import Template
import itertools
import numpy as np
from heph_densities    import *
import heph_functional

# List of fields needed 
Fields_needed = []

def initfields():
    #---------------------------------------------------------------------------
    # Go over the needed densities and the functional terms and check whether
    # we have enough derivatives to calculate the fields. 
    for term in heph_functional.Functional_terms:
        (densities, cpl) = heph_functional.ParseDensities(term)
        # Count the number of derivatives needed in this term
        totalder = 0
        totallap = 0
        for den in densities:
            (der, lap, left, right, coupling) = ParseOperators(den)
            totalder = totalder + der
            totallap = totallap + lap
        
        # Now see that for all densities in this term, the minimum number
        # of derivatives is the total one 
        for i in range(len(densities)):
            den = densities[i]
            (der, lap, left, right, coupling) = ParseOperators(den)
            for j in range(len(Densities_needed)):
                altden = Densities_needed[j]
                (altder, altlap, altleft, altright, altcoupling) = ParseOperators(altden)    
                if(altleft == left and altright == right):
                    # Set minimum derivatives
                    Der_den_needed[j] = max(totalder, Der_den_needed[j])
                    Lap_den_needed[j] = max(totallap,  Lap_den_needed[j])
        
def GenerateFields():
    #---------------------------------------------------------------------------
    # Generate a list of fields based on list of terms in the functional. 
    #
    #---------------------------------------------------------------------------
        
    global sumindices,tab

    field_decl_temp = Template(  tab + 'real(KIND=dp), allocatable :: $FIELD(:$DECLIND,:) \n')
    field_allo_temp = Template(2*tab + 'if(.not.allocated($FIELD)) then\n'+\
                           3*tab + 'allocate($FIELD(mv$ALLOCIND,2)) \n'   +\
                           2*tab + 'endif \n' + \
                           2*tab + '$FIELD = 0.0 \n')
    field_calc_temp    = Template( 2*tab + '$FIELD(:$IND,it) = $FIELD(:$IND,it)  & \n')
    
    field_calc_den_a     = Template('* $DENSITY(:$DENIND,it)  ')
    field_calc_den_b     = Template('* $DENSITY(:$DENIND,3-it)')
    field_calc_b_temp  = Template( 3*tab + '& $SIGN sum($CPLCTE(:,2)) $EXPR1 & \n') 
    field_calc_c_temp  = Template( 3*tab + '& $SIGN      $CPLCTE(2,2) $EXPR2 & \n') 

    doloop_template    = 2*tab + 'do %s = 1, 3 \n'
    enddoloop_template = 2*tab + 'enddo \n'

    cplcts    = []       
    for term in heph_functional.Functional_terms:      
        cplcts.append(term.replace('E_', 'B_'))
    #---------------------------------------------------------------------------
    # For every unique density encountered, we need to figure out the field
    # and the action of the field. 
    FIELDCALC   = ''
    declaration = ''
        
    for den in heph_functional.Densities_needed:
        # Name the field correctly
        dic = {}
        dic['FIELD'] = den.replace('D', 'F').replace('C', 'G')

        Fields_needed.append(dic['FIELD'])

        #  Get the operator structure of the density correctly                                
        (der, lap, left, right, coupling) = ParseOperators(den) 
        
        # Check all of the terms if they depend on the density
        fieldlist = []
        cpcte     = ''
        for term in heph_functional.Functional_terms: 
            (densities, cpl) = heph_functional.ParseDensities(term)
            
            # Change the coupling if the density needed is contracted
            if(OrderOfDen(den) != OrderOfDen(den, contract=False)):
                altterm = term
                for c in cpl: 
                    for i in range(len(densities)):
                        if(densities[i].count(sumindices[cpl.index(c)]) == 2):
                            # Replace internal couplings
                            altterm = altterm.replace(sumindices[cpl.index(c)],'')
                (rubbish, cpl) = heph_functional.ParseDensities(altterm)
            
            for i in range(len(densities)):
                altden = densities[i]
                (altder, altlap, altleft, altright, altcoup)=ParseOperators(altden)
                if(altleft == left and altright == right):
                    #  Add the term to the fieldlist for this density, 
                    #  and additionnally mentioning the number of external 
                    #  derivatives and laplacians
                    removed = []
                    for j in range(len(densities)):
                        if i != j :
                            removed.append(densities[j])
                    cplct = cplcts[heph_functional.Functional_terms.index(term)]
                    fieldlist.append([removed, altder, altlap, cplct, cpl])
        
        # Replace the densities in the list by the ones actually calculated
        for i in range(len(fieldlist)):
            d = fieldlist[i][0]
            for j in range(len(fieldlist[i][0])):
                (der,lap,left,right,coup) = ParseOperators(d[j])
                for altden in heph_functional.Densities_needed:
                    (altder, altlap, altleft, altright, altcoup) = ParseOperators(altden)
                    if(altleft == left and right == altright):
                        fieldlist[i][0][j] =   der*'der_' +               \
                                               lap*'Lap_' +               \
                                               altden
        # Create the expression for the field
        dic['ALLOCIND']= ''
        dic['DECLIND'] = ''
        
        for k in range(OrderOfDen(den)):
            dic['ALLOCIND'] = dic['ALLOCIND'] + ',3' 
            dic['DECLIND']  = dic['DECLIND']  + ',:'
       
        declaration  = declaration + field_decl_temp.substitute(dic)
        FIELDCALC    = FIELDCALC + field_allo_temp.substitute(dic)

    
        for fieldterm in fieldlist:
             order = 0
             for d in fieldterm[0]:
                order = order + OrderOfDen(d)
 
             for k in range(order):
                FIELDCALC = FIELDCALC + doloop_template%sumindices[k]
             
             # get the indices of the field correct
             dic['IND']     = ''
             for k in range(OrderOfDen(den)):
                for c in fieldterm[4]:
                    if k in c:
                        dic['IND']      = dic['IND']  + ',%s'%sumindices[fieldterm[4].index(c)]
       
             FIELDCALC = FIELDCALC + field_calc_temp.substitute(dic)
             
             dic['DENSITY']  = ''
             dic['EXPR1']    = ''
             dic['EXPR2']    = '' 
             for i in range(len(fieldterm[0])):
                    
                    dic['DENSITY']  =                    fieldterm[1] * 'der_' \
                                                       + fieldterm[2] * 'Lap_' \
                                                       + fieldterm[0][i] 
                    if( fieldterm[1]%2 == 0) :
                            dic['SIGN']     =  '+'
                    else:
                            dic['SIGN']     =  '-'
                    dic['CPLCTE']   =  fieldterm[3]
                    
                    # Get the indices of the density in the field
                    dic['DENIND']      = ''
                    for k in range(OrderOfDen(den), OrderOfDen(den) + OrderOfDen(dic['DENSITY'])):
                        for c in fieldterm[4]:
                            if k in c:   
                                dic['DENIND']= dic['DENIND']     + ',' \
                                             + sumindices[fieldterm[4].index(c)]
                    
                    dic['EXPR1'] = dic['EXPR1'] + field_calc_den_a.substitute(dic)
                    dic['EXPR2'] = dic['EXPR2'] + field_calc_den_b.substitute(dic)
             
             FIELDCALC = FIELDCALC + field_calc_b_temp.substitute(dic)
             FIELDCALC = FIELDCALC + field_calc_c_temp.substitute(dic)
             FIELDCALC = FIELDCALC[:-4] + '\n'
             for k in range(order):
                FIELDCALC = FIELDCALC + enddoloop_template

    return(declaration,FIELDCALC)

def GenerateAction(field):
    #---------------------------------------------------------------------------
    #
    #
    #
    #---------------------------------------------------------------------------
    
    action_final    = Template(  tab + \
    'hpsi(:,:,:,$IND) =  hpsi(:,:,:,$IND) $SIGN $TEMP(:,:,:$RIND,$RCOMP)\n')
    temp_ini        = tab + 'temp = 0.0 \n'
    action_temp    = Template(2*tab + \
    'temp(i,1,1,$IND) =  temp(i,1,1,$IND) $SIGN $FIELD(i$sFIELDIND,it) * $WF(i,1,1$RIND,$RCOMP)\n')
    
    derive_temp     = Template(tab + \
                      'call Derive_$DIR(temp(:,:,:,$RCOMP), $SYM, dtemp(:,:,:,$DIRIND,$RCOMP)) \n')
#    derive_2_temp   = Template(tab + \
#              'call Derive_grad(dtemp(:,:,:,$DERIND,$RCOMP), +1, +1, +1, dtemp(:,:,:,1,$DERIND,$RCOMP), dtemp(:,:,:,2,$DERIND,$RCOMP), dtemp(:,:,:,3,$DERIND,$RCOMP)) \n')
    lap_temp        = Template(tab + \
                      'call Derive_lap(temp(:,:,:,$RCOMP), +1, +1, +1, laptemp(:,:,:,$RCOMP)) \n')
    
    
    
    position_loop   = tab + 'do i=1,mv\n'
    position_end    = tab + 'enddo    \n'
    
    action_comment  = Template(tab + '!' + 75*'-' + '\n'\
                          +    tab + '! Action of $FIELD \n')
    
    WFNames = ['psi', 'dpsi', 'ddpsi', 'dddpsi', 'ddddpsi']
    Direction = ["X", 'Y', 'Z']
    
    (left,right,coupling) =  ParseOperatorsField(field)
    
    #---------------------------------------------------------------------------
    #Building the left and right operators
    operatordic = {}
    operatordic['I'] = Identity
    operatordic['N'] = Nabla
    operatordic['S'] = Sigma
    operatordic['C'] = Current
    
    LeftOperator = Identity
    for i in range(len(left)):
        # this needs to be done in reverse order
        l = left[len(left) - i -1 ] 
        LeftOperator  = Combine(operatordic[l], LeftOperator)
        
    RightOperator = Identity
    for i in range(len(right)):
        r = right[(len(right)) -i -1] 
        RightOperator = Combine(operatordic[r], RightOperator)
    
    #Dimension of the field without contractions
    ndim = LeftOperator.dimension + RightOperator.dimension - 2*len(coupling)
    
    dic = {}
    dic['FIELD']  = field
    
    start  = np.zeros((4,1))
    start[0,0] = 1 
    start[1,0] = 2  
    start[2,0] = 3 
    start[3,0] = 4    
        
    # Construct an iterator with all possible combinations of uncontracted 
    # indices
    expression = action_comment.substitute(dic)
    
    # So this is quite complicated. 
    # Steps:
    #    1) Construct all of the possible arguments for the left-operator
    #       consistent with possible contractions within the left-operator
    #    2) Construct all of the possible argument for the right-operator
    #       consstent with possible contractions within the right-operator
    #       AND consistent with contractions with the left-operator.
    # Repeat-----
    #    3) Write the code that calculates the action of the field for 
    #       all right arguments consistent with the left argument given
    #    4) Write the code that derives all of these right arguments for this
    #       particular left argument. 
    # End 

    ldim = LeftOperator.dimension
    
        
    lcoupl = []
    rcoupl = []
    ccoupl = []
    for c in coupling: 
       if(c[0] < LeftOperator.dimension and c[1] < LeftOperator.dimension ):
            ldim = ldim - 1
            lcoupl.append(c)
       elif(c[0] >= LeftOperator.dimension and c[1] >= LeftOperator.dimension ):
            rcoupl.append(c)
       else:
            ccoupl.append(c)
    rdim  = RightOperator.dimension - len(coupling) + len(lcoupl)

    # all possible values for the arguments of the left-operator
    largs = itertools.product(range(3), repeat=ldim)
    
    for larg in largs:
        # Construct all of the possibilites for the larg, ignoring 
        # contractions
        larg_uncontracted = []
        if(len(lcoupl) == 0):
            larg_uncontracted = [larg]
        else:
            cont = itertools.product(range(3), repeat=len(lcoupl))
            for c in cont:
                p  = ()   
                ii = 0
                for i in range(LeftOperator.dimension):
                    found = False                    
                    for combination in lcoupl:
                        if(i in combination): 
                            p = p + (c[lcoupl.index(combination)],)
                            found = True
                    if(not found): 
                        p = p + (larg[ii],)
                        ii = ii +1
                larg_uncontracted.append(p)
        
        for true_larg in larg_uncontracted:
            rargs = itertools.product(range(3), repeat=rdim)
            for rarg in rargs:
                rarg_uncontracted = []
                if(len(coupling) - len(rcoupl)== 0):
                  rarg_unconstracted=[rarg]
                else:  
                  rarg_uncontracted = []
                  cont = itertools.product(range(3), repeat=len(rcoupl))
                  for c in cont:
                    p  = ()   
                    ii = 0
                    for i in range(LeftOperator.dimension, LeftOperator.dimension+RightOperator.dimension):
                        found = False                    
                        for combination in rcoupl:
                            if(i in combination): 
                                p = p + (c[rcoupl.index(combination)],)
                                found = True
                        for combination in ccoupl:
                            if(i == combination[0]):
                                p = p + (true_larg[combination[1]],)
                                found = True
                            if(i == combination[1]):
                                p = p + (true_larg[combination[0]],)
                                found = True
                        if(not found): 
                            p = p + (rarg[ii],) #rarg[ii]
                            ii = ii +1
                    rarg_uncontracted.append(p)
                
                for true_rarg in rarg_uncontracted:
                    print field, true_larg, true_rarg
#    args = itertools.product(range(3), repeat=ndim)   
#    for arg in args:
#        expression = expression + temp_ini
#        uncontracted = []
#        if(len(coupling) == 0):
#            uncontracted = [arg]
#        else:
#            cont = itertools.product(range(3), repeat=len(coupling))
#            for c in cont:
#                p  = ()   
#                ii = 0
#                for i in range(LeftOperator.dimension + RightOperator.dimension):
#                        found = False                    
#                        for combination in coupling:
#                                if(i in combination): 
#                                        p = p + (c[coupling.index(combination)],)
#                                        found = True
#                        if(not found): 
#                                p = p + (arg[ii],)
#                                ii = ii +1
#                uncontracted.append(p)

#        for true_arg in uncontracted:
#            larg = true_arg[:LeftOperator.dimension]
#            rarg = true_arg[LeftOperator.dimension:]
#            
#            leftind  = LeftOperator(larg, start)
#            rightind = RightOperator(rarg, start)

#            #-------------------------------------------------------------------
#            # For now, we derive everything on the fly.
#            #-------------------------------------------------------------------
#            dic['WF']     = WFNames[ RightOperator.derorder ]

#            lorder = LeftOperator.derorder
#            
#            RIND = ''
#            for mu in rarg: 
#                RIND = RIND + ',' + str(mu+1) # Python indexes 0:N-1
#            dic['RIND'] = RIND
#            
#            
#            FIELDIND = ''
#            print field, rarg
#            for i in range(len(rarg)):
#                r = rarg[i]
#                Found = False
#                for c in coupling:
#                    if r in c:
#                        Found = True
#                if(not Found):
#                    FIELDIND = FIELDIND + ',' + str(r+1)               
#            dic['FIELDIND'] = FIELDIND
#            
#            expression = expression + position_loop
#            
#            #-------------------------------------------------------------------
#            # Action of the right-operator
#            for k in range(4): 
#                dic['IND']    = k + 1  
#                dic['RCOMP']  = int(abs(rightind[k,0])) 
#                
#                SIGN          = np.sign(rightind[k,0])
#                SIGN          = SIGN * (-1)**(LeftOperator.derorder)
#                if(SIGN > 0) :
#                    dic['SIGN']   = '+'
#                else :
#                    dic['SIGN']   = '-'
#                expression = expression + action_temp.substitute(dic)
#            expression = expression + position_end
#            #-------------------------------------------------------------------
#            # Calculation of the derivatives
#            if(lorder >= 1):
#                dic['DIRIND'] = larg[-1] +1 
#                dic['DIR']    = Direction[larg[-1]]
#                dic['SYM']    = '+1'       
#                for k in range(4):
#                    dic['RCOMP']  = k  +1 
#                    expression = expression + derive_temp.substitute(dic)
#            for k in range(4):
#                dic['IND']    = k + 1  
#                dic['RCOMP']  = k + 1 
#                dic['TEMP']   = lorder*'d'+'temp'  
#                expression = expression + action_final.substitute(dic)
    return (expression)
    
def ParseOperatorsField(field):    
    #---------------------------------------------------------------------------
    # Parse the operators that are used to construct the action of field from
    # its name only.
    #---------------------------------------------------------------------------
    left  = ''
    right = '' 

    split = field.split('_')

    left = split[1]
    right= split[2]
    #---------------------------------------------------------------------------
    # Don't put the couplings in the definition of left- and right-operators
    for l in sumindices:
        left  =  left.replace(l, '')
        right = right.replace(l, '')
        
    if(field[0] == 'G'):      
            right = 'C' + right
    #---------------------------------------------------------------------------
    # Find the coupling
    coupling  = []
    for l in sumindices:
        c   = ()
        ind = 0
        for i in range(len(field)):      
            if(field[i] == l):
                c = c+ (ind,)
            if(field[i] in sumindices):
                ind = ind + 1 
        if(len(c) > 0) :
            coupling.append(c)
    return(left, right, coupling)
