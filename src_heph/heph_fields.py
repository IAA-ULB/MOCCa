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
#
#   C Add option for more than one derivative and/or laplacian to the action-
#     of fields routine
#   C Add possibility for non-derivative operators in the left-operator, as 
#     currently all indices are assumed to be derivatives ones.
#   C add automatic declaration of dtemp, ddtemp, lapdtemp etc to the actions
#     to save memory in default cases.
#-------------------------------------------------------------------------------
import itertools
import numpy as np
import src_heph.heph_functional
from src_heph.heph_densities    import *
from string import Template

# List of fields needed 
Fields_needed = []
Pairing_Fields_needed = []

def initfields():
    #---------------------------------------------------------------------------
    # Go over the needed densities and the functional terms and check whether
    # we have enough derivatives to calculate the fields. 
    for term in src_heph.heph_functional.Functional_terms:
        (densities, cpl) = src_heph.heph_functional.ParseDensities(term)
        # Count the number of derivatives needed in this term
        totalder = 0
        totallap = 0
        for den in densities:
            (der, lap, left, right, coupling, cross) = ParseOperators(den)
            totalder = totalder + der
            totallap = totallap + lap
        
        # Now see that for all densities in this term, the minimum number
        # of derivatives is the total one 
        for i in range(len(densities)):
            den = densities[i]
            (der, lap, left, right, coupling, cross) = ParseOperators(den)
            for j in range(len(Densities_needed)):
                altden = Densities_needed[j]
                (altder, altlap, altleft, altright, altcoupling, altcross)     \
                                                        = ParseOperators(altden)    
                if(altleft == left and altright == right):
                    # Set minimum derivatives
                    deriv_needed[j].append((totallap, totalder))
                    
    
    src_heph.heph_functional.PruneDeriv_needed()
    
      
def GenerateFields():
    #---------------------------------------------------------------------------
    # Generate a list of fields based on list of terms in the functional. 
    #
    #---------------------------------------------------------------------------
        
    global sumindices,tab

    field_decl_temp = Template(  tab + 'real(KIND=dp), allocatable, target :: $FIELD(:$DECLIND,:) \n')
    fhist_decl_temp = Template(  tab + 'real(KIND=dp), allocatable :: ${FIELD}_hist(:$DECLIND,:) \n')
    
    field_allo_temp = Template(2*tab + 'if(.not.allocated($FIELD)) then\n'+\
                           3*tab + 'allocate($FIELD(mv$ALLOCIND,2)) \n'   +\
                           3*tab + 'allocate(${FIELD}_hist(mv$ALLOCIND,2)) \n'+\
                           3*tab + '$FIELD = 0.0 ; ${FIELD}_hist = 0.0 \n' + \
                           2*tab + 'endif \n')

    field_allo_b_temp = Template(3*tab + 'if(.not.allocated($FIELD)) then\n'+\
                           4*tab + 'allocate($FIELD(filemv$ALLOCIND,2)) \n'   +\
                           4*tab + 'allocate(${FIELD}_hist(mv$ALLOCIND,2)) \n'+\
                           4*tab + '$FIELD = 0.0 ; ${FIELD}_hist = 0.0 \n' + \
                           3*tab + 'endif \n')

    field_hist_temp = Template( 2*tab + 'if(calcall) then \n' + 
                                3*tab + '${FIELD}_hist = $FIELD \n' + 
                                3*tab + '$FIELD = 0.0 \n'           + 
                                2*tab + 'endif \n')

    field_line             = Template(2*tab+ \
    '!----------------------------------------------------------------------\n')
    field_calc_temp_a    = Template( 2*tab + '! Calculation of $FIELD \n')
    field_calc_temp_b    = Template( 3*tab + '$FIELD(:$IND,it) = $FIELD(:$IND,it)  & \n')

    field_calc_den_a     = Template('* sum($DENSITY(:$DENIND,:),$SUMIND)  ')
    field_calc_den_b     = Template('* $DENSITY(:$DENIND,it)')
    field_calc_den_c     = Template('* $DENSITY(:$DENIND,3-it)')

    # Note the convention for coupling constants for the pairing term is 
    # different. We do not deal with isoscalar and isovector coupling constants
    # but rather with the coupling constants for neutrons and protons.
    field_calc_den_pair  = Template('* $DENSITY(:$DENIND,it)')

    isoloop = Template(2*tab + 'maxit = 2 \n' + \
              2*tab + 'if((.not.calcall).and.(any(${FIELD}.ne.0.0))) maxit=0\n'+\
              2*tab + 'do it=1,maxit \n')
    
    isoloop_end = 2*tab + 'enddo\n'

    field_calc_b_temp  = Template( 3*tab + '& $SIGN $DD $CPLCTE(1,2)  $EXPR1 & \n') 
    field_calc_c_temp  = Template( 3*tab + '& $SIGN $DD $CPLCTE(2,2)  $EXPR2 & \n') 
    field_calc_d_temp  = Template( 3*tab + '& $SIGN $DD $CPLCTE(2,2)  $EXPR3 & \n') 

    field_pair_a_temp    = Template( 3*tab + '& $SIGN $DD $CPLCTE(it,1)   $EXPR2  & \n') 
    field_pair_b_temp    = Template( 3*tab + '& $SIGN $DD $CPLCTE(3-it,1) $EXPR3  & \n') 

    doloop_template    = 2*tab + 'do %s = 1, 3 \n'
    enddoloop_template = 2*tab + 'enddo \n'

    field_write_template_a = Template(tab + ('write(chan, iostat=io) "$FIELDFILLED" \n'))
    field_write_template_b = Template(tab + ('write(chan, iostat=io)  $FIELD  \n'))

    field_read_template_a  = Template(2*tab + ('case("$FIELD") \n'))
    field_read_template_b  = Template(2*tab + ( tab + 'read(chan, iostat=io)  $FIELD \n'))

    field_transfo_temp = Template(  3*tab + 'do it=1,2 \n'  \
                                  + 4*tab + '${FIELD}_hist(:$IND,it) = & \n'  
                                  + 4*tab + '&  changeboxsize_function($FIELD(:$IND,it), filenx, fileny, filenz) \n'\
                                  + 3*tab + 'enddo \n ')
    field_transfo_end  = Template(3*tab + '$FIELD = ${FIELD}_hist \n '\
                                 +3*tab + '${FIELD}_hist = 0.0d0 \n')

    # Template for cleaning fields
    clean_template   = Template(   tab+'if(allocated($FIELD)) deallocate($FIELD)')
    clean_template_b = Template(   tab+'if(allocated(${FIELD}_hist)) deallocate(${FIELD}_hist)')

    cplcts    = []       
    for term in src_heph.heph_functional.Functional_terms:      
        cplcts.append(term.replace('E_', 'B_'))
    #---------------------------------------------------------------------------
    # For every unique density encountered, we need to figure out the field
    # and the action of the field. 
    FIELDCALC   = ''
    declaration = ''

    fieldread = ''
    fieldwrite= ''
    fieldtransfo = ''
        
    fieldclean= '' 
    for den in src_heph.heph_functional.Densities_needed:
        #-----------------------------------------------------------------------
        # Name the field correctly
        dic = {}
        dic['FIELD'] = den.replace('D', 'F').replace('C', 'G')

        #Get a string of length 30 with the field name, but with extra spaces
        # at the end
        dic['FIELDFILLED'] = dic['FIELD'].ljust(30)

        if('P' not in den): 
          Fields_needed.append(dic['FIELD'])
        else:
          Pairing_Fields_needed.append(dic['FIELD'])
        #-----------------------------------------------------------------------
        #  Get the operator structure of the density correctly                                
        (der, lap, left, right, coupling,cross) = ParseOperators(den) 
        #-----------------------------------------------------------------------
        # Check all of the terms if they depend on the density
        fieldlist = []
        cpcte     = ''

        fieldclean = fieldclean + '\n' + clean_template.substitute(dic)
        fieldclean = fieldclean + '\n' + clean_template_b.substitute(dic)

        #-----------------------------------------------------------------------
        for term in src_heph.heph_functional.Functional_terms: 
            (densities, cpl) =src_heph. heph_functional.ParseDensities(term)
            #-------------------------------------------------------------------
            # Replace the densities in the list by the ones actually calculated
            for i in range(len(densities)):
                (x,y,l2,r2,c2,cr2) = ParseOperators(densities[i])
                for altden in src_heph.heph_functional.Densities_needed:
                    (altder, altlap, altleft, altright, altcoup, altcross) = ParseOperators(altden)
                    if(altleft == l2 and r2 == altright and cr2 == altcross):
                        densities[i] = y*'Lap_' +               \
                                       x*'Der_' +               \
                                       altden
            #-------------------------------------------------------------------
            # Remove all the mentions of couplings inside the density if only
            # contractions are calculated.
            # 
            altterm = term
            for c in cpl: 
                for i in range(len(densities)):
                    if(OrderOfDen(densities[i]) != OrderOfDen(densities[i], contract=False)):
                      for s in sumindices:
                          if(densities[i].count(s) == 2):
                            # Replace internal couplings
                            altterm = altterm.replace(s,'')
                        
#                     if(densities[i].count(sumindices[cpl.index(c)]) == 2):
#                          # Replace internal couplings
#                          altterm = altterm.replace(sumindices[cpl.index(c)],'')

            (rubbish, cpl) = src_heph.heph_functional.ParseDensities(altterm)
            #-------------------------------------------------------------------
            # Check if the term contains this density
            startind = 0
            for i in range(len(densities)):
                altden = densities[i]
                (altder, altlap, altleft, altright, altcoup, altcross)=ParseOperators(altden)
                if((altleft == left) and (altright == right) and (cross == altcross)):
                    #  Add the term to the fieldlist for this density, 
                    #  and additionnally mentioning the number of external 
                    #  derivatives and laplacians
                    removed    = []
                    removedsum = 0
                    for j in range(len(densities)):
                        if i != j :
                            removed.append(densities[j])
                    for j in range(i):
                            removedsum = removedsum + OrderOfDen(densities[j])
                    # Reparse the term with this density at the front, so that
                    # couplings are correctly referenced
                    newcpl = []
                    shift = OrderOfDen(den) + altder
                    for c in cpl:
                        nc = ()
                        for k in c:
                            if( k < startind): 
                                nc = nc + (k+shift,)
                            elif(k > startind + shift + 1 ):
                                nc = nc + (k,)
                            else:
                                nc = nc + (k-removedsum,)
                        newcpl.append(nc)        

                    cpl = newcpl
                    
                    ind   = src_heph.heph_functional.Functional_terms.index(term)
                    cplct = cplcts[ind]
                    dden  = src_heph.heph_functional.density_dependence[ind]
                    fieldlist.append([removed, altder, altlap, cplct, cpl,dden,0])
                
                startind = startind + OrderOfDen(altden)
                
            #-------------------------------------------------------------------
            # Now check if there are density dependences in this term that 
            # involve this density
            dd    = src_heph.heph_functional.field_DD_terms[term]
            ind   = src_heph.heph_functional.Functional_terms.index(term)
            cplct = cplcts[ind]
            if(dd[0] == den):
                fieldlist.append([densities, altder, altlap, cplct, cpl, dd[1],1])

        # Create the expression for the field
        dic['ALLOCIND']= ''
        dic['DECLIND'] = ''
        for k in range(OrderOfDen(den)):
            dic['ALLOCIND'] = dic['ALLOCIND'] + ',3' 
            dic['DECLIND']  = dic['DECLIND']  + ',:'
       
        declaration  = declaration + field_decl_temp.substitute(dic)
        declaration  = declaration + fhist_decl_temp.substitute(dic)

        fieldread    = fieldread   + field_read_template_a.substitute(dic)
        fieldread    = fieldread   + field_allo_b_temp.substitute(dic)
        fieldread    = fieldread   + field_read_template_b.substitute(dic)
        
        fieldwrite   = fieldwrite  + field_write_template_a.substitute(dic)
        fieldwrite   = fieldwrite  + field_write_template_b.substitute(dic)
        
        FIELDCALC    = FIELDCALC + field_line.substitute(dic)
        FIELDCALC    = FIELDCALC + field_calc_temp_a.substitute(dic)
        FIELDCALC    = FIELDCALC + field_allo_temp.substitute(dic)
        FIELDCALC    = FIELDCALC + field_hist_temp.substitute(dic)
        FIELDCALC    = FIELDCALC + isoloop.substitute(dic)

        args = list(itertools.product(range(3), repeat=OrderOfDen(den)))
        for arg in args:   
           # get the indices of the field correct
           dic['IND']     = ''
           for k in arg:
                  dic['IND'] = dic['IND'] + ',%s'%(k+1)

           fieldread = fieldread + field_transfo_temp.substitute(dic)
             
        for fieldterm in fieldlist:
             # Find the number of indices over which there have to be sums
             NumberOfIndices = 0
             for d2 in [den] + fieldterm[0]:
                #(der, lap, left, right, cpl, crs) = ParseOperators(d2)
                NumberOfIndices= NumberOfIndices + OrderOfDen(d2) 

             # The couplings however do not need individual indices 
             if(NumberOfIndices != 0):
               NumberOfIndices = NumberOfIndices - len(fieldterm[4]) 

             # But the external derivatives do            
             NumberOfIndices = NumberOfIndices + fieldterm[1]

             #------------------------------------------------------------------
             # arguments for all the indices
             args = list(itertools.product(range(3), repeat=NumberOfIndices))
             
             for arg in args:
                 # get the indices of the field correct
                 dic['IND']     = ''
                 for k in range(OrderOfDen(den)):
                    for c in fieldterm[4]:
                        if k in c:
                            dic['IND'] = dic['IND'] \
                                          + ',%s'%(arg[fieldterm[4].index(c)]+1)
           
                 FIELDCALC = FIELDCALC + field_calc_temp_b.substitute(dic)
                 
                 dic['DENSITY']  = ''
                 dic['EXPR1']    = ''
                 dic['EXPR2']    = ''
                 dic['EXPR3'] = ''
                 ind = src_heph.heph_functional.Functional_terms.index(term)
                 dic['DD']       = fieldterm[5]
                 
                 if(len(dic['DD']) >0 ):
                    dic['DD']       = dic['DD'] + '*'
                 
                 lastorder = OrderOfDen(den)
                 for i in range(len(fieldterm[0])):
                    dic['DENSITY']  =                    fieldterm[2] * 'Lap_' \
                                                       + fieldterm[1] * 'Der_' \
                                                       + fieldterm[0][i] 
                    #-----------------------------------------------------------
                    # Make sure the combination of laplacians and derivatives
                    # is in the right ordering 
                    dercount = dic['DENSITY'].count('Der')
                    lapcount = dic['DENSITY'].count('Lap')
                    
                    dic['DENSITY'] = dic['DENSITY'].replace('Der_', '').replace('Lap_', '')
                    dic['DENSITY'] = lapcount * 'Lap_' + dercount * 'Der_' + dic['DENSITY']
                    
                    if( fieldterm[1]%2 == 0) :
                            dic['SIGN']     =  '+'
                    else:
                            dic['SIGN']     =  '-'
                    dic['CPLCTE']   =  fieldterm[3]
                    
                    #-----------------------------------------------------------
                    # Get the indices of the density in the field
                    indices = ()
                    for k in range(lastorder, lastorder + OrderOfDen(dic['DENSITY'])):
                        for c in fieldterm[4]:
                            if k in c:   
                                indices = indices + (arg[fieldterm[4].index(c)],)
                    #-----------------------------------------------------------
                    # The first indices are necessarily external derivatives
                    if(dercount > 0):
                        derind  = Storage_Mapping(indices[:dercount])
                        indices = (derind,) + indices[dercount:]
                    
                    dic['DENIND']      = ''
                    dic['SUMIND']      = 2 
                    for l in indices:
                        dic['DENIND'] = dic['DENIND'] + ',' + str(l+1)

                    dic['EXPR1'] = dic['EXPR1'] + field_calc_den_a.substitute(dic)
                    dic['EXPR2'] = dic['EXPR2'] + field_calc_den_b.substitute(dic)
                   
                    if(len(dic['DD']) >0):
                        dic['EXPR3'] = dic['EXPR3'] + field_calc_den_c.substitute(dic)
                        lastorder = lastorder + OrderOfDen(dic['DENSITY'])

                 if('P' not in dic['DENSITY']):
                    # Ordinary mean-field densities; coupling constants are
                    # the isoscalar and isovector ones
                    FIELDCALC = FIELDCALC + field_calc_b_temp.substitute(dic)
                    FIELDCALC = FIELDCALC + field_calc_c_temp.substitute(dic)
                    if(fieldterm[6] == 1):
                       FIELDCALC = FIELDCALC + field_calc_d_temp.substitute(dic)
                    FIELDCALC = FIELDCALC[:-4] + '\n \n'
                 else:
                    # Pairing densities, coupling constants are pn ones.  
                    FIELDCALC = FIELDCALC + field_pair_a_temp.substitute(dic)       
                    if(fieldterm[6] == 1):
                      FIELDCALC = FIELDCALC + field_pair_b_temp.substitute(dic)
                    FIELDCALC = FIELDCALC[:-4] + '\n \n'
                                
        fieldread = fieldread + field_transfo_end.substitute(dic)

        FIELDCALC    = FIELDCALC + isoloop_end
        FIELDCALC    = FIELDCALC + field_line.substitute(dic)

    return(declaration, FIELDCALC, fieldwrite, fieldread, fieldclean)

def GenerateAction(field, symmetrize):
    #---------------------------------------------------------------------------
    #
    #
    #
    #---------------------------------------------------------------------------
    
    #---------------------------------------------------------------------------
    WFNames   = ['psi', 'dpsi', 'ddpsi', 'dddpsi', 'ddddpsi']
    Direction = ["X", 'Y', 'Z']
    #---------------------------------------------------------------------------
    # Templates to fill in.
    #---------------------------------------------------------------------------
    action_final            = Template(  tab + \
    'hpsi(:,$IND) =  hpsi(:,$IND) $SIGN $LMULT $TEMP(:$LIND,$RCOMP)\n')
    action_final_pairing    = Template(  tab + \
    'deltapsi(:,$IND) =  deltapsi(:,$IND) $SIGN $LMULT $TEMP(:$LIND,$RCOMP)\n')
    
    
    temp_ini        = tab + 'temp = 0.0 \n'
    action_temp    = Template(2*tab + \
    'temp(i,$IND) =  temp(i,$IND) $SIGN $RMULT $FIELD(i$FIELDIND,it) * $WF(i$RIND,$RCOMP)\n')
    
    derive_temp     = Template(tab + \
    'call Derive_$DIR($DNUMBER(:,$RCOMP), $SYM, d$DNUMBER(:,$DIRIND,$RCOMP)) \n')
    lap_temp        = Template(tab + \
    'call Derive_lap(temp(:,$RCOMP), $SYMX, $SYMY, $SYMZ, laptemp(:,$RCOMP)) \n')
    
    position_loop   = tab + 'do i=1,mv\n'
    position_end    = tab + 'enddo    \n'
    
    action_comment  = Template(tab + '!' + 75*'-' + '\n'\
                          +    tab + '! Action of $FIELD symmetrized: $SYM \n')
    
    #---------------------------------------------------------------------------
    # Parse the field under consideration
    (left,right,coupling,cross) = ParseOperatorsField(field)
    #---------------------------------------------------------------------------
    # If symmetrize is non-zero, change left <-> right and the couplings
    # accordingly
    if(symmetrize == -1 ):
        switch_field = field.split('_')
        switch_field = switch_field[0]+'_'+switch_field[2]+'_'+ switch_field[1]
        (left,right,coupling,cross) = ParseOperatorsField(switch_field)
    #---------------------------------------------------------------------------
    #Building the left and right operators
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
        if (l in crossindices):
            continue
        LeftOperator  = Combine(operatordic[l], LeftOperator)
        
    RightOperator = Identity
    for i in range(len(right)):
        r = right[(len(right)) -i -1]
        if (r in crossindices):
            continue
        RightOperator = Combine(operatordic[r], RightOperator)
    
    #Dimension of the field without contractions
    ndim = LeftOperator.dimension + RightOperator.dimension                    
    ndim = ndim - 2*len(coupling) + len(cross)
    
    dic = {}
    dic['FIELD']  = field
    dic['WF']     = WFNames[ RightOperator.derorder ]

    start  = np.zeros((4,1))
    start[0,0] = 1 
    start[1,0] = 2  
    start[2,0] = 3 
    start[3,0] = 4    
    
    #---------------------------------------------------------------------------
    # Construct an iterator with all possible combinations of uncontracted 
    # indices
    dic['SYM'] = symmetrize
    expression = action_comment.substitute(dic)
    #---------------------------------------------------------------------------
    # So this is quite complicated. 
    # Steps:
    #    1) Construct all of the possible arguments for the left-operator,
    #       for the indices that are not contracted.
    #  |-----Loop over true_larg
    #  |
    #  |  2) Construct all of the possible arguments for the right-operator
    #  |     for the noncontracted indices.
    #  | |----- Loop over rarg
    #  | | 3) Construct all possible arguments for the right-operator, given
    #  | |    contractions with its own indices and with the left-operator
    #  | ||---- Loop over true_rarg
    #  | || 4)  Add the code to add the action of the right operator to temp
    #  | ||----
    #  | |----- 
    #  | 5) Add the code to derive the array temp, for the given indices of the
    #  |    left-operator. Contractions between derivatives here are exchanged
    #  |    with calls to laplacian.
    #  | 6) Add the final result to hpsi    
    #  |-----------------
    #---------------------------------------------------------------------------

    ldim = LeftOperator.dimension

    lcoupl = []
    rcoupl = []
    ccoupl = []
    fieldind_passed = []
    
    for c in coupling: 
       if(c[0] < ldim and c[1] < ldim ):
            lcoupl.append(c)
       elif(c[0] >= ldim and c[1] >= ldim ):
            rcoupl.append(c)
       else:
            ccoupl.append(c)

    rdim  = RightOperator.dimension - len(ccoupl) - 2*len(rcoupl)
    ldim  = LeftOperator.dimension  - 2*len(lcoupl)
    
    # all possible values for the arguments of the left-operator
    largs = list(itertools.product(range(3), repeat=ldim))
    #-------------------------------------------------------------------
    # Go over the uncontracted right-indices and get the independent  
    # components, and the multiplicities. 
    larg_stor = []
    larg_new  = []
    multiplicities         = []    
    for true_larg in largs: 
        l_stor = Storage_Mapping(true_larg[:LeftOperator.derorder])
        mult   = Multiplicity(true_larg[:LeftOperator.derorder])
        l_new  = (l_stor,) + true_larg[LeftOperator.derorder:]
    
        Found = False
        for new in larg_stor:
            if(new == l_new):
                Found = True
        if(not Found):
            larg_stor.append(l_new)
            larg_new.append(true_larg)
            multiplicities.append(mult)
    
    largs = larg_new
    
    for true_larg in largs:
        # reset temp to 0
        expression = expression + temp_ini
        
        # Get the multiplicity correct
        m = multiplicities[largs.index(true_larg)]
        
        if( m != 1):
            dic['LMULT'] = str(m) + ' * '
        else:
            dic['LMULT'] = ''
            
        if(symmetrize == 1 or symmetrize== -1):
            dic['LMULT'] = dic['LMULT'] + ' 0.5d0 * '
        
        leftind  = LeftOperator(true_larg, start)
        rargs    = list(itertools.product(range(3), repeat=rdim-len(cross)))
        dic['RMULT'] = ''
        
        master_rarg = []
        for rarg in rargs:

            rarg_uncontracted = []
            if(RightOperator.dimension == 0):  
                rarg_uncontracted=[rarg]
            else:  
              rarg_uncontracted = []
              cont = itertools.product(range(3), repeat=len(rcoupl))
              
              crossind = []
              for i in range(len(cross)):
                crossind = crossind  + (Rot_ind(rarg[len(coupling) + i]))   
        
              if(len(cross) == 0):
                full_cont = cont
              else:
                full_cont = []                
                for c in cont:
                    for x in crossind:
                        full_cont.append(c + x)

              for c in full_cont:
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
                    for combination in cross:
                        if (i==combination[0]):
                                p = p + (c[cross.index(combination) + len(coupling)],)
                                found = True
                        if(i==combination[1]): 
                                p = p + (c[cross.index(combination) + len(coupling) + 1 ],)
                                found = True
                    if(not found): 
                        p = p + (rarg[ii],)
                        ii = ii +1
                rarg_uncontracted.append(p)
            #-------------------------------------------------------------------
            # Loop over right-arguments
            expression = expression + position_loop
            for true_rarg in rarg_uncontracted:
                # Action of the right operator for this indices
                rightind = RightOperator(true_rarg, start)
                
                #---------------------------------------------------------------
                # Indices of the field in the multiplication
                dic['FIELDIND'] = ''
                for l in range(LeftOperator.dimension):
                    Found = False
                    for c in coupling:
                        if(l in c):
                            Found = True
                    if(not Found):
                        dic['FIELDIND']= dic['FIELDIND'] + ',' \
                                           + str(abs(true_larg[l])+1)
                for r in range(len(rarg)):
                    dic['FIELDIND']= dic['FIELDIND'] + ',' \
                                           + str(abs(rarg[r])+1)
                                           
                # Get the packed storage-scheme index
                rarg_stor = Storage_Mapping(true_rarg[:RightOperator.derorder])
                #---------------------------------------------------------------
                # Action of the right-operator
                for k in range(4): 
                    dic['IND']     = k + 1
                    
                    dic['RIND'] = ''
                    #-----------------------------------------------------------
                    # Attention, as this if-condition was the source of some
                    # confusion. 
                    if( len(true_rarg)>0  and RightOperator.derorder != 0):
                    #-----------------------------------------------------------
                        dic['RIND'] =  dic['RIND']  + ',' + str(rarg_stor +1 )                        
                    dic['RCOMP']   = int(abs(rightind[k,0])) 
                    
                    SIGN           = np.sign(rightind[k,0])
                    SIGN           = SIGN * (-1)**(LeftOperator.derorder)
                    for l in true_rarg:
                        if( l  == 0):
                           SIGN = SIGN
                        else:
                           SIGN = SIGN * np.sign(l)
                    if(SIGN > 0) :
                        dic['SIGN']= '+'
                    else :
                        dic['SIGN']= '-'
                    expression = expression + action_temp.substitute(dic)
            expression = expression + position_end
            #---------------------------------------------------------------
            # End of true_rarg loop
        #-------------------------------------------------------------------
        # End of rarg loop
        #-----------------------------------------------------------------------
        # Calculation of the derivatives in steps:
        lasttemp = 'temp'
        for lorder in range(LeftOperator.derorder):
            Found = False
            Add   = True
            for c in lcoupl:
                if (lorder == c[0]):
                    Found = True
                elif(lorder == c[1]):
                    # Already found
                    Add = False
            if(not Add):
                continue
            if(Found):
                # Note that a laplacian never changes the quantum numbers of 
                # a function, and that the result of the laplacian needs to be
                # added to hpsi, so as long as derivatives and laplacians can 
                # not be compounded, these are the symmetries of the spwf.
                for k in range(4):
                    dic['SYMX']    = 'sx(%d)'%(k+1)
                    dic['SYMY']    = 'sy(%d)'%(k+1)
                    dic['SYMZ']    = 'sz(%d)'%(k+1)      
                
                    dic['RCOMP']  = k  + 1 
                    expression = expression + lap_temp.substitute(dic)
                lasttemp = 'laptemp'
            else:
                offset = LeftOperator.dimension - LeftOperator.derorder
                direc  = true_larg[- lorder - offset - 1] +1  # X/Y/Z derivative
                
                dic['DIRIND'] = direc
                dic['DIR']    = Direction[direc-1]
                
                dic['DNUMBER']=  lorder* 'd' + 'temp'
                
                sym = +1
                # Get the symmetries of the derivatives (including the current one)
                # that still need to be performed after this derivative.
                # Note that contractions are not considered, since they represent
                # laplacians.
                for l in range(lorder, LeftOperator.derorder):
                    if true_larg[- l - offset - 1] +1  == direc:
                        sym = -sym
                
                for k in range(4):
                    dic['RCOMP'] = ''
                    for l in range(lorder):
                        dic['RCOMP']  = str(true_larg[-l-offset-1] +1) + ',' + dic['RCOMP'] 
                    dic['RCOMP']  = dic['RCOMP'] + str(k+1) 
                    if(sym>0):   
                        dic['SYM']    = '+s' + Direction[direc-1] + '(%d)'%(k+1)
                    else:
                        dic['SYM']    = '-s' + Direction[direc-1] + '(%d)'%(k+1)   
                    
                    expression = expression + derive_temp.substitute(dic)
                lasttemp = (lorder+1) * 'd' + 'temp'
        #-----------------------------------------------------------------------
        # Add final result to hpsi
        dic['TEMP'] = lasttemp
        for k in range(4):
            dic['LIND']    = ''
            for l in true_larg[0:LeftOperator.derorder]:
                dic['LIND'] = dic['LIND'] + ',' + str(l+1)
            dic['IND']     = k + 1
            dic['RCOMP']   = int(abs(leftind[k,0])) 
            SIGN           = np.sign(leftind[k,0])
            if((symmetrize == -1) and ('C' in left or 'C' in right)):
                SIGN = - SIGN

            if('P' in field):
              SIGN = -SIGN
            if(SIGN > 0) :
                dic['SIGN']= '+'
            else :
                dic['SIGN']= '-'
            
            if('P' not in field):
              expression = expression + action_final.substitute(dic)
            else:
              expression = expression + action_final_pairing.substitute(dic)
        expression = expression + '\n'
        #-----------------------------------------------------------------------
        # End of true_larg loop
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
    
    #
    # There doesn't need to be an explicit time-reversal operator in the 
    # definition of the act of Delta
    # 
    
    #---------------------------------------------------------------------------
    # Find the coupling
    coupling  = []
    foundsums = []
    for l in sumindices:
        c   = ()
        ind = 0
        for i in range(1,len(field)):      
            if(field[i] == l):
                c = c+ (ind,)
            if(field[i-1] in ['N', 'S']):
                ind = ind + 1 
        if(len(c) > 0) :
            coupling.append(c)
        if(len(c) == 3):            
            foundsums.append(l)
    #---------------------------------------------------------------------------
    # Find the coupling over the crossindices, but only if not contracted with
    # another index
    cross  = []
    for l in crossindices:
        c   = ()
        ind = 0
        for i in range(len(field)):      
            if(field[i] == l):
                if(i == len(field) - 1):
                    c = c+ (ind,)
                elif (field[i+1] not in foundsums):
                    c = c+ (ind,)
            if(field[i] in crossindices+sumindices):
                ind = ind + 1 
        if(len(c) > 0) :
            cross.append(c)
    
    return(left, right, coupling, cross)
