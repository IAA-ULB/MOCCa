#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
#
# Module governing the treatment of the functional in Hephaestos, 
# for writing to Tantalus source files. 
#
# The type of functional is read from file, passed into the initfunctional
# procedure. 
# Example, the density dependent term from standard Skyrme
# E_D_I_I_D_I_I_DD ;  
#                 (sum(D_I_I,2)**(yt3a)) ;  D_I_I ; 
#          yt3a*(sum(D_I_I,2)**(yt3a-1)) ;   yt3a ; 
#                    3.0_dp/48.0_dp * t3 ; - 1.0_dp/24.0_dp * t3 * (0.5_dp + x3) 
#
# Once read the code will parse the input for
#   a) the densities needed to calculate these terms
#   b) the coupling between different indices in the expression of the term
#       Example: E_D_Nm_Nm_D_Nn_Nn vs E_D_Nm_Nn_D_Nm_Nn
#                tau^2             vs tau_mn tau_mn
# 
# Afterwards, the code determines which densities actually need calculating. 
# This is a non-trivial, task. For example: D_N_N is not actually needed in its
# full form with an NLO functional, as the contraction D_Nm_Nm suffices. 
#-------------------------------------------------------------------------------

from string        import Template
import itertools
import numpy as np
from heph_densities import Densities_needed, tab, sumindices, derstring
from heph_densities import lapstring, OrderOfDen, ParseOperators
from heph_densities import crossindices, Storage_Mapping, Multiplicity
from heph_densities import deriv_needed, Pair_densities_needed

import heph_linechecker
import heph_fields
#-------------------------------------------------------------------------------
# Array containing the expressions of all the functional terms. 
Functional_terms      = []
Functional_pair_terms = []
density_dependence    = []
field_DD_terms        = {}
DD_rearcoefs          = []
#-------------------------------------------------------------------------------
# Strings telling Tantalus how to calculate coupling constants from Skyrme
# force values for the isoscalar (_0) and isovector (_1) coupling. 
coupling_constants_0 = []
coupling_constants_1 = []

coupling_constants_pair_0 = []
coupling_constants_pair_1 = []

#-------------------------------------------------------------------------------
# Strings telling Tantalus which parameters are used to compute coupling 
# constants
paramparameters = []

#-------------------------------------------------------------------------------
# Switch determining what order of derivatives is needed to be computed.
#
# Note that this is not directly i => ith order.
#  derivative_order <=> derivatives of spwfs calculated
#            1          1st order + trace of 2nd order (laplacian)
#            2          2nd order + trace of 3rd order
#            3          3rd order derivatives
derivative_order = 1
#-------------------------------------------------------------------------------
# Assume whether or not the functional is local. With this == 1, the 
# script will use the following simplification
#
#  Action-of-C^{1,N} =  -i F^{1, N}  \nabla 
#
#  instead of the full
#
#  Action-of-C^{1,N} =  -i \frac{1}{2} [ F^{1, N} \nabla + \nabla F^{1, N} ]
#
# Similar relations for C densities with odd number of derivatives are used too.
#
# Note that I still fail to account for this formally (except for C^{1,N}), but 
# this seems to hold if the functional is local.
#-------------------------------------------------------------------------------
assume_locality = 1

def initfunctional(fname, fpairname):

    global Functional_terms, Densities_needed, derivative_order
    
    # Read the functional from a given file
    (description, pair_description) = ReadFunctional(fname, fpairname)
    
    # Generating the list of all densities
    tempden     = []
    rotationals = []
    for term in Functional_terms:
        (densities,coup) = ParseDensities(term)
        for den in densities:
            tempden.append(den)
    #---------------------------------------------------------------------------
    # Pruning the list
    # A) removing duplicates
    # B) removing contractions when the full density will be calculated
    Densities_needed.append(tempden[0])
    deriv_needed.append([])
    for i in range(len(tempden)):
        (deri, lapi, lefti, righti, coupi, crossi) = ParseOperators(tempden[i])
        Found = False
        for j in range(len(Densities_needed)):
            (derj, lapj, leftj, rightj, coupj, crossj) = ParseOperators(Densities_needed[j])
            #-------------------------------------------------------------------
            # Two densities are identical if the left- and right-operators
            # are the same.
            if(leftj == lefti and rightj == righti):
                # Signal that the density is already present
                Found = True
                # However, if there is a vector coupling in one, that is not in 
                # other, just calculate both. 
                if(crossj != crossi):
                    Found = False
                # However, check that we don't need any new derivatives
                # If so, add them
                deriv_needed[j].append((lapi, deri))
                # then check if the coupling of the indices is the same
                # Note that the loop starts over coupj, since that one has by
                # definition more couplings than coupi in it. Thus, this 
                # logic will not fail if coupi has no elements.
                for cj in coupj + crossj:
                    Found_coup = False
                    for ci in coupi:
                        if(cj == ci):
                            Found_coup = True
                    for ci in crossi:
                        if(cj == ci ):
                            Found_coup = True   
                    if(not Found_coup):
                        try:                        
                            l = sumindices[coupj.index(cj)]
                        except ValueError:  
                            l = sumindices[crossj.index(cj)]
                        Densities_needed[j] = Densities_needed[j].replace(l, '')
            #-------------------------------------------------------------------
        if(not Found):
            add = tempden[i] 
            # Getting the duplicates out of the derivatives
            deriv_needed.append([(lapi,deri)])
            
            # Remove all of the derivatives from the top
            add = add.replace(derstring + '_','').replace(lapstring+'_', '')
            for l in sumindices:
                add = add.replace(derstring + l + '_','')
            Densities_needed.append(add)
    #---------------------------------------------------------------------------     
    # Finding out how many derivatives we need to take of the spwfs
    derivative_order = 1    
    for den in Densities_needed:
        (der, lap, left, right, cpl, cross) = ParseOperators(den)
        ders = right.count('N')
        derivative_order = max(derivative_order, ders)

    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -' 
    print ' P-h functional taken from file %s'%fname
    print ' Description from file:'
    print  description.replace('#', tab)
    print ' Number of terms:      %d'%len(Functional_terms)
    print ' Order of derivatives: %d'%derivative_order
    print ' Locality assumed:     %d'%assume_locality
    print ' # Parameters          %d'%len(paramparameters)
    #print   paramparameters
    for i in range(len(paramparameters)/3):
        print '  ', paramparameters[3*i:3*i+3]    
    if(len(paramparameters)%3 != 0):
        print '  ', paramparameters[3*(i+1):3*(i+1)+len(paramparameters)%3]
    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -' 
    print ' P-p functional taken from file %s'%fpairname
    print ' Description from file:'
    print  pair_description.replace('#', tab)
    print ' Number of terms:      %d'%len(Functional_pair_terms)
    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

def PruneDeriv_needed():
    #---------------------------------------------------------------------------
    # Add all of the possible combinations with less derivatives and laplacians,
    # so that we can build the eventually needed combinations.
    for i in range(len(deriv_needed)):
        newderiv=[]
        for j in deriv_needed[i]:
            newderiv.append(j)
            
            options = itertools.product(range(j[0]+1), range(j[1]+1))
            for opt in options:
                newderiv.append(opt)
                
        deriv_needed[i] = newderiv
    # b) remove duplicates in the list of needed derivative combinations
    # c) Order the list in increasing level of operations
    for i in range(len(deriv_needed)):
        deriv_needed[i] = list(set(deriv_needed[i]))
        deriv_needed[i] = sorted(deriv_needed[i])
        
def ReadFunctional(fname, fpairname):
    #-------------------------------------------------------------------------
    # Read the functional terms and the coupling coefficients from the
    # ph functional file (fname) and the pp functional file (fpairname).
    #-------------------------------------------------------------------------
    global Functional_terms, Functional_pair_terms

    description      = ''
    pair_description = ''
    #-------------------------------------------------------------------------
    # Read the ph functional
    termsstart = 0
    with open(fname, 'r') as f:
        for line in f:
            try:
                if(len(line.split()) == 0):
                    print 'a'
                    continue
                elif(line[0] == '#'):
                     #-sign indicates description                   
                     description = description + line
                     continue
                elif(line[0:6] == '!TERMS'):
                    # Signal that the constants part is over
                    termsstart  = 1   
                    print paramparameters 
                    continue            
                elif(line[0] == '!'):
                    continue
                
                # Parse constants only
                if(termsstart == 0):
                    split = line.split(';')
                    for s in split:
                        paramparameters.append(s.replace(' ','').replace('\n', ''))
                elif(termsstart == 1):
                    #  Start parsing terms in the functional   
                    split = line.split(';')
                    Functional_terms.append(split[0].replace(' ', ''))
                    coupling_constants_0.append(split[1].replace(' ', ''))
                    coupling_constants_1.append(split[2].replace(' ', ''))   
                        
                    if(len(split)>3):
                        density_dependence.append(split[3].replace(' ', ''))
                        dd_den   = split[4].replace(' ', '')
                        f_dd     = split[5].replace(' ', '')
                        field_DD_terms[split[0].replace(' ', '')] =(dd_den,f_dd)
                        DD_rearcoefs.append(split[6].replace(' ', ''))
                    else:
                        density_dependence.append('')
                        field_DD_terms[split[0].replace(' ', '')] =  ('','')
                        DD_rearcoefs.append('')
            except IndexError:
                print 'Problem reading the following line in the func file.'
                print line
                exit()      
    #---------------------------------------------------------------------------
    # Read the pp functional
#    with open(fpairname, 'r') as f:
#        for line in f: 
#            try:
#                if(len(line.split()) == 0):
#                    continue
#                if(line[0] != '#' and line[0] != '!'):
#                    split = line.split(';')
#                    Functional_pair_terms.append(split[0].replace(' ', ''))                
#                    
#                    coupling_constants_0.append(split[5].replace(' ', ''))
#                    coupling_constants_1.append(split[6].replace(' ', ''))   
#                                
#                elif(line[0] == '#'):
#                    pair_description = pair_description + line
#            except IndexError:
#                print 'Problem reading the following line in the pairing-functional file.'
#                print line
#                exit()      
    return(description, pair_description)
    
def ParseDensities(term): 
    #---------------------------------------------------------------------------
    # From a term in the functional represented by a string, parse all of the
    # densities that make it up.
    #---------------------------------------------------------------------------
    densities = []
    #---------------------------------------------------------------------------
    # Split along C and D-s
    temp      = ''
    #---------------------------------------------------------------------------
    # Don't take into account any density dependence naming
    split     = term.replace('_DD', '').split('_') 
    for i in range(len(split)):
        if split[i][0:3] == derstring or split[i] == lapstring:
                temp  = temp + split[i] + '_'               
                
        if split[i] == 'D' or split[i] == 'C':
                temp = temp + split[i] + '_' + split[i+1] + '_' + split[i+2]
                densities.append(temp)
                temp = ''
    #---------------------------------------------------------------------------
    # Find the coupling
    coupling  = []
    foundx    = []
    for l in sumindices:
        c   = ()
        ind = 0
        foundx.append(0)
        for i in range(len(term)):      
            if(term[i] == l):
               if(term[i-1] != 'x'): 
                    c = c+ (ind,)
               else:    
                 if(foundx[sumindices.index(l)] == 0):
                    foundx[sumindices.index(l)] = foundx[sumindices.index(l)] +1
                    c = c+ (ind,)
                        
            if(term[i] in sumindices):
                    ind = ind + 1 
        if(len(c) > 0) :
                coupling.append(c)
    #---------------------------------------------------------------------------
    # Don't propagate couplings that are not between left and right operators
    for l in sumindices:
        for i in range(len(densities)): 
            if( derstring + l in densities[i]) :
                # Remove the coupling if it involves derivatives
                densities[i] = densities[i].replace(l, '')
                for j in range(len(densities)):
                    densities[j] = densities[j].replace(l, '')
            for j in range(len(densities)):
                if i == j:
                    pass
                elif( l in densities[i] and l in densities[j]):
                    # Remove the coupling if it is between more densities
                    densities[i] = densities[i].replace(l, '')
                    densities[j] = densities[j].replace(l, '')
                else:
                    pass
        
    return (densities, coupling)

def ProcessParameterization(fname, src, target):
        #-----------------------------------------------------------------------
        # Processing of the parameterization.f90 file to include the different 
        # parameters.
        decl_template = Template( tab + 'real(KIND=dp) :: $PARAM = -12345 \n')

        decl = ''

        for s in paramparameters:
            dic= {}
            dic['PARAM'] = s
            
            decl = decl + decl_template.substitute(dic)

        dic= {}
        dic['PARAMDECL'] = decl
        with open(src+fname, 'r') as template:
            with open(target+fname, 'w') as generated:
                for line in template:
                    generated.write(Template(line).substitute(dic))  

def ProcessFunctional(fname, src, target):
        #-----------------------------------------------------------------------
        # Master routine calling the other ones to generate a functional based
        # on the parsing done before.
        #-----------------------------------------------------------------------
    
        declaration = ''
        calculation = ''
        form        = ''
        printing    = ''
        calccoef    = ''
        printcoef   = ''
        sumtotal    = ''
        fieldcalc   = ''
        erear       = ''
        
        #-----------------------------------------------------------------------
        # Generate the terms in the functional
        for i in range(len(Functional_terms)): 
                (d,c,p,cc, pc,st,er)  = GenTermExpression(Functional_terms[i], \
                            [coupling_constants_0[i], coupling_constants_1[i]],\
                                         density_dependence[i], DD_rearcoefs[i])

                declaration = declaration + d + '\n'
                calculation = calculation + c + '\n'
                printing    = printing    + p + '\n'
                calccoef    = calccoef    + cc+ '\n'
                printcoef   = printcoef   + pc+ '\n'
                sumtotal    = sumtotal    + st+ '&\n'
                erear       = erear       + er

        sumtotal = sumtotal[:-2]

        #-----------------------------------------------------------------------
        # Generate the fields of the single-particle hamiltonian
        (fielddec, fieldcalc) = heph_fields.GenerateFields(     )
        declaration = declaration + fielddec + '\n'
        
        #-----------------------------------------------------------------------
        # Generate the expressions for the actions of the fields
        SkyrmeAction = ''
        for field in heph_fields.Fields_needed:
            #-------------------------------------------------------------------
            # Check if we need to symmetrize the action
            #-------------------------------------------------------------------
            (left,right,coupling,cross) = heph_fields.ParseOperatorsField(field)
            #-------------------------------------------------------------------
            # Generate the expression for the application of the ordinary 
            # operator structure
            if( left != right):
                if( 'C' in left or 'C' in right):   
                    # Only symmetrize non-symmetric C's if asked for
                    if(assume_locality == 1):
                        SkyrmeAction = SkyrmeAction +                          \
                                            heph_fields.GenerateAction(field, 0)
                    else:
                        SkyrmeAction = SkyrmeAction +                          \
                                            heph_fields.GenerateAction(field, 1)
                        SkyrmeAction = SkyrmeAction +                          \
                                            heph_fields.GenerateAction(field,-1)
                else:
                    # Always symmetrize non-symmetric D's
                    SkyrmeAction = SkyrmeAction + heph_fields.GenerateAction(field,+1)
                    SkyrmeAction = SkyrmeAction + heph_fields.GenerateAction(field,-1)
            else:
                SkyrmeAction = SkyrmeAction + heph_fields.GenerateAction(field, 0)
        #-----------------------------------------------------------------------
        # Now make sure all of the lines are not too long for compilation.
        declaration = heph_linechecker.LineFormat(declaration)
        calculation = heph_linechecker.LineFormat(calculation)
        printing    = heph_linechecker.LineFormat(printing)
        calccoef    = heph_linechecker.LineFormat(calccoef)
        printcoef   = heph_linechecker.LineFormat(printcoef)
        sumtotal    = heph_linechecker.LineFormat(sumtotal)
        fieldcalc   = heph_linechecker.LineFormat(fieldcalc)
        SkyrmeAction= heph_linechecker.LineFormat(SkyrmeAction)
        erear       = heph_linechecker.LineFormat(erear)
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
        # Substitute into the functional.f90 file.        
        dic={}
        dic['DECLARATION']    = declaration
        dic['CALCULATION']    = calculation
        dic['PRINT']          = printing
        dic['CALCCOEF']       = calccoef   
        dic['PRINTCOEF']      = printcoef 
        dic['TOTAL']          = sumtotal
        dic['CALCFIELDS']     = fieldcalc
        dic['SKYRMEACTION']   = SkyrmeAction
        dic['EREAR']          = erear
        
        if(derivative_order == 1):
          dic['N2'] = ' '    
          dic['N3'] = '!'
        elif(derivative_order == 2): 
          dic['N2'] = ' '
          dic['N3'] = '!'
        elif(derivative_order == 3):
          dic['N2'] = '!'
          dic['N3'] = ' '
        
        with open(src+fname, 'r') as template:
                with open(target+fname, 'w') as generated:
                    for line in template:
                        generated.write(Template(line).substitute(dic))  

def GenTermExpression( term, ccoef, DD, DDrear):
    #---------------------------------------------------------------------------
    # Generate the expressions for the terms in the functional.
    #
    #---------------------------------------------------------------------------
    global  sumindices

    declaration = ''
    calculation = ''
    printing    = ''
    
    #---------------------------------------------------------------------------
    # Templates for the declaration, calculation and printing of an energy term.
    # And, not forgetting, its contribution to the rearrangement energy
    decl_template   = Template( tab + 'real(KIND=dp) :: $CPCTE(2,2), $TERM(2,2)')
    edent_template  = Template('sum($DEN(:$IND,:),2)')
    edenq_template  = Template('$DEN(:$IND,$IT)')
    
    comment_template= Template(   tab + '!' + 38 * '- ' + '\n' +               \
                                  tab + '! Calculation of $TERM \n')
                                  
    end_comment     =             tab + '!' + 38 * '- ' + '\n'
                                  
    doloop_template    =    tab + 'do %s = 1, 3 \n'
    enddoloop_template =    tab + 'enddo \n'
   
    sumtotal_template  = Template( 2*tab + ' & + $TERM(:,1)')
    
    calc_z_template = Template(   tab + 'Edensity = 0.0_dp \n')
    calc_a_template = Template(   tab + 'EDensity(:,3) = Edensity(:,3) + $EDENT\n')
    calc_b_template = Template(   tab + 'Edensity(:,1) = Edensity(:,1) + $EDENN\n' + \
                                  tab + 'Edensity(:,2) = Edensity(:,2) + $EDENP\n'  )
    calc_DD_template= Template(   tab + 'do m=1,3 \n' +                                 \
                                2*tab + 'EDensity(:,m) = Edensity(:,m) * $DD \n' +      \
                                  tab + 'enddo \n')
    calc_c_template = Template(   tab + '$TERM(1,2) = $CPCTE(1,2) * sum( Edensity(:,3)) * dv \n')
    calc_d_template = Template(   tab + '$TERM(2,2) = $CPCTE(2,2) * dv & \n'+           \
                                  tab + '&'+6*tab+' * sum(Edensity(:,1) + Edensity(:,2) ) \n')
    calc_e_template = Template(   tab + '$TERM(1,1) = $CPCTE(1,1)/$CPCTE(1,2) * $TERM(1,2)\n')
    calc_f_template = Template(   tab + '$TERM(2,1) = sum($TERM(:,2)) - $TERM(1,1) \n')
    calc_coef_template = Template(tab + '$CPCTE(1,1) = $EXP1 \n' + \
                                  tab + '$CPCTE(2,1) = $EXP2 \n' + \
                                  tab + '$CPCTE(1,2) = $CPCTE(1,1) - $CPCTE(2,1) \n' + \
                                  tab + '$CPCTE(2,2) =             2*$CPCTE(2,1) \n')       
                                 
                                 
    write_edensity = Template( tab + ' call output_Edensity(Edensity, "$FILENAME")' )
    print_template      = Template(tab +" print('(a30 , 3f15.6)'), '$TERM', $TERM(:,1), sum($TERM(:,1)) \n")
    print_cpl_template  = Template(tab +" print('(a30 , 4f15.6)'), '$CPCTE', $CPCTE")
    
    rear_template       = Template(tab +" e_rear = e_rear $REARCOEF*sum($TERM(:,2))\n")
    
    #---------------------------------------------------------------------------
    # See how many indices are present everywhere.
    (tempden, coupling) = ParseDensities(term)
    doloops             = 0 #len(coupling)
  
    #Now see how these densities are present in the heph_densities.py module
    densities = []
    for den in tempden:
        (der,lap,left,right, coupl, cross) = ParseOperators(den)
        for i in range(len(Densities_needed)):
            (derref, lapref, leftref, rightref, couplref, crossref) = \
                                             ParseOperators(Densities_needed[i])
            if(left == leftref and right == rightref and cross == crossref):
                addden = lap*'Lap_' + der*'Der_' + Densities_needed[i]
                densities.append( addden )
                # Don't do do loops over indices that already had been contracted
                doloops = doloops + OrderOfDen(addden) 
                # Go back to the outer loop
                break
    
    orders                = []
    for i in range(len(densities)): 
        orders.append(OrderOfDen(densities[i])) 

    # Clever trick to recount the couplings of the term
    # => Find the couplings that are inside a given density
    # => Remove them from the term
    # => recount the couplings (not the densities!)
    altterm = term
    for c in coupling: 
        for i in range(len(densities)):
            if(densities[i].count(sumindices[coupling.index(c)]) == 2):
                # Replace internal couplings
                altterm = altterm.replace(sumindices[coupling.index(c)],'')
                
    (rubbish, true_coupling) = ParseDensities(altterm)
    
    doloops = doloops - len(true_coupling)
    
    dic = {}
    index_encountered=0
    name = ''
    for i in range(len(densities)):
            name = name + '_' + densities[i]

    dic ['TERM' ] = term
    dic ['CPCTE'] = 'B' + dic ['TERM'][1:]    
  
    #---------------------------------------------------------------------------
    # arguments for all the couplings
    args = list(itertools.product(range(3), repeat=len(true_coupling)))
            
    declaration = decl_template.substitute(dic) 
    calculation = comment_template.substitute(dic)
    calculation = calculation = calculation + calc_z_template.substitute(dic)
    for arg in args: 
        dic['EDENT'] = ''
        dic['EDENP'] = ''
        dic['EDENN'] = ''
        prevorder = 0 
        for i in range(len(densities)):
            isodic = {}
            isodic['DEN'] = densities[i]
            
            #-------------------------------------------------------------------
            # Get the index of the density correct
            (der,lap,left,right, coupl, cross) = ParseOperators(densities[i])
            
            isodic['IND'] = ''
            indices = ()
            for l in range(prevorder, prevorder + orders[i]):
                for c in true_coupling:
                    if( l in c ):
                       indices = indices + (arg[true_coupling.index(c)],) 
                        
            # The first indices are necessarily external derivatives
            if(der > 0):
                derind = Storage_Mapping(indices[:der])
                indices = (derind,) + indices[der:]
                
            for l in indices:
                isodic['IND'] = isodic['IND'] + ',' + str(l+1)
            
            dic['EDENT'] = dic['EDENT'] + edent_template.substitute(isodic) + '*'
            
            isodic['IT'] = 1
            dic['EDENN'] = dic['EDENN'] + edenq_template.substitute(isodic) + '*'
            isodic['IT'] = 2
            dic['EDENP'] = dic['EDENP'] + edenq_template.substitute(isodic) + '*'    
            # Take out the final '*' which should not be necessary
            prevorder = prevorder + orders[i]
      
        dic['EDENT'] = dic['EDENT'][:-1]
        dic['EDENN'] = dic['EDENN'][:-1]  
        dic['EDENP'] = dic['EDENP'][:-1]  
        
        calculation = calculation +'\n' + tab + '! indices = ' + str(arg) + '\n'
        calculation = calculation + calc_a_template.substitute(dic)
        calculation = calculation + calc_b_template.substitute(dic)

    calculation = calculation + '\n'
    if(DD != ''):
        dic['DD'] = DD
        calculation = calculation + calc_DD_template.substitute(dic)

    dic['FILENAME'] ='edensities/' + dic['TERM'] + '.dat'
    #calculation = calculation + write_edensity.substitute(dic) + '\n'
    calculation = calculation + calc_c_template.substitute(dic) 
    calculation = calculation + calc_d_template.substitute(dic) + '\n'
    calculation = calculation + calc_e_template.substitute(dic)
    calculation = calculation + calc_f_template.substitute(dic)
    
    calculation = calculation + end_comment + '\n'
    
    printing = print_template.substitute(dic) 
    
    dic['EXP1'] = ccoef[0]
    dic['EXP2'] = ccoef[1]
    
    calccoef  = calc_coef_template.substitute(dic)
    printcoef = print_cpl_template.substitute(dic)    
    sumtotal  = sumtotal_template.substitute(dic)
    
    # Getting the contribution to the rearrangement energy
    # Two-body, non-density dependent terms don't have rearrangement terms.
    if( len(densities) == 2 and DD == ''):
        erear = ''
    else:
        rearcoef = '+' + '(' + str(len(densities) - 2)
        if(DDrear != '') :
            rearcoef = rearcoef  + '+' + DDrear + ')'
        else:
            rearcoef = rearcoef + ')'
        dic['REARCOEF'] = rearcoef
        erear = rear_template.substitute(dic)
        
    return (declaration, calculation, printing, calccoef, printcoef, sumtotal, erear)
