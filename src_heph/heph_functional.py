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
# procedure. It should be built as follows
#   *) One term per line, separations by point-commas (;)
#   *) First entry is the name of the term, following the naming scheme
#   *) Second entry is a possible density dependence of the term
#   *) third entry is the isospin-scalar coupling constant, defined in 
#      function of quantities known by Tantalus.
#   *) fourth entry is the isospin-vector coupling constant, similary defined.
#   *) Lines starting with '!' will be ignored.
#   *) Lines starting with '#' will not have consequences for the code,
#      but will be printed to the output of Hephaestos.
# Example, the rho^2 term
#    E_D_I_I_D_I_I ;   ; +3.0/8.0*t0+0*t0*x0+0*t3+0*t3*x3 ; 
#                                            -1.0/8.0*t0-0.25*t0*x0+0*t3+0*t3*x3 
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
from heph_densities import Der_den_needed, Lap_den_needed
#-------------------------------------------------------------------------------
# Array containing the expressions of all the functional terms. 
Functional_terms     = []
density_dependence   = []
#-------------------------------------------------------------------------------
# Strings telling Tantalus how to calculate coupling constants from Skyrme
# force values for the isoscalar (_0) and isovector (_1) coupling. 
coupling_constants_0 = []
coupling_constants_1 = []

def initfunctional(fname):

    global Functional_terms, Densities_needed
    
    # Read the functional from a given file
    description = ReadFunctional(fname)
    
    # Generating the list of all densities
    tempden     = []
    rotationals = []
    for term in Functional_terms:
        (densities,coup) = ParseDensities(term)
        for den in densities:
            tempden.append(den)
                    
    # Pruning the list
    # A) removing duplicates
    # B) removing contractions when the full density will be calculates
    Densities_needed.append(tempden[0])
    Der_den_needed.append(0)
    Lap_den_needed.append(0)
    for i in range(len(tempden)):
        (deri, lapi, lefti, righti, coupi) = ParseOperators(tempden[i])
        Found = False
        for j in range(len(Densities_needed)):
            (derj, lapj, leftj, rightj, coupj) = ParseOperators(Densities_needed[j])
            #-------------------------------------------------------------------
            # Two densities are identical if the left- and right-operatos
            # are the same.
            if(leftj == lefti and rightj == righti):
                # Signal that the density is already present
                Found = True
                # However, check that we don't need any new derivatives
                # If so, add them
                Der_den_needed[j] = max(Der_den_needed[j], deri)
                Lap_den_needed[j] = max(Lap_den_needed[j], lapi)
                
                # then check if the coupling of the indices is the same
                # Note that the loop starts over coupj, since that one has by
                # definition more couplings than coupi in it. Thus, this 
                # logic will not fail if coupi has no elements.
                for cj in coupj:
                    Found_coup = False
                    for ci in coupi:
                        if(cj[0] == ci[0] and cj[1] == ci[1]):
                            Found_coup = True
                        elif(cj[1] == ci[1] and cj[0] == ci[0]):
                            Found_coup = True
                    if(not Found_coup):
                        l = sumindices[coupj.index(cj)]
                        Densities_needed[j] = Densities_needed[j].replace(l, '')
            #-------------------------------------------------------------------
        if(not Found):
            Densities_needed.append(tempden[i].replace(derstring + '_','')\
                                              .replace(lapstring + '_', ''))
            Der_den_needed.append(deri)
            Lap_den_needed.append(lapi) 
    
        
    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -' 
    print ' Functional taken from file %s'%fname
    print ' Description from file:'
    print  description.replace('#', tab)
    print ' Number of terms:     %d'%len(Functional_terms)
    print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'
    
    
def ReadFunctional(fname):
    #
    # Read the functional terms and the coupling coefficients from file fname.
    #
    global Functional_terms

    description=''
    with open(fname, 'r') as f:
        for line in f: 
            if(line[0] != '#' and line[0] != '!'):
                split = line.split(';')
                Functional_terms.append(split[0].replace(' ', ''))
                density_dependence.append(split[1].replace(' ', ''))
                coupling_constants_0.append(split[2].replace(' ', ''))
                coupling_constants_1.append(split[3].replace(' ', ''))   
            elif(line[0] == '#'):
                description = description + line     

    return(description)
    
def ParseDensities(term): 
    #---------------------------------------------------------------------------
    # From a term in the functional represented by a string, parse all of the
    # densities that make it up.
    #---------------------------------------------------------------------------
    densities = []
    # Split along C and D-s
    temp      = ''
    # Don't take into account any Density dependence naming
    split     = term.replace('_DD', '').split('_') 
    for i in range(len(split)):
        if split[i][0:3] == derstring or split[i] == lapstring:
                temp  = temp + split[i] + '_'               
                
        if split[i] == 'D' or split[i] == 'C':
                temp = temp + split[i] + '_' + split[i+1] + '_' + split[i+2]
                densities.append(temp)
                temp = ''

    # Don't propagate couplings that are not internal to the density                
    for i in range(len(densities)):
        for l in sumindices: 
            if(densities[i].count(l) != 2):
                densities[i] = densities[i].replace(l, '')
        
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
    return (densities, coupling)

def ProcessFunctional(fname, src, target):
        declaration = ''
        calculation = ''
        form        = ''
        printing    = ''
        calccoef    = ''
        printcoef   = ''
        sumtotal    = ''
        fieldcalc   = ''

        for i in range(len(Functional_terms)): 
                (d,c,p,cc, pc,st) = GenTermExpression(Functional_terms[i], \
                            [coupling_constants_0[i], coupling_constants_1[i]],\
                                                   density_dependence[i])

                declaration = declaration + d + '\n'
                calculation = calculation + c + '\n'
                printing    = printing    + p + '\n'
                calccoef    = calccoef    + cc+ '\n'
                printcoef   = printcoef   + pc+ '\n'
                sumtotal    = sumtotal    + st+ '&\n&'

#        (fielddec, fieldcalc) = heph_fields.GenerateFields(     )
#        declaration = declaration + fielddec + '\n'
        # - - - - - - - - - - - - - - - - - - - - -
        # Substitute into the functional.f90 file.        
        dic={}
        dic['DECLARATION']    = declaration
        dic['CALCULATION']    = calculation
        dic['PRINT']          = printing
        dic['CALCCOEF']       = calccoef   
        dic['PRINTCOEF']      = printcoef 
        dic['TOTAL']          = sumtotal[:-3]
        dic['CALCFIELDS']     = fieldcalc
        with open(src+fname, 'r') as template:
                with open(target+fname, 'w') as generated:
                    for line in template:
                        generated.write(Template(line).substitute(dic))  

def GenTermExpression( term, ccoef, DD):

    global  sumindices

    declaration = ''
    calculation = ''
    printing    = ''
    
    #---------------------------------------------------------------------------
    # Templates for the declaration, calculation and printing of an energy term.
    # 
    decl_template   = Template('real(KIND=dp) :: $CPCTE(2,2), $TERM(2,2)')
    edent_template  = Template('sum($DEN(:$IND,:),2)')
    edenq_template  = Template('$DEN(:$IND,it)')
    
    comment_template= Template(   tab + '!' + 38 * '- ' + '\n' +               \
                                  tab + '! Calculation of $TERM \n')
                                  
    doloop_template    =    tab + 'do %s = 1, 3 \n'
    enddoloop_template =    tab + 'enddo \n'
   
    sumtotal_template  = Template( 2*tab + ' + $TERM(:,1)')
    
    calc_z_template = Template(   tab + 'Edensity = 0.0_dp \n')
    calc_a_template = Template(   tab + 'EDensity(:,3) = Edensity(:,3) + $EDENT\n')
    calc_b_template = Template(   tab + '$TERM(2,2) = 0.0_dp \n' +                      \
                                  tab + 'do it=1,2 \n' +                                \
                                2*tab + 'Edensity(:,it) = Edensity(:,it) + $EDENQ \n' + \
                                  tab + 'enddo \n')
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
                                 
    print_template      = Template(" print('(a30 , 3f15.6)'), '$TERM', $TERM(:,1), sum($TERM(:,1))")
    print_cpl_template  = Template(" print('(a30 , 4f15.6)'), '$CPCTE', $CPCTE")
    #---------------------------------------------------------------------------
    # See how many indices are present everywhere.
    (tempden, coupling) = ParseDensities(term)
    doloops             = len(coupling)
    #Now see how these densities are present in the heph_densities.py module
    densities = []
    for den in tempden:
        (der,lap,left,right, coupl) = ParseOperators(den)
        for i in range(len(Densities_needed)):
            (derref, lapref, leftref, rightref, couplref) = \
                                             ParseOperators(Densities_needed[i])
            if(left == leftref and right == rightref):
                addden = lap*'Lap_' + der*'der_' + Densities_needed[i]
                densities.append( addden )
                # Don't do do loops over indices that already had been contracted
                doloops = doloops + \
                        OrderOfDen(addden) - OrderOfDen(addden, contract=False)
                

    orders                = []
    for i in range(len(densities)): 
        orders.append(OrderOfDen(densities[i])) 
    
    dic = {}
    index_encountered=0
    name = ''
    for i in range(len(densities)):
            name = name + '_' + densities[i]

    dic ['TERM' ] = term
    dic ['CPCTE'] = 'B' + dic ['TERM'][1:]    
  
    dic['EDENT'] = ''
    dic['EDENQ'] = ''
    prevorder = 0 
    for i in range(len(densities)):
        isodic = {}
        isodic['DEN'] = densities[i]
        
        isodic['IND'] = ''
        for l in range(prevorder, prevorder + orders[i]):
            for c in coupling:
                if( l in c ):
                    isodic['IND'] = isodic['IND'] + ',' + sumindices[coupling.index(c)]
        
        dic['EDENT'] = dic['EDENT'] + edent_template.substitute(isodic) + '*'
        dic['EDENQ'] = dic['EDENQ'] + edenq_template.substitute(isodic) + '*'
            
        # Take out the final '*' which should not be necessary
        prevorder = prevorder + orders[i]
    dic['EDENT'] = dic['EDENT'][:-1]
    dic['EDENQ'] = dic['EDENQ'][:-1]  
    
    declaration = decl_template.substitute(dic) 

    calculation = comment_template.substitute(dic)
    calculation = calculation = calculation + calc_z_template.substitute(dic)
    for i in range(doloops):
        calculation  = calculation + doloop_template%sumindices[i]
        
    calculation = calculation + calc_a_template.substitute(dic)
    calculation = calculation + calc_b_template.substitute(dic)
    for i in range(doloops):
        calculation = calculation + enddoloop_template
    if(DD != ''):
        dic['DD'] = DD
        calculation = calculation + calc_DD_template.substitute(dic)
    calculation = calculation + calc_c_template.substitute(dic) 
    calculation = calculation + calc_d_template.substitute(dic) + '\n'
    calculation = calculation + calc_e_template.substitute(dic)
    calculation = calculation + calc_f_template.substitute(dic)
    
    printing = print_template.substitute(dic) 
    
    dic['EXP1'] = ccoef[0]
    dic['EXP2'] = ccoef[1]
    
    calccoef = calc_coef_template.substitute(dic)
    printcoef = print_cpl_template.substitute(dic)    
    sumtotal  = sumtotal_template.substitute(dic)
    
    return (declaration, calculation, printing, calccoef, printcoef, sumtotal)
