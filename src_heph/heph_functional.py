#============================================================================================
#          _______  _______           _______  _______  _______ _________ _______  _______ 
#|\     /|(  ____ \(  ____ )|\     /|(  ___  )(  ____ \(  ____ \\__   __/(  ___  )(  ____ \
#| )   ( || (    \/| (    )|| )   ( || (   ) || (    \/| (    \/   ) (   | (   ) || (    \/
#| (___) || (__    | (____)|| (___) || (___) || (__    | (_____    | |   | |   | || (_____ 
#|  ___  ||  __)   |  _____)|  ___  ||  ___  ||  __)   (_____  )   | |   | |   | |(_____  )
#| (   ) || (      | (      | (   ) || (   ) || (            ) |   | |   | |   | |      ) |
#| )   ( || (____/\| )      | )   ( || )   ( || (____/\/\____) |   | |   | (___) |/\____) |
#|/     \|(_______/|/       |/     \||/     \|(_______/\_______)   )_(   (_______)\_______)
#============================================================================================
from string import Template
import itertools
import numpy as np
################################################################################
#
#
#
#
#
# TODO
#    C Add options for density dependent terms
#    C add option  for derivative of density
#    C add option  for rotational 
################################################################################

def CPL_string(cplct, coefs_0, coefs_1):
  coefstrings = []
  for i in range(len(coefs_0[:,0])):
    coef0 = ''
    coef1 = ''
    for j in range(len(coefs_0[0,:])):
        coef0 = coef0 + str(coefs_0[i,j]) + '*' + cplct[j] + '+'
        coef1 = coef1 + str(coefs_1[i,j]) + '*' + cplct[j] + '+'
    coef0 = coef0[:-1]
    coef1 = coef1[:-1]
    coefstrings.append((coef0,coef1))
  return(coefstrings)


#-------------------------------------------------------------------------------
# Tab-character for the fortran code.
# 4 spaces for W.R., but I can imagine other people have different standards.
tab           = '    '
#-------------------------------------------------------------------------------
# Check how many format statements for printing the energy contributions have 
# been made. Currently there are already two format statements in the
# functional.f90 subroutine printEnergy
formatcounter = 3

#-------------------------------------------------------------------------------
# Indices over which sums are supposed to go in both the FORTRAN code and the 
# naming scheme.
sumindices    = ['m', 'n', 'k', 'l']

#-------------------------------------------------------------------------------
#
LO_skyrme  =              ['t0', 't0*x0']
LO_coefs_0 = np.array(   [[' 3.0/8.0',  0]                                       # rho^2
                        , [        0 ,  0]]  )                  
LO_coefs_1 = np.array(   [['-1.0/8.0', '-0.25']                                  # rho^2
                        , [       0  ,  0     ]] )              
#
#                                     
NLO_skyrme  =            ['t1',           't1*x1',           't2',    't2*x2']
NLO_coefs_0 = np.array(  [[' 3.0/16.0',          0,     '5.0/16.0', '  1/4.0'],  # rho   tau
                          ['-9.0/64.0',          0,     '5.0/64.0', ' 1/16.0']]) # rho D rho
NLO_coefs_1 = np.array(  [['-1.0/16.0', '-1.0/8.0',     '1.0/16.0', '1.0/8.0'],  # rho   tau
                          [' 3.0/64.0', ' 3.0/32 ',     '1.0/64.0', '1.0/32 ']]) # rho D rho                     


LO_terms    = [('D_I_I', 'D_I_I')]
LO_coupling = [[()]              ]
LO_coefs = CPL_string(LO_skyrme, LO_coefs_0, LO_coefs_1)

NLO_terms    = [('D_I_I', 'D_N_N'), ('lap_D_I_I', 'D_I_I')]
NLO_coupling = [[(0,1)]           , []                    ]
NLO_coefs    = CPL_string(NLO_skyrme, NLO_coefs_0, NLO_coefs_1)



#    CnablaJ(1) = B9 + 0.5_dp * B9q
#    CNablaJ(2) =      0.5_dp * B9q

#    Ct(1)   =-(B14 + 0.5_dp * B15)
#    Ct(2)   =-       0.5_dp * B15

#    Cf(1)   =-2.0_dp *(B16 + 0.5_dp*B17) ! additional factor -2 sign as the C
#    Cf(2)   =-2.0_dp *(      0.5_dp*B17) ! refer to s*F, the b to J_ij J_ij

#    Cnablas(1) = B20 + 0.5_dp * B21
#    Cnablas(2) =       0.5_dp * B21

#    CJ0    = -1.0_dp/3.0_dp * (Ct - 2.0_dp * Cf)
#    CJ1    = -0.5_dp        * (Ct - 0.5_dp * Cf)
#    CJ2    = -                (Ct + 0.5_dp * Cf)



def ProcessFunctional(fname, src, target):
    declaration = ''
    calculation = ''
    form        = ''
    printing    = ''
    calccoef    = ''
    printcoef   = ''

    for i in range(len(LO_terms)): 
        (d,c,p,cc, pc) = GenTermExpression( LO_terms[i] , LO_coupling[i], LO_coefs[i])

        declaration = declaration + d + '\n'
        calculation = calculation + c + '\n'
        printing    = printing    + p + '\n'
        calccoef    = calccoef    + cc+ '\n'
        printcoef   = printcoef   + pc+ '\n'
        
    for i in range(len(NLO_terms)): 
        (d,c,p,cc, pc) = GenTermExpression(NLO_terms[i] , NLO_coupling[i], NLO_coefs[i])

        declaration = declaration + d + '\n'
        calculation = calculation + c + '\n'
        printing    = printing    + p + '\n'
        calccoef    = calccoef    + cc+ '\n'
        printcoef   = printcoef   + pc+ '\n'
    # - - - - - - - - - - - - - - - - - - - - -
    # Substitute into the functional.f90 file.        
    dic={}
    dic['DECLARATION']    = declaration
    dic['CALCULATION']    = calculation
    dic['PRINT']          = printing
    dic['CALCCOEF']       = calccoef   
    dic['PRINTCOEF']       = printcoef    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))  

def GenTermExpression( densities, coupling, ccoef):
    #
    #
    # 
    global formatcounter, sumindices

    declaration = ''
    calculation = ''
    printing    = ''
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - 
    # Templates for the declaration, calculation and printing of an energy term.
    # 
    decl_template   = Template('real(KIND=dp) :: $CPCTE(2,2), $TERM(2,2)')
    edent_template  = Template('sum($DEN(:$IND,:),2)')
    edenq_template  = Template('$DEN(:$IND,it)')
    
    comment_template= Template(   tab + '!' + 38 * '- ' + '\n' +               \
                                  tab + '! Calculation of $TERM \n')
                                  
    doloop_template    =    tab + 'do %s = 1, 3 \n'
    enddoloop_template =    tab + 'enddo \n'
    
    calc_z_template = Template(   tab + 'Edensity = 0.0_dp \n')
    calc_a_template = Template(   tab + 'EDensity(:,3) = Edensity(:,3) + $EDENT\n')
    calc_b_template = Template(   tab + '$TERM(2,2) = 0.0_dp \n' +             \
                                  tab + 'do it=1,2 \n' +                       \
                                2*tab + 'Edensity(:,it) = Edensity(:,it) + $EDENQ \n' +          \
                                  tab + 'enddo \n')
    calc_c_template = Template(   tab + '$TERM(1,2) = $CPCTE(1,2) * sum(Edensity(:,3)) * dv \n')
    calc_d_template = Template(   tab + '$TERM(2,2) = $CPCTE(2,2) * dv & \n'+   \
                                  tab + '&'+6*tab+' * sum(Edensity(:,1) +  Edensity(:,2) ) \n')
    calc_e_template = Template(   tab + '$TERM(1,1) = $CPCTE(1,1)/$CPCTE(1,2) * $TERM(1,2)\n')
    calc_f_template = Template(   tab + '$TERM(2,1) = sum($TERM(:,2)) - $TERM(1,1) \n')
    calc_coef_template = Template(tab + '$CPCTE(1,1) = $EXP1 \n' + \
                                  tab + '$CPCTE(2,1) = $EXP2 \n' + \
                                  tab + '$CPCTE(1,2) = $CPCTE(1,1) - $CPCTE(2,1) \n' + \
                                  tab + '$CPCTE(2,2) =             2*$CPCTE(2,1) \n')       
                                 
    print_template      = Template(" print('(a20 , 3f12.3)'), '$TERM', $TERM(:,1), sum($TERM(:,1))")
    print_cpl_template  = Template(" print('(a20 , 4f12.3)'), '$CPCTE', $CPCTE")
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - 
    # See how many indices are present everywhere.
    orders      = []
    for i in range(len(densities)): 
        orders.append(OrderOfDen(densities[i])) 
    
    dic = {}
    index_encountered=0
    name = ''
    for i in range(len(densities)):
        name = name + '_' + densities[i]
    
    dic ['TERM' ] = 'E'

    
    for i in range(len(name)):
        l = name[i]
        
        if( name[i-1:i+1] == 'der'):
            dic ['TERM' ] = dic ['TERM' ] + 'd'
            for c in coupling: 
                if index_encountered in c:
                    dic ['TERM' ] = dic ['TERM' ] + sumindices[coupling.index(c)]
            index_encountered=index_encountered+1
        elif( (l.isupper() and l != 'I' and l != 'D' and l!= 'C')):
            dic ['TERM' ] = dic ['TERM' ] + l
            for c in coupling: 
                if index_encountered in c:
                    dic ['TERM' ] = dic ['TERM' ] + sumindices[coupling.index(c)]
            index_encountered=index_encountered+1
        else:   
            dic ['TERM' ] = dic ['TERM' ] + l
    
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
    for i in range(len(coupling)):
        calculation  = calculation + doloop_template%sumindices[i]
        
    calculation = calculation + calc_a_template.substitute(dic)
    calculation = calculation + calc_b_template.substitute(dic)
    for i in range(len(coupling)):
        calculation  = calculation + enddoloop_template
    calculation = calculation + calc_c_template.substitute(dic) 
    calculation = calculation + calc_d_template.substitute(dic) + '\n'
    calculation = calculation + calc_e_template.substitute(dic)
    calculation = calculation + calc_f_template.substitute(dic)
    
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # 
    dic['N'] = formatcounter
    printing = print_template.substitute(dic) 
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # 
    dic['EXP1'] = ccoef[0]
    dic['EXP2'] = ccoef[1]
    
    calccoef = calc_coef_template.substitute(dic)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # 
    printcoef = print_cpl_template.substitute(dic)    
    
    return (declaration, calculation, printing, calccoef, printcoef)
    
def OrderOfDen(density):
    #---------------------------------------------------------------------------
    # Simply checks the number of capital letters N and S in the name
    # that are not contracted.
    
    # add a dimension for every derivative and don't count the capital D or C
    test  = density.replace('_', '')
    order = test.count("der") - 1 
    test  = test.replace('der', '')
    test  = test.replace('lap', '')
    
    for letter in test:
        if(letter.isupper() and letter != 'I'):
            # Every capital letter that is not I adds an index
            order = order + 1
        if(not letter.isupper()):
            # But if the a non-capital letter follows, the index is contracted.
            order = order - 1
    return order
