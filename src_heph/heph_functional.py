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
        coef0 = coef0 + str(coefs_0[i,j]) + '*' + cplct[j] 
        coef1 = coef1 + str(coefs_1[i,j]) + '*' + cplct[j] 
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
LO_skyrme  =              ['t0', 't0*x0', 't3', 't3*x3']
LO_coefs_0 = np.array(   [['+3.0/8.0', '+0','+0','+0']                         # rho^2
                        , [       '+0' , '+0', '+3/48.0','+0']]  )                  
LO_coefs_1 = np.array(   [['-1.0/8.0', '-0.25','+0','+0']                    # rho^2
                        , [      '+0'  , '+0','-1/48.0', '-1/24.0']])              
#
#                                     
NLO_skyrme  =            ['t1',           't1*x1',           't2',    't2*x2']
NLO_coefs_0 = np.array(  [['+3.0/16.0',         '+0',    '+5.0/16.0', '+1/4.0'], # rho   tau
                          ['-9.0/64.0',         '+0',    '+5.0/64.0', '+1/16.0']])# rho D rho
                          
NLO_coefs_1 = np.array(  [['-1.0/16.0', '-1.0/8.0',     '+1.0/16.0', '+1.0/8.0'], # rho   tau
                          ['+3.0/64.0', '+3.0/32 ',     '+1.0/64.0', '+1.0/32 ']])# rho D rho                     

#
#
SO_skyrme  =              ['wso',  'wsoq'] 
SO_coefs_0 = np.array(   [['-0.5', '+0.25']                                      
                        , [   '+0' , '+0']] )
SO_coefs_1 = np.array(   [[   '+0' , '-0.5']                                      
                        , [   '+0' , '+0']] )  
                        
#
#
N2LO_skyrme  = [    't1n2', 't1n2 * x1n2',     't2n2', 't2n2 * x2n2']
N2LO_coefs_0 = np.array([['+9/128.0','+0', '-5/128.0', '-4/128.0'   ],  #DrhoDrho
                         ['+3/32.0' ,'+0',  '+5/32.0', '+1/8.0'     ],  # RhoQ
                         ['+3/32.0' ,'+0',  '+5/32.0', '+1/8.0'     ],  #tau^2
                         ['+6/32.0' ,'+0',  '+10/32.0', '+2/8.0'     ],  #tau_mn^2 
                         ['-6/32.0' ,'+0',  '-10/32.0', '-2/8.0'     ],
                         ['-1/32.0' , '+1/16.0', '+1/32.0' ,  '+1/16.0']]) 
N2LO_coefs_1 = np.array([['-3/128.0', '-3/64.0', '-1/128.0',  '-1/64.0'],  #DrhoDrho
                         [ '-1/32.0', '-1/16.0', '+1/32.0' ,  '+1/16.0'],  #tau^2
                         [ '-1/32.0', '-1/16.0', '+1/32.0' ,  '+1/16.0'],  #rhoQ
                         [ '-2/32.0', '-2/16.0', '+2/32.0' ,  '+2/16.0'],  #tau_mn^2
                         [ '+2/32.0', '+2/16.0', '-2/32.0' ,  '-2/16.0'],
                         [ '-1/32.0',        '+0', '+1/32.0' ,         '+0']]) 

LO_terms    = [('D_I_I', 'D_I_I'), ('D_I_I', 'D_I_I')]
LO_coupling = [[]                , []                ]
LO_coefs    = CPL_string(LO_skyrme, LO_coefs_0, LO_coefs_1)
LO_DD       = [ ''               , 'sum(D_I_I,2)**yt3a']

NLO_terms    = [('D_I_I', 'D_N_N'), ('lap_D_I_I', 'D_I_I')]
NLO_coupling = [[(0,1)]           , []                    ]
NLO_coefs    = CPL_string(NLO_skyrme, NLO_coefs_0, NLO_coefs_1)

SO_terms     = [('D_I_I', 'der_C_I_Nx1Sx1')]
SO_coupling  = [[(0,1)]                    ]
SO_coefs     = CPL_string(SO_skyrme, SO_coefs_0, SO_coefs_1)

N2LO_terms     = [('lap_D_I_I', 'lap_D_I_I'), ('D_I_I', 'D_NmNm_NnNn' )]
N2LO_coupling  = [[]                        , []]

N2LO_terms     = N2LO_terms    + [('D_N_N', 'D_N_N'),('D_N_N', 'D_N_N'), ('D_N_N', 'der_der_D_I_I') ]
N2LO_coupling  = N2LO_coupling + [    [(0,1), (2,3)],[(0,2), (1,3)]    , [(0,2), (1,3)]]

N2LO_terms     = N2LO_terms    + [('C_N_NS', 'C_N_NS') ]
N2LO_coupling  = N2LO_coupling + [[(0,3), (1,4), (2,5)]]

N2LO_coefs     = CPL_string(N2LO_skyrme, N2LO_coefs_0, N2LO_coefs_1)

def ProcessFunctional(fname, src, target):
    declaration = ''
    calculation = ''
    form        = ''
    printing    = ''
    calccoef    = ''
    printcoef   = ''
    sumtotal    = ''
    for i in range(len(LO_terms)): 
        (d,c,p,cc, pc,st) = GenTermExpression( LO_terms[i] , LO_coupling[i], LO_coefs[i], LO_DD[i])

        declaration = declaration + d + '\n'
        calculation = calculation + c + '\n'
        printing    = printing    + p + '\n'
        calccoef    = calccoef    + cc+ '\n'
        printcoef   = printcoef   + pc+ '\n'
        sumtotal    = sumtotal    + st+ '&\n&'
        
    for i in range(len(NLO_terms)): 
        (d,c,p,cc, pc,st) = GenTermExpression(NLO_terms[i] , NLO_coupling[i], NLO_coefs[i])

        declaration = declaration + d + '\n'
        calculation = calculation + c + '\n'
        printing    = printing    + p + '\n'
        calccoef    = calccoef    + cc+ '\n'
        printcoef   = printcoef   + pc+ '\n'
        sumtotal    = sumtotal    + st+ '&\n&'
    for i in range(len(SO_terms)): 
        (d,c,p,cc, pc,st) = GenTermExpression(SO_terms[i] , SO_coupling[i], SO_coefs[i])

        declaration = declaration + d + '\n'
        calculation = calculation + c + '\n'
        printing    = printing    + p + '\n'
        calccoef    = calccoef    + cc+ '\n'
        printcoef   = printcoef   + pc+ '\n'
        sumtotal    = sumtotal    + st+ '&\n&'
        
    printing = printing + "print * \n print *, 'N2LO Terms' \n print * \n"
    for i in range(len(N2LO_terms)): 
        (d,c,p,cc, pc,st) = GenTermExpression(N2LO_terms[i] , N2LO_coupling[i], N2LO_coefs[i])

        declaration = declaration + d + '\n'
        calculation = calculation + c + '\n'
        printing    = printing    + p + '\n'
        calccoef    = calccoef    + cc+ '\n'
        printcoef   = printcoef   + pc+ '\n'
        sumtotal    = sumtotal    + st+ '&\n&'
    # - - - - - - - - - - - - - - - - - - - - -
    # Substitute into the functional.f90 file.        
    dic={}
    dic['DECLARATION']    = declaration
    dic['CALCULATION']    = calculation
    dic['PRINT']          = printing
    dic['CALCCOEF']       = calccoef   
    dic['PRINTCOEF']      = printcoef 
    dic['TOTAL']          = sumtotal[:-3]
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))  

def GenTermExpression( densities, coupling, ccoef, DD=''):
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
    
    if(DD != ''):
        dic ['TERM' ] = dic ['TERM' ] + '_DD' 
  
    
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
        calculation = calculation + enddoloop_template
    if(DD != ''):
        dic['DD'] = DD
        calculation = calculation + calc_DD_template.substitute(dic)
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
    sumtotal  = sumtotal_template.substitute(dic)
    
    return (declaration, calculation, printing, calccoef, printcoef, sumtotal)
    
def OrderOfDen(density):
    #---------------------------------------------------------------------------
    # Simply checks the number of capital letters N and S in the name
    # that are not contracted.
    
    # add a dimension for every derivative and don't count the capital D or C
    test  = density.replace('_', '')
    order = test.count("der") - 1 
    test  = test.replace('der', '')
    test  = test.replace('lap', '')
    
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
