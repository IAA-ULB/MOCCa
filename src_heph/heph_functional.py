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
#
#

from string import Template


#-------------------------------------------------------------------------------
# Tab-character for the fortran code.
# 4 spaces for W.R., but I can imagine other people have different standards.
tab       ='    '
#-------------------------------------------------------------------------------
# Check how many format statements for printing the energy contributions have 
# been made. Currently there are already two format statements in the
# functional.f90 subroutine printEnergy
formatcounter = 3

def ProcessFunctional(fname, src, target):
    declaration = ''
    calculation = ''
    form        = ''
    printing    = ''
    calccoef    = ''
    printcoef   = ''

#    densities = ['D_I_I', 'D_I_I']    
#    coupling  = []
#    ccoef = [ '(3.0_dp / 8.0_dp) * t0', '-(1.0_dp / 4.0_dp) * t0 * (1.0_dp/2.0_dp + x0)']
#    
    
    densities = ['D_I_I', 'lap_D_I_I']    
    coupling  = []
    Cdrho_a='-(9.0_dp/64.0_dp) * t1 + 1.0_dp/16.0_dp * t2 * (5.0_dp/4.0_dp + x2)'
    Cdrho_b=' (1.0_dp/32.0_dp) * (3 * t1 *(0.5_dp + x1) + t2 * (0.5_dp + x2))   '
    ccoef = [Cdrho_a, Cdrho_b]

#    
    (d,c,f,p,cc, pc) = GenTermExpression( densities, coupling, ccoef)

    declaration = declaration + d
    calculation = calculation + c
    form        = form        + f
    printing    = printing    + p
    calccoef    = calccoef    + cc
    printcoef   = printcoef   + pc
    # - - - - - - - - - - - - - - - - - - - - -
    # Substitute into the functional.f90 file.        
    dic={}
    dic['DECLARATION']    = declaration
    dic['CALCULATION']    = calculation
    dic['FORM']           = form
    dic['PRINT']          = printing
    dic['CALCCOEF']       = calccoef   
    dic['PRINTCOEF']       = printcoef    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))  

def GenTermExpression( densities, coupling, ccoef):

    global formatcounter

    declaration = ''
    calculation = ''
    form        = ''
    printing    = ''
    
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - 
    # Templates for the declaration, calculation and printing of an energy term.
    # 
    decl_template   = Template('real(KIND=dp) :: $CPCTE(2,2), $TERM(2,2)')
    isosum_template = Template('sum($DEN,$IND)')
    nosum_template  = Template('$DEN(:$IND,it)')
    calc_a_template = Template(   tab + '$TERM(1,2) = $CPCTE(1,2) * sum($ISOSUM) * dv \n')
    calc_b_template = Template(   tab + '$TERM(2,2) = 0.0_dp \n ' +                          \
                                  tab + 'do it=1,2 \n' +                                     \
                                2*tab + '$TERM(2,2) = $TERM(2,2) + $CPCTE(2,2) * dv & \n'+   \
                                2*tab + '& * sum($NOSUM) \n'                             +   \
                                  tab + 'enddo \n')
    calc_c_template = Template(   tab + '$TERM(1,1) = $CPCTE(1,1)/$CPCTE(1,2) * $TERM(1,2)\n')
    calc_d_template = Template(   tab + '$TERM(2,1) = sum($TERM(:,2)) - $TERM(1,1) \n')
    calc_coef_template = Template(tab + '$CPCTE(1,1) = $EXP1 \n' + \
                                  tab + '$CPCTE(2,1) = $EXP2 \n' + \
                                  tab + '$CPCTE(1,2) = $CPCTE(1,1) - $CPCTE(2,1) \n' + \
                                  tab + '$CPCTE(2,2) =             2*$CPCTE(2,1) \n')       
                                 
    form_template       = Template(" $N  format('$TERM', 3f12.3)" )
    print_template      = Template(" print $N, $TERM(:,1), sum($TERM(:,1))")
    print_cpl_template  = Template(" print('(a20 , 4f12.3)'), '$CPCTE', $CPCTE")
    #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - 
                                 
    # See how many indices are present everywhere.
    orders      = []
    for den in densities: 
        orders.append(OrderOfDen(den)) 
    
    dic = {}
    dic ['TERM' ] = 'E'
    dic ['CPCTE'] = 'B'
    for den in densities:
         dic ['TERM' ] = dic ['TERM' ] + '_' + den
         dic ['CPCTE'] = dic ['CPCTE'] + '_' + den

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Declaration of the coupling constant and the term
    declaration = decl_template.substitute(dic) 

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Coding the isospin sums for the first (n+p) term. 
    dic['ISOSUM'] = ''
    for den in densities:
        isodic = {}
        isodic['DEN'] = den
        isodic['IND'] = orders[densities.index(den)] + 2 # +2 because of the spatial indices 
        
        dic['ISOSUM'] = dic['ISOSUM'] + isosum_template.substitute(isodic) + '*'
        
    # Take out the final '*' which should not be necessary
    dic['ISOSUM'] = dic['ISOSUM'][:-1]
    
    calculation = calc_a_template.substitute(dic)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Second term with a loop over isospin
    dic['NOSUM'] = ''
    for den in densities: 
        nosumdic = {}
        nosumdic['DEN'] = den
        
        order = orders[densities.index(den)]
        nosumdic['IND'] = ''
        for l in range(order):  
            nosumdic['IND'] = nosumdic['IND'] + ',:'
        dic['NOSUM'] = dic['NOSUM'] + nosum_template.substitute(nosumdic) + '*'
    # Take out the final '*' which should not be necessary
    dic['NOSUM'] = dic['NOSUM'][:-1]    
    
    calculation = calculation + calc_b_template.substitute(dic)
    #
    #
    calculation = calculation + calc_c_template.substitute(dic)
    calculation = calculation + calc_d_template.substitute(dic)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # 
    formatcounter = formatcounter +1 
    dic['N'] = formatcounter
    form     = form_template.substitute(dic) 
    printing = print_template.substitute(dic) 
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # 
    dic['EXP1'] = ccoef[0]
    dic['EXP2'] = ccoef[1]
    
    calccoef = calc_coef_template.substitute(dic)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # 
    printcoef = print_cpl_template.substitute(dic)    
    
    return (declaration, calculation, form, printing, calccoef, printcoef)
    
def OrderOfDen(density):
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
