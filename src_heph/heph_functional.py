#--------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#--------------------------------------------------------------------

from string import Template
import itertools
import numpy as np
import heph_fields              # The file for the fields
import heph_densities           # The file for the densities

#------------------------------------------------------------------------------
# TODO
#
#
#-------------------------------------------------------------------------------
# Indices over which sums are supposed to go in both the FORTRAN code and the 
# naming scheme.
sumindices    = ['m', 'n', 'k', 'l']
#-------------------------------------------------------------------------------
# Tab-character for the fortran code.
# 4 spaces for W.R., but I can imagine other people have different standards.
tab           = '    '

#
#
#
Functional_terms     = []
Functional_coupling  = []
coupling_constants_0 = [] 
coupling_constants_1 = [] 

#-------------------------------------------------------------------------------
# Densities_needed contains all of the densities that Hephaestos deems are 
# necessary to calculate, based on the terms in the functional. 
Densities_needed = []
Deriv_den_needed = []
Lapla_den_needed = []

def initfunctional(fname):
        global Functional_terms, Densities_needed

        # Read the functional from a given file
        description = ReadFunctional(fname)
        
        # Parse all of the densities needed
        temp      = []
        for term in Functional_terms:
                (dens, coup) = ParseDensities(term)
                for d in dens: 
                        temp.append(d)
                        Functional_coupling.append(coup)

        # Now scan the list for duplicates:     
        for i in range(len(temp)): 
                Found = False
                (deri, lapi, lefti, righti) = heph_densities.ParseOperators(temp[i])
                                        
                for j in range(len(Densities_needed)):
                        (x, y, leftj, rightj) = heph_densities.ParseOperators(Densities_needed[j])  
                        derj = Deriv_den_needed[j]
                        lapj = Lapla_den_needed[j]
                        if(lefti == leftj and righti == rightj): 
                                # Already added this density   
                                deri =   max(deri, derj)
                                Deriv_den_needed[j] = deri
                                lapi = max(lapi, lapj)
                                Lapla_den_needed[j] = lapi
                                Found = True
                                continue
                if(not Found):
                        Densities_needed.append(temp[i])
                        Deriv_den_needed.append(deri)
                        Lapla_den_needed.append(lapi)  
                                       
                        
        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -' 
        print ' Functional taken from file %s'%fname
        print ' Description from file:'
        print  description.replace('#', tab)
        print ' Number of terms:     %d'%len(Functional_terms)
        print ' Number of densities: %d'%len(Densities_needed)
        print '- - - - - - - - - - - - - - - - - - - - - - - - - - - - -'

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
                (d,c,p,cc, pc,st) = GenTermExpression(Functional_terms[i], [coupling_constants_0[i], coupling_constants_1[i]])

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
        dic['CALCFIELDS']     = fieldcalc
        with open(src+fname, 'r') as template:
                with open(target+fname, 'w') as generated:
                    for line in template:
                        generated.write(Template(line).substitute(dic))  

def ReadFunctional(fname):
        #
        # Read the functional terms and the coupling coefficients from file fname.
        #
        global Functional_terms

        description=''
        with open(fname, 'r') as f:
                for line in f: 
                        if(line[0] != '#'):
                                split = line.split(',')
                                Functional_terms.append(split[0].replace(' ', ''))        
                                coupling_constants_0.append(split[1].replace(' ', ''))
                                coupling_constants_1.append(split[2].replace(' ', ''))   
                        else:
                                description = description + line     
        

        return(description)

def ParseDensities(term): 
        #------------------------------------------------------------------------
        # From a term in the functional represented by a string, parse all of the
        # densities that make it up.
        #-----------------------------------------------------------------------
        densities = ()
        # Split along C and D-s
        nocoupling= term
        for l in sumindices:
                nocoupling = nocoupling.replace(l, '')        
        split     = nocoupling.split('_')
        temp      = ''
        for i in range(len(split)):
                if split[i] == 'der' or split[i] == 'Lap':
                        temp  = temp + split[i] + '_'               
                        
                if split[i] == 'D' or split[i] == 'C':
                        temp = temp + split[i] + '_' + split[i+1] + '_' + split[i+2]
                        densities = densities + (temp,)
                        temp = ''

        # Find the coupling
        coupling  = []
        for l in sumindices:
                c   = ()
                ind = 0
                for i in range(len(term)):      
                        if(term[i] == l):
                                c = c+ (ind,)
                                     
                        if(term[i] in sumindices):
                                ind = ind + 1 
                if(len(c) > 0) :
                        coupling.append(c)
        return (densities, coupling)


def GenTermExpression( term, ccoef, DD=''):

    global  sumindices

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
    
    
    (densities, coupling) = ParseDensities(term)

    orders                = []
    
    for i in range(len(densities)): 
        orders.append(heph_densities.OrderOfDen(densities[i])) 
    
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
