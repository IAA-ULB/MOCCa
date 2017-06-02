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

#############################################################################################
# TODO
#     C  Add support for traces and vector products: might be better from the top-down.
#############################################################################################

#-------------------------------------------------------------------------------
# Tab-character for the fortran code.
# 4 spaces for W.R., but I can imagine other people have different standards.
tab           = '    '

#-------------------------------------------------------------------------------
# Indices over which sums are supposed to go in both the FORTRAN code and the 
# naming scheme.
sumindices    = ['m', 'n', 'k', 'l']

def ParseDensities(term): 
        densities = ()
        # Split along C and D-s
        split = term.split('_')
        
        temp = ''
        for i in range(len(split)):
                if split[i] == 'der' or split[i] == 'lap':
                        temp  = temp + split[i] + '_'               
                        
                if split[i] == 'D' or split[i] == 'C':
                        temp = temp + split[i] + '_' + split[i+1] + '_' + split[i+2]
                        densities = densities + (temp,)
                        temp = ''
        return densities

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

def ParseOperators(density):    
        
        left  = ''
        right = '' 

        der   = density.count('der_')
        lap   = density.count('lap_')

        split = density.replace('der_', '').replace('lap_', '').split('_')
        
        left = split[1]
        right= split[2]

        if(density[0] == 'C'):      
                right = 'C' + right

        return(der, lap, left, right)

def GenerateFields(Eterms):
        #---------------------------------------------------------------------
        # Generate a list of fields based on list of terms in the functional. 
        #
        #---------------------------------------------------------------------
              
        global sumindices,tab

        field_decl_temp = Template(  tab + 'real(KIND=dp), allocatable :: $FIELD(:$DECLIND,:) \n')
        field_allo_temp = Template(  tab + 'if(.not.allocated($FIELD)) then\n'+\
                                   2*tab + 'allocate($FIELD(mv$ALLOCIND,2)) \n'   +\
                                     tab + 'endif \n' + \
                                     tab + '$FIELD = 0.0 \n')
        field_calc_temp    = Template(   tab + '$FIELD(:$IND,it) =  & \n')
        field_calc_b_temp  = Template( 3*tab + '& $SIGN sum($CPLCTE(:,2)) * $DENSITY(:$IND,it)   & \n') 
        field_calc_c_temp  = Template( 3*tab + '& $SIGN      $CPLCTE(2,2) * $DENSITY(:$IND,3-it) & \n') 

        doloop_template    = tab + 'do %s = 1, 3 \n'
        enddoloop_template = tab + 'enddo \n'
        
        densities = []
        cplcts    = []
        for term in Eterms:
                temp = ParseDensities(term)
                densities.append(temp)
                cplcts.append(term.replace('E_', 'B_'))

        noder_densities=[]
        unique_densities=[]
        for t in densities:
                dt = []
                for den in t: 
                       d = den.replace('der_', '').replace('lap_', '')
                       dt.append(d)
                       unique_densities.append(d)
                noder_densities.append(dt)
                
        unique_densities = list(set(unique_densities))
        
        #----------------------------------------------------------------------
        # For every unique density encountered, we need to figure out the field
        # and the action of the field. 
        FIELDCALC   = ''
        declaration = ''
                
        for den in unique_densities:
                # Name the field correctly
                dic = {}
                dic['FIELD'] = den.replace('D', 'F').replace('C', 'G')
                #  Get the operator structure of the density correctly                                
                (der, lap, left, right) = ParseOperators(den) 
                
                # Check all of the terms if they depend on the density
                fieldlist = []
                cpcte     = ''
                for term in densities: 
                        for i in range(len(term)):
                                altden = term[i]
                                (altder, altlap, altleft, altright) = ParseOperators(altden)
                                if(altleft == left and altright == right):
                                        #  Add the term to the fieldlist for this density, 
                                        #  and additionnally mentioning the number of external 
                                        #  derivatives and laplacians
                                        removed = []
                                        for j in range(len(term)):
                                                if i != j :
                                                        removed.append(term[j])
                                        cplct = cplcts[densities.index(term)]
                                        fieldlist.append((removed, altder, altlap,cplct))
                                                                      
                # Create the expression for the field
                dic['IND']     = ''
                dic['ALLOCIND'] = ''
                dic['DECLIND'] = ''
                for k in range(OrderOfDen(den)):
                        dic['IND']      = dic['IND']      + ',' + sumindices[k]
                        dic['ALLOCIND'] = dic['ALLOCIND'] + ',3' 
                        dic['DECLIND']  = dic['DECLIND']  + ',:' 

                declaration  = declaration + field_decl_temp.substitute(dic)
                FIELDCALC = FIELDCALC + field_allo_temp.substitute(dic)
                for k in range(OrderOfDen(den)):
                        FIELDCALC = FIELDCALC + doloop_template%sumindices[k]
                        
                FIELDCALC = FIELDCALC + field_calc_temp.substitute(dic)
                for fieldterm in fieldlist:
                                                
                        dic['DENSITY']  = ''
                        for i in range(len(fieldterm[0])):
                                dic['DENSITY']  =  dic['DENSITY']  + fieldterm[1] * 'der_' + fieldterm[0][i]
                        if( fieldterm[1]%2 == 0) :
                                dic['SIGN']     =  '+'
                        else:
                                dic['SIGN']     =  '-'
                        dic['CPLCTE']   =  fieldterm[3]
                        FIELDCALC = FIELDCALC + field_calc_b_temp.substitute(dic)
                        FIELDCALC = FIELDCALC + field_calc_c_temp.substitute(dic)
                FIELDCALC = FIELDCALC[:-4] + '\n'
                # Create the expression for the field
                for k in range(OrderOfDen(den)):
                        FIELDCALC = FIELDCALC + enddoloop_template
                        
        return(declaration,FIELDCALC)



