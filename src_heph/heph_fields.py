#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
from string import Template
import itertools
import numpy as np
from heph_densities  import tab, ParseOperators, OrderOfDen, sumindices
import heph_functional

def GenerateFields():
    #---------------------------------------------------------------------------
    # Generate a list of fields based on list of terms in the functional. 
    #
    #---------------------------------------------------------------------------
        
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

        #  Get the operator structure of the density correctly                                
        (der, lap, left, right, coupling) = ParseOperators(den) 
        
        # Check all of the terms if they depend on the density
        fieldlist = []
        cpcte     = ''
        for term in heph_functional.Functional_terms: 
            (densities, cpl) = heph_functional.ParseDensities(term)
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
                    fieldlist.append([removed, altder, altlap,cplct])
        
        # Replace the densities in the list by the ones actually calculated
        for i in range(len(fieldlist)):
            d = fieldlist[i][0]
            (x,y,left,right,coup) = ParseOperators(d[0])
            for altden in heph_functional.Densities_needed:
                (altder, altlap, altleft, altright, altcoup) = ParseOperators(altden)
                if(altleft == left and right == altright):
                    fieldlist[i][0] = [altden]
        
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
            print fieldterm
            dic['DENSITY']  = ''
            for i in range(len(fieldterm[0])):
                    dic['DENSITY']  =  dic['DENSITY']  + fieldterm[1] * 'der_' \
                                                       + fieldterm[2] * 'Lap_' \
                                                       + fieldterm[0][i]
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



