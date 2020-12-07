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
# This is a non-trivial task. For example: D_N_N is not actually needed in its
# full form with an NLO functional, as the contraction D_Nm_Nm suffices. 
#-------------------------------------------------------------------------------

import itertools
import numpy as np

from string                  import Template
from src_heph.heph_densities import Densities_needed, tab, sumindices, derstring
from src_heph.heph_densities import lapstring, OrderOfDen, ParseOperators
from src_heph.heph_densities import crossindices, Storage_Mapping, Multiplicity
from src_heph.heph_densities import deriv_needed
from src_heph.heph_linechecker import *
from src_heph.heph_fields      import *

#-------------------------------------------------------------------------------
# Name of the functional_file
func_name = ''

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
assume_locality = 0

def initfunctional(fname, so):
    """
     Read the functional form from a file, populating on the way the list of 
     densities that we need to calculate. 
    """
    global Functional_terms, Densities_needed, derivative_order, func_name
    
    # Read the functional from a given file
    description = ReadFunctional(fname)

    # Check if time-reversal (or time-parity) is conserved
    if(so.timelike):
      RemoveTimeOddTerms()

    # Set the name of the functional file, without the directory structure
    func_name = '"%s"'%fname.split('/')[-1].upper()
    
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
        (deri, lapi, lefti, righti, coupi, crossi) = \
                                          ParseOperators(tempden[i],so.timelike)
        Found = False
        for j in range(len(Densities_needed)):
            (derj, lapj, leftj, rightj, coupj, crossj) = \
                                 ParseOperators(Densities_needed[j],so.timelike)
            #-------------------------------------------------------------------
            # Two densities are identical if the left- and right-operators
            # are the same.
            if(leftj == lefti and rightj == righti):
                # Signal that the density is already present
                Found = True
                # However, if there is a vector coupling in one, that is not in 
                # other, just calculate both. 
#                if(crossj != crossi):
#                    Found = False
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
        (der, lap, left, right, cpl, cross) = ParseOperators(den,so.timelike)
        ders = right.count('N')
        derivative_order = max(derivative_order, ders)

    #---------------------------------------------------------------------------
    print (' Functional form taken from file %s'%fname)
    print (' Description from file:')
    print ( description.replace('#', tab))
    print (' Number of terms     : %d'%len(Functional_terms))
    if(so.timelike):
      print (' ! Attention: terms with time-odd densities dropped. ' )
    print (' Order of derivatives: %d'%derivative_order)
    print (' Locality assumed    : %d'%assume_locality)
    print (' # Parameters        : %d'%len(paramparameters))
    #print   paramparameters
    for i in range(int(len(paramparameters)/3)):
        print ('  ', paramparameters[3*i:3*i+3])
    if(len(paramparameters)%3 != 0):
        print ('  ', paramparameters[3*(i+1):3*(i+1)+len(paramparameters)%3])

    return (description)

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
        
def ReadFunctional(fname):
    """
     Read the functional terms and the coupling coefficients from the
     functional file (fname)
    """
    global Functional_terms, Functional_pair_terms

    description      = ''
    pair_description = ''
    termsstart = 0
    with open(fname, 'r') as f:
      for line in f:
        try:
            if(len(line.split()) == 0):
                continue
            elif(line[0] == '#'):
                 #-sign indicates description of the functional                  
                 description = description + line
                 continue
            elif(line[0:6] == '!TERMS'):
                # Signal that the parameters part of the functional is over.
                termsstart  = 1   
                continue            
            elif(line[0] == '!'):
                continue
            
            if(termsstart == 0):
                # Parse the parameters of the functional
                split = line.split(';')
                for s in split:
                    paramparameters.append(s.replace(' ','').replace('\n', ''))
            elif(termsstart == 1):
                #  Start the actual terms of the functional
                split = line.split(';')
                Functional_terms.append(split[0].replace(' ', ''))
                coupling_constants_0.append(split[1].replace(' ', ''))
                coupling_constants_1.append(split[2].replace(' ', '').replace('\n', ''))   
                    
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
            print ('Problem reading the following line in the func file.')
            print (line)
            exit()      

    return description

def RemoveTimeOddTerms():
    """
      We scan through all the terms, removing all those containing time-odd 
      densities. We profit from the opportunity to check the correctness of 
      the functional, that all individual terms are time-even.
    """
    global Functional_terms, coupling_constants_0, coupling_constants_1
    global field_DD_terms, DD_rearcoefs, density_dependence

    toremove = []
    for i,term in enumerate(Functional_terms):
      (densities,coup) = ParseDensities(term)      

      timeodd= False
      totalt = +1
      for den in densities:
          t = TimeDen(den)
          totalt = totalt * t
          if(t == -1):
            timeodd = True
            
      if(totalt != +1):
        print (" A term in your functional is not time-even.")
        print ( term)
        exit()
      
      if(timeodd):
        toremove.append(i)

    tempterms   = Functional_terms
    tempcc0     = coupling_constants_0
    tempcc1     = coupling_constants_1
    tempdd      = density_dependence
    tempddrear  = DD_rearcoefs
    tempfieldDD = field_DD_terms

    Functional_terms      = []
    coupling_constants_0  = []
    coupling_constants_1  = []
    density_dependence    = []
    field_DD_terms        = {}
    DD_rearcoefs          = []
     
    for j in range(len(tempterms)):
       if (j not in toremove):
          Functional_terms.append(tempterms[j])
          coupling_constants_0.append(tempcc0[j])
          coupling_constants_1.append(tempcc1[j])
          density_dependence.append(tempdd[j])
          DD_rearcoefs.append(tempddrear[j])
          field_DD_terms[tempterms[j]] = tempfieldDD[tempterms[j]]

def ParseDensities(term): 
    """
      We deconstruct a term in the functional. 
      
      Example:
                E_D_I_Sm_Derxm_C_I_Nxm
      
      leads to
          densities : D_I_S, Der_C_I_N
          couplings : [(0,1,2)]
          
    """
    densities = []
    #---------------------------------------------------------------------------
    #  Split the input string along the underscores, removing any "_DD" 
    #  suffixes
    temp      = ''
    split     = term.replace('_DD', '').split('_') 
    for i in range(len(split)):
        if split[i][0:3] == derstring or split[i] == lapstring:
                temp  = temp + split[i] + '_'               
                
        if split[i] == 'D' or split[i] == 'C' \
                           or split[i] == 'DP' or split[i] =='CP':
                temp = temp + split[i] + '_' + split[i+1] + '_' + split[i+2]
                densities.append(temp)
                temp = ''
    #---------------------------------------------------------------------------
    # Find all couplings by looping over all possible accepted summation letters
    coupling  = []
    foundx    = []
    for l in sumindices:
        c   = ()
        cc  = ()
        ind = 0
        foundx.append(0)
        for i in range(len(term)):      
            if(term[i] == l):
                    c = c+ (ind,)
            if(term[i] in sumindices):
                    ind = ind + 1 
        if(len(c) > 0) :
                coupling.append(c)

    #---------------------------------------------------------------------------
    # Don't propagate couplings into the name that are not between left and 
    # right operators
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
    # Remove vector coupling indices
    for l in crossindices:
        for i in range(len(densities)): 
            densities[i] = densities[i].replace(l, '')
    #---------------------------------------------------------------------------
    return (densities, coupling)

def ProcessParameterization(fname, src, target):
    """
     Processing of the parameterization.f90 file to include the different 
     parameters.
    """
    decl_template = Template(  tab + 'real(KIND=dp) :: $PARAM = -123456789 \n')
    read_template = Template(        ' $PARAM,')
    print_template= Template(2*tab + 'print "(a6,2x, f10.3)", "$PARAM", $PARAM \n ')
    
    check_a_template = Template(  tab + 'if($PARAM .eq. -123456789) then \n')
    check_b_template = Template(2*tab + '   print *, "$PARAM not read from .param file." \n')
    check_c_template = Template(2*tab + '   stop \n')        
    check_d_template = Template(  tab + 'endif \n')

    decl      = ''
    readparam = ''
    printparam= ''
    checkparam= ''
    for s in paramparameters:
        dic= {}
        dic['PARAM']     = s
        
        decl      = decl + decl_template.substitute(dic)
        readparam = readparam + read_template.substitute(dic)
        printparam= printparam+ print_template.substitute(dic)
        checkparam= checkparam+ check_a_template.substitute(dic)
        checkparam= checkparam+ check_b_template.substitute(dic)
        checkparam= checkparam+ check_c_template.substitute(dic)
        checkparam= checkparam+ check_d_template.substitute(dic)
        checkparam= checkparam+ '\n'
        
    # Remove the trailing comma and add line-end
    readparam = readparam[:-1] + '\n'    
    
    dic= {}
    dic['PARAMDECL']   = decl
    dic['READPARAMS']  = readparam
    dic['PRINTPARAMS'] = printparam
    dic['CHECKPARAMS'] = checkparam
    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))  

def ProcessFunctional(fname, src, target, so):
    """
     Master routine calling the other ones to generate a functional based
     on the parsing done before.
    """

    declaration   = ''
    calculation   = ''
    form          = ''
    printing      = ''
    calccoef      = ''
    printcoef_iso = ''
    printcoef_pn  = ''
    printcoef_pair= ''
    sumtotal_even = ''
    sumtotal_odd  = ''
    pairtotal     = ''
    fieldcalc     = ''
    erear         = ''
    writing       = ''
    reading       = ''
    cleaning      = ''

    #---------------------------------------------------------------------------
    # Generate the terms in the functional
    for i in range(len(Functional_terms)): 
        (d,c,p,cc, pc_iso, pc_pn,pc_pair, st,pt,er, T)  = \
         GenTermExpression(Functional_terms[i], [coupling_constants_0[i], \
         coupling_constants_1[i]],density_dependence[i], DD_rearcoefs[i],so)

        declaration = declaration + d + '\n'
        calculation = calculation + c + '\n'
        printing    = printing    + p + '\n'
        calccoef    = calccoef    + cc+ '\n'
        printcoef_iso = printcoef_iso   + pc_iso + '\n'
        printcoef_pn  = printcoef_pn    + pc_pn  + '\n'
        printcoef_pair= printcoef_pair  + pc_pair+ '\n'
        if(T):
          # Term to be added to the time-even subtotal
          sumtotal_even = sumtotal_even    + st+ '&\n'
        else:
          sumtotal_odd  = sumtotal_odd     + st+ '&\n'

        if('P' in Functional_terms[i]):    
          pairtotal   = pairtotal   + pt+ '&\n'
        erear       = erear       + er
    
    sumtotal_even  = rreplace( sumtotal_even, '&\n', '', 1)
    sumtotal_odd   = rreplace( sumtotal_odd , '&\n', '', 1)
    pairtotal      = rreplace( pairtotal, '&\n', '', 1)

    if(len(pairtotal) == 0):    
        pairtotal = '0'

    #---------------------------------------------------------------------------
    # Generate the fields of the single-particle hamiltonian
    (fielddec, fieldcalc, fieldwrite,fieldread, fieldclean) =              \
                                                          GenerateFields(so)
    declaration = declaration + fielddec   + '\n'
    writing     = writing     + fieldwrite 
    reading     = reading     + fieldread 
    cleaning    = cleaning    + fieldclean
    #---------------------------------------------------------------------------
    # Generate the expressions for the actions of the Skyrme fields
    SkyrmeAction = ''
    for field in  Fields_needed:
      #-------------------------------------------------------------------------
      # Check if we need to symmetrize the action
      #-------------------------------------------------------------------------
      (left,right,coupling,cross) = ParseOperatorsField(field, so.timelike)
      #-------------------------------------------------------------------------
      # Generate the expression for the application of the ordinary 
      # operator structure
      # Disregard T's that are present
      if( left != right) : 
        if( 'C' in left or 'C' in right):   
          # Only symmetrize non-symmetric C's if asked for
          if(assume_locality == 1):
              SkyrmeAction = SkyrmeAction +                          \
                                  GenerateAction(field, 0, so)
          else:
              SkyrmeAction = SkyrmeAction +                          \
                                           GenerateAction(field, 1, so)
              SkyrmeAction = SkyrmeAction +                          \
                                           GenerateAction(field,-1, so)
        else:
          # Always symmetrize non-symmetric D's
          SkyrmeAction = SkyrmeAction + GenerateAction(field,+1, so)
          SkyrmeAction = SkyrmeAction + GenerateAction(field,-1, so)
      else:
          SkyrmeAction = SkyrmeAction + GenerateAction(field, 0, so)

    #---------------------------------------------------------------------------
    # Generate the expressions for the actions of the pairing fields.
    PairingAction = ''
    for field in  Pairing_Fields_needed:
      (left,right,coupling,cross) = ParseOperatorsField(field, so.timelike)
      PairingAction =PairingAction + GenerateAction(field, 0, so)

    #---------------------------------------------------------------------------
    # Now make sure all of the lines are not too long for compilation.
    declaration   = LineFormat(declaration)
    calculation   = LineFormat(calculation)
    printing      = LineFormat(printing)
    calccoef      = LineFormat(calccoef)
    printcoef_iso = LineFormat(printcoef_iso)
    printcoef_pn  = LineFormat(printcoef_pn)
    sumtotal_even = LineFormat(sumtotal_even)
    sumtotal_odd  = LineFormat(sumtotal_odd)
    fieldcalc     = LineFormat(fieldcalc)
    SkyrmeAction  = LineFormat(SkyrmeAction)
    PairingAction = LineFormat(PairingAction)
    erear         = LineFormat(erear)
    reading       = LineFormat(reading)
    writing       = LineFormat(writing)
    cleaning      = LineFormat(cleaning)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Substitute into the functional.f90 file.  
    dic={}

    dic['DECLARATION']    = declaration
    dic['CALCULATION']    = calculation
    dic['PRINT']          = printing
    dic['CALCCOEF']       = calccoef   
    dic['PRINTCOEF_ISO']  = printcoef_iso
    dic['PRINTCOEF_PN']   = printcoef_pn
    dic['PRINTCOEF_PAIR'] = printcoef_pair

    dic['TOTAL_EVEN']     = sumtotal_even
    if(len(sumtotal_odd)>1):
      dic['TOTAL_ODD']      = sumtotal_odd
    else:
      dic['TOTAL_ODD']      = '0.0d0'

    dic['TOTALPAIR']      = pairtotal
    dic['CALCFIELDS']     = fieldcalc
    dic['SKYRMEACTION']   = SkyrmeAction
    dic['PAIRINGACTION']  = PairingAction
    dic['EREAR']          = erear
    dic['FUNC_NAME']      = func_name
    dic['FIELDNUMBER']    = len(Densities_needed)
    dic['WRITEPOTENTIALS']= writing
    dic['READPOTENTIALS'] = reading
    dic['CLEANING']       = cleaning

    
    if(derivative_order == 1):
      dic['N2'] = ' '    
      dic['N3'] = '!'
    elif(derivative_order == 2): 
      dic['N2'] = ' '
      dic['N3'] = '!'
    elif(derivative_order == 3):
      dic['N2'] = '!'
      dic['N3'] = ' '
   
    if(so.timelike):
      dic['NTR'] = '!'
    else:
      dic['NTR'] = ''
    
    with open(src+fname, 'r') as template:
            with open(target+fname, 'w') as generated:
                for line in template:
                    generated.write(Template(line).substitute(dic))  

def GenTermExpression( term, ccoef, DD, DDrear, so):
    """
     Generate the expressions for the terms in the functional.
    
    """
    global sumindices

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
    pairtotal_template  = Template( 2*tab + ' & + $TERM(:,1)')
    
    calc_z_template = Template(   tab + 'Edensity = 0.0_dp \n')
    calc_a_template = Template(   tab + 'EDensity(:,3) = Edensity(:,3) $SIGN $EDENT\n')
    calc_b_template = Template(   tab + 'Edensity(:,1) = Edensity(:,1) $SIGN $EDENN\n' + \
                                  tab + 'Edensity(:,2) = Edensity(:,2) $SIGN $EDENP\n'  )
    calc_DD_template= Template(   tab + 'do m=1,3 \n' +                                 \
                                2*tab + 'EDensity(:,m) = Edensity(:,m) * $DD \n' +      \
                                  tab + 'enddo \n')
    calc_c_template = Template(   tab + '$TERM(1,2) = $CPCTE(1,2) * sum( Edensity(:,3)) * dv \n')
    calc_d_template = Template(   tab + '$TERM(2,2) = $CPCTE(2,2) * dv & \n'+           \
                                  tab + '&'+6*tab+' * sum(Edensity(:,1) + Edensity(:,2) ) \n')

    # This new division of contributions is only correct for bilinear terms.
    # If we want to include trilinear terms, we should probably explicitly 
    # construct isoscalar and isovector densities.    
    calc_e_template = Template(   tab + '$TERM(1,1) = dv * ( & \n' + 
                                  tab + '& ($CPCTE(1,2)+0.5*$CPCTE(2,2)) * sum(Edensity(:,3))) \n')
    calc_f_template = Template(   tab + '$TERM(2,1) = dv * & \n' +
                                  tab + '& (         $CPCTE(2,2) * sum(Edensity(:,1) + Edensity(:,2)) & \n' + 
                                  tab + '&   - 0.5 * $CPCTE(2,2) * sum(Edensity(:,3))               )   \n')
  
    calc_pair_a_temp= Template(   tab + 'Edensity(:,1) = Edensity(:,1) + $EDENN\n' + \
                                  tab + 'Edensity(:,2) = Edensity(:,2) + $EDENP\n'  )
    calc_pair_c_temp = Template(  tab + '$TERM(1,1) = $CPCTE(1,1) * sum( Edensity(:,1)) * dv \n')
    calc_pair_d_temp = Template(  tab + '$TERM(2,1) = $CPCTE(2,1) * sum( Edensity(:,2)) * dv \n')

    calc_coef_template  = Template(tab + '$CPCTE(1,1) = $EXP1 \n' + \
                                   tab + '$CPCTE(2,1) = $EXP2 \n' + \
                                   tab + '$CPCTE(1,2) = $CPCTE(1,1) - $CPCTE(2,1) \n' + \
                                   tab + '$CPCTE(2,2) =             2*$CPCTE(2,1) \n')       

    calc_coef_pair_temp = Template(tab + '$CPCTE(1,1) = $EXP1 \n' + \
                                   tab + '$CPCTE(2,1) = $EXP2 \n')    
                                
    write_edensity = Template( tab + ' call output_Edensity(Edensity, "$FILENAME")' )
    print_template = Template(tab +" print('(a30 , 3f15.6)'), '$TERM',    & \n"+ 
                              tab +"                           $TERM(:,1),& \n"+
                              tab +"                       sum($TERM(:,1))  \n")
    print_P_templ  = Template(tab +" print('(a30 , 30x, f15.6)'), '$TERM',& \n"+
                              tab +"                       sum($TERM(:,1))  \n")

    
    print_cpl_pn_template   = Template(tab +" print('(a30 , 4f15.6)'), '$CPCTE', $CPCTE(:,1)")
    print_cpl_iso_template  = Template(tab +" print('(a30 , 4f15.6)'), '$CPCTE', $CPCTE(:,2)")
    print_cpl_pair_template = Template(tab +" print('(a30 , 4f15.6)'), '$CPCTE', $CPCTE(:,1)")
    
    rear_template       = Template(tab +" e_rear = e_rear $REARCOEF*sum($TERM(:,2))\n")
    rear_p_template     = Template(tab +" e_rear = e_rear $REARCOEF*sum($TERM(:,1))\n")
    
    #---------------------------------------------------------------------------
    # Parse the term of the functional.
    # The result is 
    #   tempden : a list of densities that make up the term
    #   coupling: the set of scalar couplings in the term
    #
    #   Example:
    #           E_D_I_Sm_Derxm_C_I_Nxm 
    #
    #   =>  tempden :  D_I_S, Der_C_I_N
    #       coupling:  [(0,1,2)], i.e. the summation over Sm and the curl of the 
    #                  current 
    (tempden, coupling) = ParseDensities(term)
  
    # We scan the list of actually calculated densities (Densities_needed)
    # to see what contractions we have/can use.
    densities = []
    for den in tempden:
        (der,lap,left,right, coupl, cross) = ParseOperators(den,so.timelike)
        for i in range(len(Densities_needed)):
            (derref, lapref, leftref, rightref, couplref, crossref) = \
                                 ParseOperators(Densities_needed[i],so.timelike)
            if(left == leftref and right == rightref and cross == crossref):
                addden = lap*'Lap_' + der*'Der_' + Densities_needed[i]
                densities.append( addden )
                # Go back to the outer loop
                break

    # Signal back about whether this term is built out of time-odd or time-even
    # densities. Note that we don't do any checking of consistency between 
    # symmetries and terms, this is achieve somewhere else in Hephaestos.
    timerev = True
    for den in densities:
      if( TimeDen(den) < 0):
        timerev = False
    
    # We calculate the total order of all densities in the term, i.e. the total
    # number of indices.
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
    
    dic = {}
    index_encountered=0
    name = ''
    for i in range(len(densities)):
            name = name + '_' + densities[i]

    dic ['TERM' ] = term
    dic ['CPCTE'] = 'B' + dic ['TERM'][1:]    
    #---------------------------------------------------------------------------
    # We construct all possible values for all indices
    # 1. We count the number of three-length contractions in true_coupling
    nthree = 0
    for c in true_coupling:
      if(len(c) == 3):
        nthree = nthree + 1
    # We sort the coupling on length, i.e. all 2-length couplings first
    true_coupling = sorted(true_coupling, key = lambda x:len(x))

    args     = list(itertools.product(range(3), repeat=len(true_coupling)-nthree))
    vec_args = list(itertools.product(range(6), repeat=nthree))

    if(nthree == 0):
      true_args = args
    elif(len(true_coupling) == nthree):
      true_args = vec_args
    else:
      true_args = list(itertools.product(args, vec_args))
            
    declaration = decl_template.substitute(dic) 
    calculation = comment_template.substitute(dic)
    calculation = calculation = calculation + calc_z_template.substitute(dic)
    for arg in true_args: 
        dic['EDENT'] = ''
        dic['EDENP'] = ''
        dic['EDENN'] = ''
        
        sign      = +1
        prevorder =  0 
        for i in range(len(densities)):
            isodic = {}
            isodic['DEN'] = densities[i]
            
            #-------------------------------------------------------------------
            # Get the index of the density correct
            (der,lap,left,right, coupl, cross) = \
                                        ParseOperators(densities[i],so.timelike)
            
            isodic['IND'] = ''
            indices = ()
            for l in range(prevorder, prevorder + orders[i]):
                for c in true_coupling:
                    if( l in c ):
                       if(len(c) == 2):
                          mu = arg[true_coupling.index(c)]
                          indices = indices + (mu,) 
                       elif(len(c) == 3):
                          # Integer division
                          mu   = int(arg[true_coupling.index(c)]/2) 
                          # Which term of two?
                          t = arg[true_coupling.index(c)] - 2*mu
                          if(t == 1):
                            sign = sign * -1
                          nuka = Rot_ind(mu)[t] 
                              
                          temp = (mu,abs(nuka[0]), abs(nuka[1]))
                          indices = indices + (temp[c.index(l)],)
                         

            # The first indices are necessarily external derivatives
            if(der > 0):
                derind = Storage_Mapping(indices[:der])
                indices = (derind,) + indices[der:]
                
            for l in indices:
                isodic['IND'] = isodic['IND'] + ',%d'%(l+1)
           
            dic['EDENT'] = dic['EDENT'] + edent_template.substitute(isodic)+ '*'
            
            isodic['IT'] = 1
            dic['EDENN'] = dic['EDENN'] + edenq_template.substitute(isodic)+ '*'
            isodic['IT'] = 2
            dic['EDENP'] = dic['EDENP'] + edenq_template.substitute(isodic)+ '*'    
            # Take out the final '*' which should not be necessary
            prevorder = prevorder + orders[i]
      
        if(sign == +1):
          dic['SIGN'] = '+'
        else:
          dic['SIGN'] = '-'
            
        dic['EDENT'] = dic['EDENT'][:-1]
        dic['EDENN'] = dic['EDENN'][:-1]  
        dic['EDENP'] = dic['EDENP'][:-1]  
        
        calculation = calculation +'\n' + tab + '! indices = ' + str(arg) + '\n'

        if('P' not in term ):
          # Ordinary mean-field densities only in the term
          calculation = calculation + calc_a_template.substitute(dic)
          calculation = calculation + calc_b_template.substitute(dic)
        else:
          # Pairing mean-field densities in the term.
          calculation = calculation + calc_pair_a_temp.substitute(dic)
            
    calculation = calculation + '\n'
    if(DD != ''):
        dic['DD'] = DD  
        calculation = calculation + calc_DD_template.substitute(dic)

    dic['FILENAME'] ='edensities/' + dic['TERM'] + '.dat'
    #calculation = calculation + write_edensity.substitute(dic) + '\n'
 
    if('P' not in term):  
      # Ordinary mean-field densities
      calculation = calculation + calc_c_template.substitute(dic) 
      calculation = calculation + calc_d_template.substitute(dic) + '\n'
      calculation = calculation + calc_e_template.substitute(dic)
      calculation = calculation + calc_f_template.substitute(dic)

      printing = print_template.substitute(dic) 

    else:
      # Pairing mean-field densities.
      calculation = calculation + calc_pair_c_temp.substitute(dic) 
      calculation = calculation + calc_pair_d_temp.substitute(dic) + '\n'

      printing = print_P_templ.substitute(dic) 

    calculation = calculation + end_comment 
    
    dic['EXP1'] = ccoef[0]
    dic['EXP2'] = ccoef[1]
    
    calccoef  = calc_coef_template.substitute(dic)
    
    if('P' not in term):
        # Ordinary mean-field densities
        printcoef_pn  = print_cpl_pn_template.substitute(dic)    
        printcoef_iso = print_cpl_iso_template.substitute(dic)
        printcoef_pair= ''
    else:
        printcoef_pn  = ''
        printcoef_iso = ''
        printcoef_pair= print_cpl_pair_template.substitute(dic) 
    
    sumtotal  = sumtotal_template.substitute(dic)
    pairtotal = pairtotal_template.substitute(dic)
    
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

        if('P' not in term):
          erear = rear_template.substitute(dic)
        else:
          erear = rear_p_template.substitute(dic)
        
    return (declaration, calculation, printing, calccoef, printcoef_iso, 
              printcoef_pn, printcoef_pair, sumtotal, pairtotal, erear, timerev)    
    
def rreplace(s, old, new, occurrence):
     li = s.rsplit(old, occurrence)  
     return new.join(li)
