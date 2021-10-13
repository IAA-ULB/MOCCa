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

# List of fields needed 
Fields_needed         = []
Pairing_Fields_needed = []

def initfields(so):

  #-----------------------------------------------------------------------------
  # Go over the needed densities and the functional terms and check whether
  # we have enough derivatives to calculate the fields. 
  for term in src_heph.heph_functional.Functional_terms:
      (densities, cpl) = src_heph.heph_functional.ParseDensities(term)
      # Count the number of derivatives needed in this term
      totalder = 0
      totallap = 0
      for den in densities:
          (der,lap,left,right, coupling, cross) = \
                                                 ParseOperators(den,so.timelike)
          totalder = totalder + der
          totallap = totallap + lap
      
      # Now see that for all densities in this term, the minimum number
      # of derivatives is the total one 
      for i in range(len(densities)):
          den = densities[i]
          (der,lap,left,right,coupling,cross) = ParseOperators(den,so.timelike)
          for j in range(len(Densities_needed)):
              altden = Densities_needed[j]
              (altder, altlap, altleft, altright, altcoupling, altcross)     \
                                            = ParseOperators(altden,so.timelike)    
              if(altleft == left and altright == right):
                  # Set minimum derivatives
                  deriv_needed[j].append((totallap, totalder))
                  
  
  src_heph.heph_functional.PruneDeriv_needed()
  
      
def GenerateFields(so, oldso):
    """
     Generate a list of fields based on list of terms in the functional. 
    """
        
    global sumindices,tab
    
    import src_heph.fortran_templates.GenerateFields_templates as ts 

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

        # We check if there is a change of spatial symmetries between the 
        # new (so) and old (oldso) symmetry options        
        spatial = False
        for k in range(3):
          if(so.ReduceAxes[k] != oldso.ReduceAxes[k]):
            spatial = True
        
        if(spatial):
          # Tantalus will essentially not read potentials from file
          dic['UNDOREAD'] = ' '
        else:
          # Tantalus will read potentials correctly from file
          dic['UNDOREAD'] = '!'
        
        if('P' not in den): 
          Fields_needed.append(dic['FIELD'])
        else:
          Pairing_Fields_needed.append(dic['FIELD'])
        #-----------------------------------------------------------------------
        #  Get the operator structure of the density correctly                                
        (der,lap,left, right, coupling,cross) = ParseOperators(den,so.timelike) 
        #-----------------------------------------------------------------------
        # Check all of the terms if they depend on the density
        fieldlist = []
        cpcte     = ''

        fieldclean = fieldclean + '\n' + ts.clean.substitute(dic)
        fieldclean = fieldclean + '\n' + ts.clean_b.substitute(dic)

        #-----------------------------------------------------------------------
        for term in src_heph.heph_functional.Functional_terms: 
            (densities, cpl)=src_heph.heph_functional.ParseDensities(term)
            #-------------------------------------------------------------------
            # Replace the densities in the list by the ones actually calculated
            for i in range(len(densities)):
                (x,y,l2,r2,c2,cr2) = ParseOperators(densities[i],so.timelike)
                for altden in src_heph.heph_functional.Densities_needed:
                    (altder, altlap, altleft, altright, altcoup, altcross) = \
                                              ParseOperators(altden,so.timelike)
                    if(altleft == l2 and r2 == altright and cr2 == altcross):
                        densities[i] = y*'Lap_' +               \
                                       x*'Der_' +               \
                                       altden
            #-------------------------------------------------------------------
            # Remove all the mentions of couplings inside the density if only
            # contractions are calculated.
            altterm = term
            for c in cpl: 
                for i in range(len(densities)):
                    if(OrderOfDen(densities[i]) != OrderOfDen(densities[i], contract=False)):
                      for s in sumindices:
                          if(densities[i].count(s) == 2):
                            # Replace internal couplings
                            altterm = altterm.replace(s,'')

            (rubbish, cpl) = src_heph.heph_functional.ParseDensities(altterm)
            #-------------------------------------------------------------------
            # Check if the term contains this density
            startind = 0
            for i in range(len(densities)):
                altden = densities[i]
                (altder, altlap, altleft, altright, altcoup, altcross) \
                                            = ParseOperators(altden,so.timelike)
                
                if((altleft == left) and (altright == right)):
                    #  Add the term to the fieldlist for this density, 
                    #  and additionally mentioning the number of external 
                    #  derivatives and laplacians

                    # First we count how much indices are accounted for by the 
                    # other densities in the term
                    removed    = []
                    removedsum = 0
                    for j in range(len(densities)):
                        if i != j :
                            removed.append(densities[j])
                    for j in range(i):
                            removedsum = removedsum + OrderOfDen(densities[j])
                    #-----------------------------------------------------------
                    # From a term
                    # 
                    #      E_C/D_L'_R'_C/D_L_R
                    #
                    # with some prescribed coupling between the indices of L,R
                    #  and L', R' (and possibly more densities...)
                    #
                    # Now we derive this term with respect to a density, which 
                    # is not necessarily the first in the list, i.e. C/D_L_R.
                    # We are then looking for the contribution to the 
                    # corresponding 
                    #
                    #   F/G_L_R  = ... +  CC  C/D_L'_R' + ...
                    # 
                    # where CC is a coupling constant. 
                    #-----------------------------------------------------------
                    newcpl = []
                    shift = OrderOfDen(den) + altder
                    for c in cpl:
                        nc = ()
                        for k in c:     
                            if( k < startind): 
                                nc = nc + (k+shift,)
                            elif(k >= startind + OrderOfDen(altden) ):
                                nc = nc + (k,)
                            elif(k >= startind + altder):
                                nc = nc + (k-removedsum-altder,)
                            else:
                                nc = nc + (k-removedsum+OrderOfDen(den),)
                        newcpl.append(nc)        

                    ind   = src_heph.heph_functional.Functional_terms.index(term)
                    cplct = cplcts[ind]
                    dden  = src_heph.heph_functional.density_dependence[ind]

                    fieldlist.append([removed, altder,altlap,cplct,newcpl,dden,0])
                startind = startind + OrderOfDen(altden) 
                
            #-------------------------------------------------------------------
            # Now check if there are density dependences in this term that 
            # involve this density
            dd    = src_heph.heph_functional.field_DD_terms[term]
            ind   = src_heph.heph_functional.Functional_terms.index(term)
            cplct = cplcts[ind]
            if(dd[0] == den):
                fieldlist.append([densities, altder, altlap,cplct,cpl,dd[1],1])

        # Create the expression for the field
        dic['ALLOCIND']= ''
        dic['DECLIND'] = ''
        for k in range(OrderOfDen(den)):
            dic['ALLOCIND'] = dic['ALLOCIND'] + ',3' 
            dic['DECLIND']  = dic['DECLIND']  + ',:'
       
        declaration  = declaration + ts.field_decl.substitute(dic)
        declaration  = declaration + ts.fhist_decl.substitute(dic)

        fieldread    = fieldread   + ts.field_read_a.substitute(dic)
        fieldread    = fieldread   + ts.field_allo_b.substitute(dic)
        fieldread    = fieldread   + ts.field_read_b.substitute(dic)
        fieldread    = fieldread   + ts.field_read_c.substitute(dic)
        fieldread    = fieldread   + ts.field_read_d.substitute(dic)
        fieldread    = fieldread   + ts.field_read_e.substitute(dic)
        
        fieldwrite   = fieldwrite  + ts.field_write_a.substitute(dic)
        fieldwrite   = fieldwrite  + ts.field_write_b.substitute(dic)
        
        FIELDCALC    = FIELDCALC + ts.field_line.substitute(dic)
        FIELDCALC    = FIELDCALC + ts.field_calc_a_start.substitute(dic)
        FIELDCALC    = FIELDCALC + ts.field_allo.substitute(dic)
        FIELDCALC    = FIELDCALC + ts.field_hist.substitute(dic)
        FIELDCALC    = FIELDCALC + ts.isoloop.substitute(dic)

        args = list(itertools.product(range(3), repeat=OrderOfDen(den)))
        for arg in args:   
           # get the indices of the field correct
           dic['IND']     = ''
           for k in arg:
                  dic['IND'] = dic['IND'] + ',%d'%(k+1)

           fieldread = fieldread + ts.field_transfo.substitute(dic)

        fieldread = fieldread + ts.field_transfo.substitute(dic)
        fieldread = fieldread  + ts.field_read_f.substitute(dic)

        for fieldterm in fieldlist:
             #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
             # We start constructing individual contributions to the fields  
             #
             # A statement in the code will have the following form formally
             #
             #   F/G_{mu nu .... kappa } =  F/G_{mu nu .... kappa }
             #    +/- C sum_{a,b,c} D_{mu nu ... kappa a b c  }
             #
             # where the we have dropped the isospin index and the C is some
             # coupling constant. 
             # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
             # The first order of business is to determine the number of free
             # indices there are in the field F/G, i.e. how many INDEPENDENT 
             # indices there are in total in the entire expression above
             NumberOfIndices = 0
             nthree = 0
             for c in fieldterm[4]:
                if(len(c) == 2):
                  NumberOfIndices = NumberOfIndices + 1
                elif(len(c) == 3):
                  nthree          = nthree +1  
             
             #------------------------------------------------------------------
             # arguments for all the indices
             args = list(itertools.product(range(3), repeat=NumberOfIndices))
             vec_args = list(itertools.product(range(6), repeat=nthree))

             if(nthree == 0):
                true_args = args
             elif(NumberOfIndices == 0 ):
                true_args = vec_args
             else:
                true_args = itertools.product(args, vec_args)

             for arg in true_args:
                 # get the indices of the field (i.e. the lhs above) correct
                 dic['IND']     = ''
                 sign           = +1

                 #--------------------------------------------------------------
                 # Now we do something wacky: we detect if   
                 # we have changed the ordering of indices in a 
                 # vector product by partial integration
                 #
                 # For example:
                 #
                 #      E_C_I_Nm_Derxm_D_I_Sxm
                 # 
                 # contributes to F_I_S. It gives rise to a fieldterm
                 # like 
                 #    F_I_Sxm  \sim Derxm_C_I_Nm
                 # 
                 # I detect this here in a simple way: if in a three-coupling 
                 # ( a,b,c ) the two largest elements are ordered correctly, 
                 # things are okay. If not, we have changed the ordering of a 
                 # vector product  by partial integration.
                 #--------------------------------------------------------------

#                 for c in fieldterm[4]:
#                  if(len(c) == 3):
#                    nc = list(c)
#                    sign = sign * perm_parity(nc)
#                    print (den, sign)
                 #--------------------------------------------------------------
                 for k in range(OrderOfDen(den)):
                    for c in fieldterm[4]:
                        if k in c:
                            if(len(c) == 2):
                              mu = arg[fieldterm[4].index(c)]
                              dic['IND'] = dic['IND'] + ',%d'%(mu+1)
                            elif(len(c) == 3):
                              # Integer division
                              mu   = int(arg[fieldterm[4].index(c)]/2) 
                              # Which term of two?
                              t = arg[fieldterm[4].index(c)] - 2*mu
                              nuka = Rot_ind(mu)[t] 
                              if(t == 1):
                                sign = sign * -1                              
                              temp = (mu,abs(nuka[0]), abs(nuka[1]))
                              dic['IND'] = dic['IND'] + ',%d'%(temp[c.index(k)]+1)

                 FIELDCALC = FIELDCALC + ts.field_calc_b_start.substitute(dic)
                 
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
                    
                    dic['DENSITY'] = dic['DENSITY'].replace('Der_', '')
                    dic['DENSITY'] = dic['DENSITY'].replace('Lap_', '')
                    dic['DENSITY'] = lapcount * 'Lap_' \
                                   + dercount * 'Der_' + dic['DENSITY']

                    dic['CPLCTE']   =  fieldterm[3]
                    #-----------------------------------------------------------
                    # Get the indices of the density on the rhs.
                    indices = ()
                    for k in range(lastorder,lastorder+OrderOfDen(dic['DENSITY'])):
                      for c in fieldterm[4]:
                        if k in c:   
                          if(len(c) == 2):    
                            mu      = arg[fieldterm[4].index(c)]
                            indices = indices + (mu,)
                          elif(len(c) == 3):
                            # Integer division
                            mu   = int(arg[fieldterm[4].index(c)]/2) 
                            # Which term of two?
                            t = arg[fieldterm[4].index(c)] - 2*mu
                            if(t == 1):
                              sign = sign * -1

                            nuka = Rot_ind(mu)[t] 
                                  
                            temp = (mu,abs(nuka[0]), abs(nuka[1]))
                            indices = indices + (temp[c.index(k)],)                                 
                    #-----------------------------------------------------------                      
                    # Put an extra sign for every partial integration of a nabla
                    # needed
                    if( fieldterm[1]%2 != 0):
                        localsign = sign * (-1)
                    else:
                        localsign = sign
  
                    if(localsign > 0):
                      dic['SIGN']     =  '+'
                    else:
                      dic['SIGN']     =  '-'

                    #-----------------------------------------------------------
                    # The first indices are necessarily external derivatives
                    if(dercount > 0):
                        derind  = Storage_Mapping(indices[:dercount])
                        indices = (derind,) + indices[dercount:]
                    
                    dic['DENIND']      = ''
                    dic['SUMIND']      = 2 
                    for l in indices:
                        dic['DENIND'] = dic['DENIND'] + ',%d'%int(l+1)

                    dic['EXPR1'] = dic['EXPR1'] + ts.field_calc_den_a.substitute(dic)
                    dic['EXPR2'] = dic['EXPR2'] + ts.field_calc_den_b.substitute(dic)
                   
                    if(len(dic['DD']) >0):
                        dic['EXPR3'] = dic['EXPR3'] + ts.field_calc_den_c.substitute(dic)
                        lastorder = lastorder + OrderOfDen(dic['DENSITY'])

                 if('P' not in dic['DENSITY']):
                    # Ordinary mean-field densities; coupling constants are
                    # the isoscalar and isovector ones
                    FIELDCALC = FIELDCALC + ts.field_calc_b.substitute(dic)
                    FIELDCALC = FIELDCALC + ts.field_calc_c.substitute(dic)
                    if(fieldterm[6] == 1):
                       FIELDCALC = FIELDCALC + ts.field_calc_d.substitute(dic)
                    FIELDCALC = FIELDCALC[:-4] + '\n \n'
                 else:
                    # Pairing densities, coupling constants are pn ones.  
                    FIELDCALC = FIELDCALC + ts.field_pair_a.substitute(dic)       
                    if(fieldterm[6] == 1):
                      FIELDCALC = FIELDCALC + ts.field_pair_b.substitute(dic)
                    FIELDCALC = FIELDCALC[:-4] + '\n \n'
                                
        
        FIELDCALC    = FIELDCALC + ts.isoloop_end
        FIELDCALC    = FIELDCALC + ts.field_line.substitute(dic)

    return(declaration, FIELDCALC, fieldwrite, fieldread, fieldclean)

def GenerateAction(field, symmetrize, so):
    """
     Generate the action of a field. 
    
        field       : name of the field, for example F_I_I
    
        symmetrize  : whether or not to generate a term for a symmetrised field  
                      (  0) generate action "as is"  
                      (+-1) generate one out of two terms of the action
    
        so          : a set of symmetry options 
    """
    
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
    # First, figure out whether there is an antilinear, antihermitian symmetry
    # that is conserved.

    #---------------------------------------------------------------------------
    # Parse the field under consideration
    (left,right,coupling,cross) = ParseOperatorsField(field, so.timelike)
    #---------------------------------------------------------------------------
    # If symmetrize is non-zero, change left <-> right and the couplings
    # accordingly
    if(symmetrize == -1 ):
        switch_field = field.split('_')
        R = switch_field[2]
        L = switch_field[1]
        # We put all spin indices on the r.h.s.
        # This allows us to 
        # 1. Be more efficient with "on the fly" derivatives
        #    (They can be grouped more)   
        # 2. make it easier to figure out the symmetries of the intermediate
        #    functions for the "on the fly" derivatives
        if('S' in R): 
          for l in sumindices:
            if('S' + l in R):
               R = R.replace('S'+l, '')
               L = L + 'S' + l
          if('S' in R):
            R = R.replace('S', '')
            L = L + 'S'
        switch_field = switch_field[0]+'_'+R+'_'+ L
        (left,right,coupling,cross) \
                                 = ParseOperatorsField(switch_field,so.timelike)
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

#-------------------------------------------------------------------------------
# These multiplicities are a good idea and often useful, but as implemented is
# not always correct.
    #-------------------------------------------------------------------
    # Go over the uncontracted right-indices and get the independent  
    # components, and the multiplicities. 
    #larg_stor = []
    #larg_new  = []
    #multiplicities         = []    
    #for true_larg in largs: 
    #    l_stor = Storage_Mapping(true_larg[:LeftOperator.derorder])
    #    mult   = Multiplicity(true_larg[:LeftOperator.derorder])
    #    l_new  = (l_stor,) + true_larg[LeftOperator.derorder:]
   # 
   #     Found = False
   #     for new in larg_stor:
   #         if(new == l_new):
   #             Found = True
   #     if(not Found):
   #         larg_stor.append(l_new)
   #        larg_new.append(true_larg)
   #         multiplicities.append(mult)
    
   # largs = larg_new
#-------------------------------------------------------------------------------
    
    for true_larg in largs:
        
        # reset temp to 0
        expression = expression + temp_ini
        
        # Get the multiplicity correct
        #m = multiplicities[largs.index(true_larg)]
        
        #if( m != 1):
        #    dic['LMULT'] = str(m) + ' * '
        #else:
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

                if(symmetrize != -1):
                    # Original ordering of indices
                    for l in range(LeftOperator.dimension):
                        Found = False
                        for c in coupling:
                            if(l in c):
                                Found = True
                        if(not Found):
                            dic['FIELDIND']= dic['FIELDIND']  \
                                               + ',%d'%int(abs(true_larg[l])+1)
                    for r in range(len(rarg)):
                        dic['FIELDIND']= dic['FIELDIND']      \
                                               + ',%d'%int(abs(rarg[r])+1)
                else:
                    # Switching (L<->R) of all derivative indices, except for 
                    # the spin operator, which we always take to the right
                    if('S' in field):
                        offset = 1
                    else:
                        offset = 0
                    for r in range(len(rarg)-offset):
                        dic['FIELDIND']= dic['FIELDIND']      \
                                               + ',%d'%int(abs(rarg[r])+1)
                    for l in range(LeftOperator.dimension):
                        Found = False
                        for c in coupling:
                            if(l in c):
                                Found = True
                        if(not Found):
                            dic['FIELDIND']= dic['FIELDIND']  \
                                               + ',%d'%int(abs(true_larg[l])+1)
                    # Reinsert the spin index
                    if(offset == 1):
                        dic['FIELDIND']= dic['FIELDIND']      \
                                               + ',%d'%int(abs(rarg[-1])+1)
                #---------------------------------------------------------------                     
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
                        dic['RIND'] =  dic['RIND']  + ',%d'%int(rarg_stor +1 )                        
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
                # that still need to be performed after this derivative. This
                # works, since we have forced the application of the spin 
                # operators to be part of the right operator
                for l in range(lorder, LeftOperator.derorder):
                    if true_larg[- l - offset - 1] +1  == direc:
                        sym = -sym
                
                for k in range(4):
                    dic['RCOMP'] = ''
                    for l in range(lorder):
                        dic['RCOMP']  = str(true_larg[-l-offset-1] +1) \
                                                            + ',' + dic['RCOMP'] 
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
    
def ParseOperatorsField(field, timelike):    
    """
     Parse the operators that are used to construct the action of field.
     This routine is very close to its analogue for the densities in 
     heph_densities.py, but has a subtle twist to it concerning calculations
     with conserved T or PT.

     timelike indicates if there is a conserved antilinear, antihermitian 
     symmetry that can be used to simplify the calculation. 
    
     For ordinary fields, this has no impact. 

    """    
    
    #
    # timelike indicates if there is a conserved antilinear, antihermitian 
    # symmetry that can be used to simplify the calculation. 
    #
    # For ordinary fields, this has no impact. 
    #
    # For pairing fields however, this information is vital to generate 
    # correct code. Without using any symmetry, a pairing density is represented
    # here with an explicit time-reversal operator
    #  
    #   \tilde{\rho} = DP_L_R 
    #   =\sum_{ij,s} \kappa_{ij} s \int d^3 r [L  \psi_j] (r,-s) [R \psi_i](r,s)
    #   =\sum_{ij,s} \kappa_{ij}   \int d^3 r [TL \psi_j]^*(r,s) [R \psi]_i(r,s)
    #
    # so we could have written
    #
    #     DP_L_R => D_TL_R
    #
    # and hence, the corresponding field is 
    #
    #     FP_L_R => F_TL_R
    #
    # Now, pairing fields are only needed in a mean-field code when calculating 
    # the pairing gaps. If a timelike symmetry is not conserved, their
    # calculation can be achieved in exactly the same way as for ordinary
    # fields. If a timelike symmetry is conserved however, we only calculate
    # part of the pairing gaps, i.e. only gaps of the form
    #
    #         Delta_{i \bar{j}} = < i \bar{j} | Delta(\vec{r}) | 0 \rangle
    #
    # where both \bar{j} is the timelike partner of j, and both i and j range
    # over only half of the basis. In that case, the state \bar{j} is 
    # not in memory. Combining the time-reversal needed to construct \bar{j}
    # with the intrinsic one in the definition of the pairing density/field, 
    # we obtain T^2 = -1, i.e. NO remaining time-reversal and a sign.
    #
    # Below, we take out the time-reversal if needed. The extra sign is not 
    # added here, but in the place where the actual action is constructed 
    # above.
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

    if('P' in field):
        # Add a timereversal operator on the left for pairing densities
        left  = left + 'T'
        if(timelike == 1):
            left = left + 'T'
    
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

def perm_parity(lst):
    '''
    Given a permutation of the digits 0..N in order as a list, 
    returns its parity (or sign): +1 for even parity; -1 for odd.
    '''
    parity = 1
    for i in range(0,len(lst)-1):
        if lst[i] != i:
            parity *= -1
            mn = min(range(i,len(lst)), key=lst.__getitem__)
            lst[i],lst[mn] = lst[mn],lst[i]
    return parity    
