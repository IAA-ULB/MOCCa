#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
#
# Module governing the calculation of mean-field potentials. 
#
#-------------------------------------------------------------------------------
# TODO
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

#def initfields(so):
#  """
#    Go over the needed densities and the functional terms and check whether
#    we have enough derivatives of the spwfs to calculate the fields.
#  """ 
#  for term in src_heph.heph_functional.Functional_terms:
#      (densities, cpl) = src_heph.heph_functional.ParseDensities(term)
#      # Count the number of derivatives needed in this term
#      totalder = 0
#      totallap = 0
#      for den in densities:
#          (der,lap,left,right, coupling, cross) = \
#                                                 ParseOperators(den,so.timelike)
#          totalder = totalder + der
#          totallap = totallap + lap
#      
#      # Now see that for all densities in this term, the minimum number
#      # of derivatives is the total one 
#      for i in range(len(densities)):
#          den = densities[i]
#          (der,lap,left,right,coupling,cross) = ParseOperators(den,so.timelike)
#          for j in range(len(Densities_needed)):
#              altden = Densities_needed[j]
#              (altder, altlap, altleft, altright, altcoupling, altcross)     \
#                                            = ParseOperators(altden,so.timelike)    
#              if(altleft == left and altright == right):
#                  # Set minimum derivatives
#                  deriv_needed[j].append((totallap, totalder))
#                  if(len(densities)>2):
#                    # For trilinear and quadrilinear terms, we will need 
#                    # more derivatives, as the laplacians can "uncouple"
#                    # for the calculation of the fields
#                    deriv_needed[j].append((0, totalder+totallap))

#  src_heph.heph_functional.PruneDeriv_needed()
      
def GenerateFields(so, oldso, ph_pp_decoupl, fam_active):
  """
   Generate a list of fields based on list of terms in the functional. 
   
   Attention!
    The output of this routine depends on whether or not PH_PP_DECOUPL is True
    or not. If True, ALL possible contributions to the fields are taken into
    account. If False, all contributions to fields associated with normal 
    densities by pairing densities are omitted.
   
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   Input: 
    so   :  symmetry options for the CURRENT executable being constructed
    oldso:  symmetry options for the executable that wrote the .wf file
    ph_pp_decoupl: Boolean. If True, drop all contributions to the normal
                   potentials that arise from density-dependent pairing 
                   terms.
    fam_active: Boolean. If True, declare fields as complex for the FAM code.
                If False, declare fields as real for the mean-field code.
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
   Output: 
   
    declaration:
    FIELDCALC  : 
    fieldprecon:
    fieldwrite :
    fieldread  :
    fieldclean :
    fieldINMk2 :
    fieldINMk4 :
   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  """
      
  global sumindices,tab
  
  import src_heph.fortran_templates.GenerateFields_templates as ts 

  #---------------------------------------------------------------------------
  # For every unique density encountered, we need to figure out the field
  # and the action of the field. 
  FIELDCALC           = ''
  FIELDCALC_perturbed = ''
  fieldprecon = ''
  declaration = ''

  fieldread       = ''
  fieldread_hdf5  = ''
  fieldwrite      = ''
  fieldwrite_hdf5 = ''

  fieldtransfo = ''

  fieldINMk2 = ''
  fieldINMk4 = ''
      
  fieldini     = ''
  fieldadd     = ''
  fieldmultiply= ''
 
  fieldinproduct = ''
  fieldclean = ''
  #---------------------------------------------------------------------------
  for den in src_heph.heph_functional.Densities_needed:
      #-----------------------------------------------------------------------
      #  Get the operator structure of the density correctly                    
      (der,lap,left, right, coupling,cross) = ParseOperators(den,so.timelike) 

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
        dic['POTREAD']  = '!'
      else:
        # Tantalus will read potentials correctly from file
        dic['UNDOREAD'] = '!'
        dic['POTREAD']  = ' '

      #-------------------------------------------------------------------------
      # Check if this particular density is truly a relevant one, i.e. does it
      # occur in at least one term of the functional? This is not guaranteed, 
      # since setting DENSITY_SPWF_SUMMATION generates additional intermediate
      # densities.
      found = False
      for term in src_heph.heph_functional.Functional_terms:
        (densities,coup) = src_heph.heph_functional.ParseDensities(term)
        for test in densities:
          (derj, lapj, leftj, rightj, coupj, crossj) = \
                                                ParseOperators(test,so.timelike)
          if(leftj == left and rightj == right):
            found = True

      if(not found):
        # If the density is not a relevant one, we don't define the associated
        # field and hence save a ton of CPU time.
        continue
      #-------------------------------------------------------------------------
      if('P' not in den): 
        Fields_needed.append(dic['FIELD'])
      else:
        Pairing_Fields_needed.append(dic['FIELD'])
      #-----------------------------------------------------------------------
      # Construct the left/right operators
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
          LeftOperator  = Combine(operatordic[l], LeftOperator)
          
      RightOperator = Identity
      for i in range(len(right)):
          r = right[(len(right)) -i -1] 
          RightOperator = Combine(operatordic[r], RightOperator)

      #-------------------------------------------------------------------------
      # Find the correct dimensions
      ndim = LeftOperator.dimension + RightOperator.dimension 
      # Subtract two dimensions for every summation, and only one for every 
      # vector product    
      ndim = ndim - 2*len(coupling) - len(cross)
  
      # But we have double-counted possible contractions of vector products,
      # couplings of length 3
      threedim=0
      threecoupl = []
      for c in coupling:
          if(len(c) == 3):
              threedim = threedim + 1
              threecoupl.append(c)
      ndim = ndim - threedim
      
      #-----------------------------------------------------------------------
      # Check all of the terms if they depend on the density
      fieldlist = {}
      fieldlist['0'] = []
      fieldlist['1'] = []        
      fieldlist['p'] = []
      fieldlist['n'] = []        
      cpcte     = ''

      #-----------------------------------------------------------------------
      for nterm, term in enumerate(src_heph.heph_functional.Functional_terms): 
          (densities, cpl)=src_heph.heph_functional.ParseDensities(term)
          iso_ind = src_heph.heph_functional.isospin_indices[nterm]
          extra_call = src_heph.heph_functional.extra_calls[nterm]
          
          #-------------------------------------------------------------------
          # Replace the densities in the list by the ones actually calculated
          newden = densities.copy()
          for i, d in enumerate(densities):
              (x,y,l2,r2,c2,cr2) = ParseOperators(d,so.timelike)
              for altden in src_heph.heph_functional.Densities_needed:
                  (altder, altlap, altleft, altright, altcoup, altcross) = \
                                            ParseOperators(altden,so.timelike)
                  if(altleft == l2 and r2 == altright and cr2 == altcross):
                      newden[i]    = y*'Lap_' +               \
                                     x*'Der_' +               \
                                     altden
          densities = newden
          #-------------------------------------------------------------------
          # Remove all the mentions of couplings inside the density if only
          # contractions are calculated.
          altterm = term
          for c in cpl: 
            for i in range(len(densities)):
              if(OrderOfDen(densities[i]) != OrderOfDen(densities[i], contract=False)):
               for s1 in sumindices:
                if(densities[i].count(s1) == 2):
                  altterm = altterm.replace(s1, "")
                  
            # WR 14/04/22: I'm unsure why I made things this 
            #              complicated. For future reference, 
            #              the below code fails when multiple
            #              couplings are present in one density
            #              Example : D_I_I_D_NmNm_NkNk term in
            #              N2LO functionals.
#                  for s2 in sumindices:
#                    tryout = densities[i].replace(s1, s2)
#                    nosum  = densities[i].replace(s1, '')
#                    altterm = altterm.replace(tryout, nosum)
          (rubbish, cpl) = src_heph.heph_functional.ParseDensities(altterm)

          #---------------------------------------------------------------------
          # We drop contributions to the normal potentials from pairing 
          # densities if PH_PP_DECOUPL is true
          mixedterm = False  
          for d in densities:
            mixedterm = mixedterm or ('P' in d)
          mixedterm = mixedterm and ('P' not in den)
          if(mixedterm and ph_pp_decoupl):
            continue
          #-------------------------------------------------------------------
          # Loop over the possible isospin components of this density
          if ('P' in den):
            isorange = ['n', 'p'] # Pairing densities are treated in 
                                  # proton-neutron formalism
          else:
            isorange = ['0', '1'] # Normal densities are treated in isospin
                                  # formalism
          for iso in isorange:
            # Check if the term of the EDF contains this isospin component 
            # of this particular density
            startind = 0
            for i, altden  in enumerate(densities):
                (altder, altlap, altleft, altright, altcoup, altcross) \
                                            = ParseOperators(altden,so.timelike)
                
                if(     (altleft == left) 
                    and (altright == right) 
                    and (iso == iso_ind[i])):
                  #  Add the term to the fieldlist for this density andd isospin
                  #  and additionally mentioning the number of external 
                  #  derivatives and laplacians
                  
                  dden  = src_heph.heph_functional.density_dependence[nterm]

                  if(dden == '1' or i != 0):
                    # Either
                    #  a) the term is not density-dependent; or
                    #  b) the functional derivative did not fall on the 
                    #     density dependence
                    # Proceed as standard
                    # - - - - - - - - - - -- - - - - - - - - - - - - - - - - -

                    # First we count how much indices are accounted for by the 
                    # other densities in the term
                    removed    = []
                    removedsum = 0
                    for j in range(len(densities)):
                      if i != j :
                        removed.append(densities[j])
                    for j in range(i):
                        removedsum = removedsum + OrderOfDen(densities[j])

                    cplct = 'coupl_constant(%2d)'%(nterm+1)
                    isoc  = iso_ind.copy()
                    isoc.pop(i)                                          
                                        
                    # When terms are trilinear or quadrilinear, every partial 
                    # integration generates MORE THAN ONE TERM in the 
                    # calculation of the field. Example:
                    #   fg Delta (h) => F_h ~ Delta(f) g + f Delta(g) 
                    #                                   + 2 nabla(f) nabla(g)  
                    nden = len(removed)
                    #    = Total number of densities the derivatives can fall on
                    indiv_lap = itertools.product(range(nden), repeat=2)
                    lap_combinations = itertools.product(indiv_lap, repeat=altlap)
                    #    = All possibilities to redistribute the laplacian operators
                 
                    for lcmb in lap_combinations:
                      # Need to regenerate this everytime, due to the structure of Python iterators
                      der_combinations = itertools.product(range(nden),repeat=altder)
                      #    = All possibilities to redistribute the individual derivatives
                      for dcmb in der_combinations:
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

                        if(dcmb.count(0) == 1 and dden != '1'):
                          if(len(dcmb) !=  1):
                            print ("Hephaestos cannot yet deal with derivatives upon derivatives of density dependent terms!")
                            exit(1)
                          # The derivative falls on the density dependence!
                          # nabla_m rho^alpha = alpha rho^(alpha-1) nabla_m rho
                          cc      = '(%s) * '%(dden) + cplct
                          dd      = '(' + dden + ' - 1)'
                          
                          a = removed.copy()
                          a.insert(0, removed[0])
                          b = isoc.copy()
                          b.insert(0,isoc[0])

                          newd = ()
                          for i in range(len(dcmb)):
                             newd = newd + (dcmb[i]+1,)
                          derc = newd
                          lapc = lcmb
                        
                        elif(dden != '1' and len(lcmb) != 0):
                          # There is (maybe) (part of ) a Laplacian falling on 
                          # the density dependent factor
                          
                          if(len(lcmb)>1):
                            print ("Repeated derivatives of density dependence.")
                            exit(1)

                          lc = lcmb[0]                          
                          if( lc.count(0) == 1):
                              # Single derivative
                              cc      = '(%s) * '%(dden) + cplct
                              dd      = '(' + dden + ' - 1)'
                              
                              a = removed.copy()
                              a.insert(0, removed[0])
                              b = isoc.copy()
                              b.insert(0,isoc[0])

                              newd = ()
                              for i in range(len(lc)):
                                 newd = newd + (lc[i]+1,)
                              derc = dcmb
                              lapc = (newd,)
                              
                          elif( lc.count(0) == 2):
                              # Complete Laplacian acting on rho^alpha
                              continue
                          else:
                            # No partial derivatives fall on the density dependent term
                            a   = removed
                            b   = isoc
                            derc= dcmb
                            lapc= lcmb
                            cc  = cplct
                            dd  = dden
                        else:
                          # No partial derivatives fall on the density dependent term
                          a   = removed
                          b   = isoc
                          derc= dcmb
                          lapc= lcmb
                          cc  = cplct
                          dd  = dden
                          
                        # We add a term to the expression for the field
                        fieldlist[iso].append([    a, # all the densities in the contribution to the field
                                                derc, # derivative combination
                                                lapc, # laplacian combination
                                                  cc, # Coupling constant
                                              newcpl, # Coupling of the indices 
                                                   b, # Isospin indices
                                                  dd, # Density dependence power
                                          extra_call # Extra call to a subroutine
                                                ])
                  else:
                    #-----------------------------------------------------------
                    # The term is density dependent AND the functional 
                    # derivative fell on the density in the density-dependence
                    cplct = 'coupl_constant(%d)'%(nterm+1)
                    altdden = '(' + dden + ' - 1)'                          
                    altcpl  = cplct + ' * (' +  dden + ')' 
                    isoc    = iso_ind.copy()
                    #-----------------------------------------------------------
                    # When terms are trilinear or quadrilinear, every partial 
                    # integration generates MORE THAN ONE TERM in the 
                    # calculation of the field. Example:
                    #   fg Delta (h) => F_h ~ Delta(f) g + f Delta(g) 
                    #                                   + 2 nabla(f) nabla(g)  
                    nden = len(densities)
                    #    = Total number of densities the derivatives can fall on
                    indiv_lap = itertools.product(range(nden), repeat=2)
                    lap_combinations = itertools.product(indiv_lap, repeat=altlap)
                    #    = All possibilities to redistribute the laplacian operators
                    for lcmb in lap_combinations:
                      der_combinations = itertools.combinations_with_replacement(range(nden),   altder)

                      for dcmb in der_combinations:
                        fieldlist[iso].append([densities, 
                                                    dcmb, 
                                                    lcmb,
                                                  altcpl,
                                                     cpl,
                                                    isoc,
                                                altdden, 
                                                extra_call])
                  
                startind = startind + OrderOfDen(altden) 
      # End of the costruction of all terms in the fields
      # Now we BUILD the code that calculates all these terms
      #-----------------------------------------------------------------------
      # Create the expression for the field
      dic['ALLOCIND']     = ''
      dic['ALLOCINDHDF5'] = ''
      dic['DECLIND']      = ''
      for k in range(OrderOfDen(den)):
          dic['ALLOCIND']     = dic['ALLOCIND']     + ',3'
          dic['ALLOCINDHDF5'] = dic['ALLOCINDHDF5'] + '*3' 
          dic['DECLIND']      = dic['DECLIND']      + ',:'
      if('P' in den):
        dic['ISOSIZE'] = 2
      else:
        dic['ISOSIZE'] = 4
        
      if(fam_active):
        declaration  = declaration + ts.field_decl_complex.substitute(dic)
      else:
        declaration  = declaration + ts.field_decl_real.substitute(dic)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # Explicit declaration of the history is no longer needed
      #declaration  = declaration + ts.fhist_decl.substitute(dic)
      
      # ... and the expression for reading/writing the fields from file
      #     NOTE THAT ONLY NEUTRON/PROTON FIELDS ARE WRITTEN/READ FROM FILE
      fieldread    = fieldread   + ts.field_read_a.substitute(dic)
      #fieldread    = fieldread   + ts.field_allo_b.substitute(dic)
      fieldread    = fieldread   + ts.field_read_b.substitute(dic)
      fieldread    = fieldread   + ts.field_read_c.substitute(dic)
      fieldread    = fieldread   + ts.field_read_d.substitute(dic)
      fieldread    = fieldread   + ts.field_read_e.substitute(dic)
 
      args = list(itertools.product(range(3), repeat=OrderOfDen(den)))
      for arg in args:
           # get the indices of the field correct
           dic['IND']     = ''
           for k in arg:
              dic['IND'] = dic['IND'] + ',%d'%(k+1)

           fieldread = fieldread + ts.field_transfo.substitute(dic)
      if('P' not in den):
        # Only recombine normal fields
        fieldread = fieldread + ts.field_transfo_recomb.substitute(dic)
      fieldread = fieldread + ts.field_read_f.substitute(dic)
      
      # HDF5 option
      fieldread_hdf5 = fieldread_hdf5 + ts.field_read_hdf5_a.substitute(dic)
      fieldread_hdf5 = fieldread_hdf5 + ts.field_read_hdf5_b.substitute(dic)
      fieldread_hdf5 = fieldread_hdf5 + ts.field_read_hdf5_c.substitute(dic)
      fieldread_hdf5 = fieldread_hdf5 + ts.field_read_hdf5_d.substitute(dic)
      fieldread_hdf5 = fieldread_hdf5 + ts.field_read_hdf5_e.substitute(dic)
 
      args = list(itertools.product(range(3), repeat=OrderOfDen(den)))
      for arg in args:   
           # get the indices of the field correct
           dic['IND']     = ''
           for k in arg:
              dic['IND'] = dic['IND'] + ',%d'%(k+1)
              
           fieldread_hdf5 = fieldread_hdf5 + ts.field_transfo_hdf5.substitute(dic)
      if('P' not in den):
        # Only recombine normal fields
        fieldread_hdf5 = fieldread_hdf5 + ts.field_transfo_recomb_hdf5.substitute(dic)
      fieldread_hdf5 = fieldread_hdf5 + ts.field_read_hdf5_f.substitute(dic)
        
      fieldwrite   = fieldwrite  + ts.field_write_a.substitute(dic)
      fieldwrite   = fieldwrite  + ts.field_write_b.substitute(dic)

      fieldwrite_hdf5   = fieldwrite_hdf5  + ts.field_write_hdf5.substitute(dic)

      fieldini     = fieldini  + ts.field_allo.substitute(dic)

      FIELDCALC    = FIELDCALC + ts.field_line.substitute(dic)
      FIELDCALC    = FIELDCALC + ts.field_calc_a_start.substitute(dic)
      
      FIELDCALC_perturbed = FIELDCALC_perturbed + ts.field_line.substitute(dic)
      FIELDCALC_perturbed = FIELDCALC_perturbed + ts.field_calc_a_start.substitute(dic)
      
      FIELDCALC           = FIELDCALC           + ts.field_condition_start.substitute(dic)
      FIELDCALC_perturbed = FIELDCALC_perturbed + ts.field_condition_start.substitute(dic)
      
      fieldinproduct = fieldinproduct + ts.field_inproduct.substitute(dic)

      #------------------------------------------------------------------------
      #
      fieldINM = ""
      #-------------------------------------------------------------------------
      # Code generation for the calculation of the fields. 
      #
      # Isospin loop:
      for iso in ['0','1', 'n', 'p']:        
      
        if(len(fieldlist[iso]) == 0):
          continue
      
        dic['ISO']    = str(iso)
        dic['ISOIND'] = Isospinindices(iso)
        
        FIELDCALC            = FIELDCALC           + ts.field_calc_iso_start.substitute(dic)
        FIELDCALC_perturbed  = FIELDCALC_perturbed + ts.field_calc_iso_start.substitute(dic)
        for fieldterm in fieldlist[iso]:
           dic['CPLCTE']   =  fieldterm[3]

           #-----------------------------------------------------------------
           # We rearrange things in this separate routine
           NumberOfIndices, nthree, densities, newcpl, pisigns = \
            Adaptdensities(fieldterm[0],fieldterm[4],fieldterm[1],fieldterm[2])

           # We construct arguments for all latin indices
           args = list(itertools.product(range(3), repeat=NumberOfIndices))
           vec_args = list(itertools.product(range(6), repeat=nthree))

           if(nthree == 0):
              true_args = args
           elif(NumberOfIndices == 0 ):
              true_args = vec_args
           else:
              true_args = []
              for a in args:
                for va in vec_args:
                  true_args.append(a+va)

           # This quantity counts how many index combinations of this term
           # give non-zero contributions to the potentials in homogeneous,
           # unpolarised infinite nuclear matter.
           INM_term_count = 0
           for arg in true_args:
             # get the indices of the field (i.e. the lhs above) correct
             dic['IND']     = ''
             globalsign     = +1

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
                    globalsign = globalsign * -1
                  temp = (mu,abs(nuka[0]), abs(nuka[1]))
                  dic['IND'] = dic['IND'] + ',%d'%(temp[c.index(k)]+1)

             #----------------------------------------------------------------
             # Before we start doing complicated stuff for the index coupling
             # of fields, we first look into the calculation of potentials
             # in homogeneous unpolarised infinite nuclear matter.
             #
             # We select on terms for which
             #  - no Pauli sigma matrices are involved
             #  - either two or four derivatives in the density
             #  - no external derivatives of densities are involved
             #  - indices couple to a scalar, i.e. only terms that
             #    contain the diagonal part of the tensor densities
             #      D_NNN...._NNN_{aaaaaa..., aaaaa...}
             #
             # If all of this is true, we know that this particular term
             # with this particular combination of indices in non-zero in INM.
             # We could add a line in the FORTRAN code, but this would result
             # in a large amount of repeated lines.
             #
             # Note: I'm not entirely sure (yet) whether this selection procedure
             # is entirely general; I've double checked NLO and standard N2LO
             # EDFS but more general things might not work
             if(   LeftOperator.derorder + RightOperator.derorder == 2 \
                 or LeftOperator.derorder + RightOperator.derorder == 4 ):
                 # Pauli spin matrices?
                 nosigma=True
                 if('S' in left or 'S' in right):
                     nosigma=False
                  # Any external derivatives?
                 noderivatives=True
                 if(len(fieldterm[1])>0 or len(fieldterm[2])>0):
                   noderivatives = False
                 for d in fieldterm[0]:
                   if('Der' in d or 'Lap' in d):
                       noderivatives=False

                 # - - - - - - - - - - - - - - - - - - - - - - - - - - -
                 # Poor mans reverse engineering of the indices to
                 # determine if they couple to a scalar in position space
                 if(len(dic['IND']) > 0):
                    scalar = True
                    inds = []
                    for k in dic['IND']:
                     try:
                      inds.append(int(k))
                     except:
                      continue
                    for k in range(1,max(inds)+1):
                     if(inds.count(k)%2 != 0):
                      scalar = False
                 else:
                    scalar = True

                 if(nosigma and noderivatives and scalar):
                   INM_term_count = INM_term_count + 1
             #----------------------------------------------------------------
              
             dic['DENSITY']  = ''
             dic['EXPR1']    = ''
             dic['EXPR_PERT'] = ''
             
             # TODO: refactor this into a function  
             # Build the expression for the traditional mean-field densities
             lastorder = OrderOfDen(den)
             FIELDCALC = FIELDCALC + ts.field_calc_b_start.substitute(dic)
             sign      = globalsign
             for i,d in enumerate(densities):
                dic['DENSITY'] = d
                #-----------------------------------------------------------
                # Get the indices of the density on the rhs.
                indices = ()
                for k in range(lastorder,lastorder+OrderOfDen(d)):
                  for c in newcpl: #fieldterm[4]:
                    if k in c:   
                      if(len(c) == 2):    
                        mu      = arg[newcpl.index(c)]  #fieldterm[4].index(c)]
                        indices = indices + (mu,)
                      elif(len(c) == 3):
                        # Integer division
                        mu   = int(arg[newcpl.index(c)]/2) 
                        # Which term of two?
                        t = arg[newcpl.index(c)] - 2*mu
                        if(t == 1):
                          sign = sign * -1

                        nuka = Rot_ind(mu)[t] 
                              
                        temp = (mu,abs(nuka[0]), abs(nuka[1]))
                        indices = indices + (temp[c.index(k)],)                                 

                #-----------------------------------------------------------
                # The first indices are necessarily external derivatives
                dercount = d.count('Der')
                if(dercount > 0):
                    derind  = Storage_Mapping(indices[:dercount])
                    indices = (derind,) + indices[dercount:]
                
                dic['DENIND']      = ''
                for l in indices:
                    dic['DENIND'] = dic['DENIND'] + ',%d'%int(l+1)
                dic['ISOALT']= Isospinindices(fieldterm[5][i])
                if(fieldterm[+6] != '1' and i == 0): 
                  dic['DD'] = fieldterm[+6]
                  # Expression with a call to 'pow'
                  dic['EXPR1'] = dic['EXPR1'] + ts.field_calc_DD.substitute(dic)
                else:
                  # Ordinary expression
                  dic['EXPR1'] = dic['EXPR1'] + ts.field_calc_den.substitute(dic)

                # Increment the starting point of indices
                lastorder = lastorder + OrderOfDen(dic['DENSITY'])

             #------------------------------------------------------------------
             # Put an extra sign for every partial integration of a nabla
             if( len(fieldterm[1])%2 != 0):
                 localsign = sign * (-1)
             else:
                 localsign = sign

             if(localsign > 0):
               dic['SIGN']    =  '+'
             else:
               dic['SIGN']    =  '-'

             if(fieldterm[-1] == ''):
               dic['EXTRA'] = ''
             else:
               dic['EXTRA'] = '*(%s)'%fieldterm[-1]

             FIELDCALC = FIELDCALC + ts.field_calc_full.substitute(dic)
             FIELDCALC = FIELDCALC[:-4] + '\n \n'
             
             # Now do it again, repeatedly, but taking one of the densities
             # from the perturbed density vector at any one time
             for j in range(len(densities)):
               lastorder = OrderOfDen(den)
               FIELDCALC_perturbed = FIELDCALC_perturbed + ts.field_calc_b_start.substitute(dic)
               dic['EXPR_PERT']    = ''
               sign = globalsign
               # Pick the J-th density to be a perturbation; otherwise do the 
               # same thing as above...
               for i,d in enumerate(densities):
                  dic['DENSITY'] = d
                  #-----------------------------------------------------------
                  # Get the indices of the density on the rhs.
                  indices = ()
                  for k in range(lastorder,lastorder+OrderOfDen(d)):
                    for c in newcpl: #fieldterm[4]:
                      if k in c:   
                        if(len(c) == 2):    
                          mu      = arg[newcpl.index(c)]  #fieldterm[4].index(c)]
                          indices = indices + (mu,)
                        elif(len(c) == 3):
                          # Integer division
                          mu   = int(arg[newcpl.index(c)]/2) 
                          # Which term of two?
                          t = arg[newcpl.index(c)] - 2*mu
                          if(t == 1):
                            sign = sign * -1

                          nuka = Rot_ind(mu)[t] 
                                
                          temp = (mu,abs(nuka[0]), abs(nuka[1]))
                          indices = indices + (temp[c.index(k)],)                                 

                  #-----------------------------------------------------------
                  # The first indices are necessarily external derivatives
                  dercount = d.count('Der')
                  if(dercount > 0):
                      derind  = Storage_Mapping(indices[:dercount])
                      indices = (derind,) + indices[dercount:]
                  dic['DENIND']      = ''
                  for l in indices:
                      dic['DENIND'] = dic['DENIND'] + ',%d'%int(l+1)
                  dic['ISOALT']= Isospinindices(fieldterm[5][i])
                  if(fieldterm[+6] != '1' and i == 0):
                    dic['DD']      = fieldterm[+6]
                    dic['DD_pert'] = fieldterm[+6]  #+ '-1'
                    # Expression with a call to 'pow'
                    if(i != j):
                      dic['EXPR_PERT'] = dic['EXPR_PERT'] + ts.field_calc_DD.substitute(dic)
                    else:
                      dic['EXPR_PERT'] = dic['EXPR_PERT'] + ts.field_calc_DD_pert.substitute(dic)
                  else:
                    # Ordinary expression
                    if (i != j):
                      dic['EXPR_PERT'] = dic['EXPR_PERT'] + ts.field_calc_den.substitute(dic)
                    else:
                      if('P' in densities[j]):
                        # Extraordinary dirty trick: if we pick a pairing density to vary, we add an extra factor
                        #   0.5 in order to offset the double-counting.
                        #
                        # Note: this will ONLY work if pairing terms are of the simple form (pairing density) * (pairing density)^*
                        #       i.e. one cannot have EDF terms with different pairing densities!
                        dic['EXPR_PERT'] = dic['EXPR_PERT'] + ' * 0.5d0' + ts.field_calc_den_pert.substitute(dic)
                        print (densities[j], j,i, dic['EXPR_PERT'])
                      else:
                        dic['EXPR_PERT'] = dic['EXPR_PERT'] + ts.field_calc_den_pert.substitute(dic)

                  # Increment the starting point of indices
                  lastorder = lastorder + OrderOfDen(dic['DENSITY'])

               #------------------------------------------------------------------
               # Put an extra sign for every partial integration of a nabla
               if( len(fieldterm[1])%2 != 0):
                  localsign = sign * (-1)
               else:
                  localsign = sign

               if(localsign > 0):
                dic['SIGN']    =  '+'
               else:
                dic['SIGN']    =  '-'

               if(fieldterm[-1] == ''):
                dic['EXTRA'] = ''
               else:
                dic['EXTRA'] = '*(%s)'%fieldterm[-1]

               FIELDCALC_perturbed = FIELDCALC_perturbed + ts.field_calc_pert.substitute(dic)
               FIELDCALC_perturbed = FIELDCALC_perturbed[:-4] + '\n \n'

           #------------------------------------------------------------------
           # This term contributes to potentials in homogeneous unpolarised
           # infinite matter, so we count how exactly it does this.
           if(INM_term_count > 0):
            dic['EXPR1']  = ''
            dic['SIGN']   = '+'

            for i,d in enumerate(densities):
              # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
              # For now, this implementation is limited to densities of the form
              #       D^NNNNN.....,NNNNN.....
              # with 2*l nabla operators. They translate to INM as
              #
              #   isospin 0 : [(3+2*l) pi^2 ]^{-1}  (kfn^{3+2*l} + kfp^{3+2*l})
              #   isospin 1 : [(3+2*l) pi^2 ]^{-1}  (kfn^{3+2*l} + kfp^{3+2*l})
              #
              power         = 3 + OrderOfDen(d)
              # Isospin indices
              if(fieldterm[5][i] == '0'):
                isosign = '+'
              else:
                isosign = '-'
              # There is additional degeneracies to account for:
              #    consider the general equation
              #
              #      F =  C * D_1 * D_2 * ...
              #    where F is a field and the Ds are possibly different densities.
              #
              #    - If F has indices, there is a factor three degeneracy since only
              #      its scalar part contributes to INM.
              #    - For each of the densities that has indices, there is an additional
              #      factor of three degeneracy because only their scalar part will be
              #      non-zero.
              #
              degeneracy = 1
              if(OrderOfDen(d) > 0):
                degeneracy = degeneracy * 3
              if(OrderOfDen(den) > 0):
                degeneracy = degeneracy * 3
              # We multiply by the number of valid index couplings we counted (INM_term_count)
              temp_string = " %2d/ (%2d * pi**2) * ( kfn**(%2d) %s kfp**(%2d) )"%(INM_term_count, degeneracy*power, power, isosign, power)
              # ... and we don't forget to account for density-dependencies
              if(fieldterm[6] != '1' and i == 0):
                temp_string = "pow( %s , %s)"%(temp_string, fieldterm[6])
              dic['EXPR1']  = dic['EXPR1'] + "*" + temp_string
            fieldINM = fieldINM + ts.field_calc_INM.substitute(dic)
        fieldINM = fieldINM + '\n'

      # Add the recombination statements from isospin representation to 
      # proton-neutron representation, but only for normal densities
      if('P' not in den):
        FIELDCALC           = FIELDCALC           + ts.field_recombination.substitute(dic)
        FIELDCALC_perturbed = FIELDCALC_perturbed + ts.field_recombination.substitute(dic)

      FIELDCALC    = FIELDCALC + ts.field_condition_end.substitute(dic)
      FIELDCALC    = FIELDCALC + ts.field_line.substitute(dic) + '\n'

      FIELDCALC_perturbed    = FIELDCALC_perturbed + ts.field_condition_end.substitute(dic)
      FIELDCALC_perturbed    = FIELDCALC_perturbed + ts.field_line.substitute(dic) + '\n'

      # Add some lines for the multiplication and addition of potentialvectors!
      fieldadd      = fieldadd      + '\n' +  ts.Add.substitute(dic)
      fieldmultiply = fieldmultiply + '\n' +  ts.Multiply.substitute(dic)

      if(LeftOperator.derorder + RightOperator.derorder == 2):
        fieldINMk2 = fieldINMk2 + fieldINM
      elif(LeftOperator.derorder + RightOperator.derorder == 4):
        fieldINMk4 = fieldINMk4 + fieldINM

      #-------------------------------------------------------------------------
      # Code generation for the preconditioning of the fields
      #
      # Current strategy: precondition ALL potentials, irrespective of their 
      #                   content in terms of gradients. This should remove
      #                   all potential for surprises and makes sure that 
      #                   the parameter preconfactor influences ALL potentials.
      #                   The only exception (for now) are the pairing potentials.
      # 
      if('P' not in den):
        fieldprecon  = fieldprecon + ts.field_precon_start.substitute(dic)

        # TODO: define a dedicated function for this next particular piece of 
        #       code. It is highly complex, yet identical to the one in the 
        #       src_heph/heph_densities module. It should be put in one spot
        #       and abstracted.

        args = list(itertools.product(range(3), repeat=ndim))
        # Outer loop over all possible indices

        for arg in args:
          # We have the uncontracted indices. Now construct the combinations of
          # indices, including contracted ones, that correspond to this. 
          uncontracted = []
          if(len(coupling) + len(cross) == 0):
              # Nothing to do if no couplings needed
              uncontracted = [arg]
          else:
              #-------------------------------------------------------------------
              # These are all of the combinations needed for the summation indices
              cont  = itertools.product(range(3), repeat=len(coupling))
              # All of the possibilities for the vector products
              crossind = []
              for i in range(len(cross)):
                  crossind = crossind  + (Rot_ind(arg[len(coupling) + i]))    

              # Combine all possibilities
              if(len(crossind) != 0):
                  # all the combinations , including scalar and vector products
                  fullcont = []
                  for s in cont: 
                      for r in crossind:            
                          fullcont.append(s + r)
              elif (threedim != 0 and len(crossind) == 0):
                  fullcont = []
                  for s in cont: 
                      # Modify the number of possibilities for couplings between vector
                      # products and scalar products. 
                      threeopt = itertools.product(range(2), repeat=threedim)
                      for r in threeopt:            
                          fullcont.append(s + r)
              elif (threedim != 0 and len(crossind) != 0):
                  fullcont = []
                  for s in cont: 
                     for r in crossind: 
                      threeopt = itertools.product(range(2), repeat=threedim)
                      for t in threeopt:            
                          fullcont.append(s + r + t)
              else:
                  # No vector indices, and no scalar-vector
                  fullcont = cont
                
              #-------------------------------------------------------------------
              # Note that now fullcont contains all of the terms needed for the
              # particular argument of the left-hand side. 
              #
              #  The ordering of the indices is:
              #    
              #  ( mu, nu, ...., xsi , mx, nx, ....,zx  ,  ex, ey, ....., ez )  
              #   < scalar indices >  < vector indices >  < contracted vectors)
              #
              #  corresponding to things of the form
              #
              #  coupling, (0,1)        cross (0,1)         coupling (0,1,2)
              #
              #  meaning 
              #
              #  the value of the   | the values of the  |  whether it is the 
              #  indices in the     | vector indices     |  first term or the 
              #  summation          |                    |  second in the vector
              #                     |                    |  product
              #-------------------------------------------------------------------
              for c in fullcont:
                  p  = ()   
                  ii = 0
                  for i in range(LeftOperator.dimension + RightOperator.dimension):
                      found = False                   
                      for combination in coupling:
                          if(i in combination): 
                              if(len(combination) == 2):
                                  p = p + (c[coupling.index(combination)],)
                                  found = True
                              elif(len(combination) == 3):    
                                  found = True
                                  if( i == combination[0] ):
                                      p = p + (c[coupling.index(combination)],)
                                  elif( i == combination[1]):
                                      rot = Rot_ind(c[coupling.index(combination)])
                                      o = threecoupl.index(combination)
                                      p = p + (rot[c[-1 -o]][0],)
                                  elif( i == combination[2]):
                                      rot = Rot_ind(c[coupling.index(combination)])
                                      o   = threecoupl.index(combination)
                                      p = p + (rot[c[-1 -o]][1],)
                      for combination in cross:
                          if(i==combination[0]): 
                                  p = p + (c[cross.index(combination) + len(coupling)],)
                                  found = True
                          if(i==combination[1]): 
                                  p = p + (c[cross.index(combination) + len(coupling) +1 ],)
                                  found = True    
                      if(not found): 
                              p = p + (arg[ii],)
                              ii = ii +1
                  uncontracted.append(p)

          IND = ''
          for mu in arg: 
              IND = IND + ',' + str(int(abs(mu)+1)) # Python indexes 0:N-1

          dic['IND'] = IND

          #--------------------------------------------------------------------
          # Add the actual call(s) to precondition this potential to the code.
          # First we coun the number of external derivatives that occurs in the 
          # field, such that we know how often to apply the preconditioner.
          extder = 0
          for iso in ['0','1', 'n', 'p']:        
            if(len(fieldlist[iso]) == 0):
              continue
            for fieldterm in fieldlist[iso]:
              nder = 0 ; nlap = 0
              # Number of external derivatives in the "density", i.e. not incurred by 
              # partial integration
              for tempden in fieldterm[0]:
                  nder = nder + tempden.count('Der')
                  nlap = nlap + tempden.count('Lap')
              # extra derivatives incurred by partial integration
              nder = nder + len(fieldterm[1])
              nlap = nlap + len(fieldterm[2])
              # Every laplacian counts for two derivatives of course!
              extder = max(extder, nder + 2 * nlap)
          
          if(extder != 0):
            fieldprecon  = fieldprecon + ts.field_precon_update.substitute(dic)

            for true_arg in uncontracted: 
              # We do the whole loop but only use the values for the final set of
              # indices: all terms in a given contraction should have the same 
              # behaviour under symmetry.

              # The ugly tuple(np.abs( construction is simply because abs doesn't 
              # accept tuples as arguments, for whatever reasons.
              larg = tuple(np.abs(true_arg[:LeftOperator.dimension]))
              rarg = tuple(np.abs(true_arg[LeftOperator.dimension:]))

            # The ugly tuple(np.abs( construction is simply because abs doesn't 
            # accept tuples as arguments, for whatever reasons.
            larg = tuple(np.abs(true_arg[:LeftOperator.dimension]))
            rarg = tuple(np.abs(true_arg[LeftOperator.dimension:]))

            (px,py,pz)   = AxisReflection(LeftOperator, RightOperator,larg,rarg,so,'P' in den)
            dic['PX'] = str(px)
            dic['PY'] = str(py)
            dic['PZ'] = str(pz)
            # Every two external derivatives get one preconditioning run.
            dic['DIVISOR'] = (extder//2)**2 # ansatz for the appropriate preconfactor
            fieldprecon  = fieldprecon + + (extder//2) * ts.field_precon_call.substitute(dic)
            fieldprecon  = fieldprecon + ts.field_precon_add.substitute(dic)
  #-----------------------------------------------------------------------------

  return(declaration, fieldini, FIELDCALC, FIELDCALC_perturbed, fieldprecon, fieldwrite,fieldwrite_hdf5, \
   fieldread,fieldread_hdf5,fieldadd,fieldmultiply,fieldinproduct,fieldINMk2,fieldINMk4)

def Adaptdensities( dens, cpl, dcmb, lcmb):
  """
  
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    Input  : 
      dens      :
      cpl       : 
      lcmb      : 
      dcmb      :
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    Output : 
      NumberOfIndices :
      nthree          :
      densities       : 
      newcpl          : 
      pisigns         : 
  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  """
  
  densities = []
  newcpl    = []
  pisigns   = []
  
  NumberOfIndices = 0
  nthree = 0
  
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # We have here all ingredients to construct one contribution to one field 
  #
  # A statement in the code will have the following form formally
  #
  #   F/G_{mu nu .... kappa } =  F/G_{mu nu .... kappa }
  #    +/- C sum_{a,b,c,....}[\sum_{l=1}^{3} D_{l, mu nu ... kappa a b c .... }]
  #  
  #  where 
  #  * we have dropped the isospin index 
  #  * C is some coupling constant
  #  * l is the index of the density in the terms. We permit quadrilinear terms,
  #    so it is maximum three
  #  * D_l is some density; which includes the possibility of a DERIVATIVE of a
  #                         density
  #  * mu,nu,kappa : greek indices are indices of the FIELD being calculated
  #  * a,b;c; ...  : latin indices are indices of the r.h.s., to be summed over
  #                  for a given set of greek indices
  #
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The first order of business is to determine the number of latin indices!
  #
  # For this, we first check the coupling array
  for c in cpl:
    if(len(c) == 2):
      NumberOfIndices = NumberOfIndices + 1
    elif(len(c) == 3):
      nthree          = nthree +1  
  
  # But this is not sufficient, as Laplacian operators that need to be 
  # integrated by parts can "uncouple", i.e. act on different densities for
  # terms that are trilinear or quadrilinear.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Note: this kind of construction NEVER appears for bilinear terms
  #       and can always be eliminated by partial integration for trilinear
  #       terms. For quadrilinear terms however, this cannot be avoided.....
  uncoupled_laplacians = 0            
  for m, d in enumerate(dens):
   for lc in lcmb:
    nl = lc.count(m)
    if(nl%2 != 0):
      uncoupled_laplacians = uncoupled_laplacians + 1     

  uncoupled_laplacians = uncoupled_laplacians//2
  NumberOfIndices = NumberOfIndices + uncoupled_laplacians
  #- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Next, we construct the full expressions for the densities with the 
  # integration by parts of all kinds of derivatives
  for m,d in enumerate(dens):
    denstring  = ''
    if(len(lcmb) != 0):
      for lc in lcmb:
        nl = lc.count(m)
        
        if(nl == 2):
          # this laplacian remains coupled
          denstring = denstring + 'Lap_'
        elif(nl == 0):
          # No extra derivatives
          pass
        else:
          # this laplacian is uncoupled
          denstring = denstring + 'Der_'
           
    if(len(dcmb) != 0):
       nder = dcmb.count(m)
       denstring = denstring +  nder * 'Der_'
      
    # and we add the name of the density itself
    denstring = denstring + d 
   
    # Note that it is possible that the above code has mixed up the ordering 
    # of the "Lap_" and "Der_" expressions. We reorder them for surety
    dercount = denstring.count('Der')
    lapcount = denstring.count('Lap')
    
    denstring = denstring.replace('Der_', '')
    denstring = denstring.replace('Lap_', '')
    denstring = lapcount * 'Lap_' + dercount * 'Der_' + denstring
   
    # Adding the constructed density to the array
    densities.append(denstring)
    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Next, we modify the couplings of the indices as 
  #  (a) modify existing couplings to reflect the insertion of the derivatives
  #      from uncoupled laplacians
  for c in cpl:
    newc = ()
    for k in c:
      ind = 0
      inc = 0
      for m, d in enumerate(dens):
        km = 0
        for lc in lcmb:
          km = km + lc.count(m)%2
        if(k >= ind):
          inc = inc + km - 2*(km//2)
        
        ind = ind + OrderOfDen(d)
      newc = newc + (k+inc,)
    newcpl.append(newc)

  #  (b) insert the couplings for the uncoupled laplacians
  for lc in lcmb:
    newc = ()
    ind  = 0 
    for m, d in enumerate(densities):
      if(lc.count(m)%2 == 1):        
        newc = newc + (ind,)
      
      ind = ind + OrderOfDen(d)
    newcpl.append(newc)
  
  return NumberOfIndices, nthree, densities, newcpl, pisigns


def GenerateAction(field, symmetrize, so):
    """
     Generate the action of a field. 
    
        field       : name of the field, for example F_I_I
    
        symmetrize  : whether or not to generate a term for a symmetrised field  
                      (  0) generate action "as is"  
                      (+-1) generate one out of two terms of the action
    
        so          : a set of symmetry options 
    """
    
    import src_heph.fortran_templates.GenerateAction_templates as ts
    
    #---------------------------------------------------------------------------
    WFNames   = ['psi', 'dpsi', 'ddpsi', 'dddpsi', 'ddddpsi']
    Direction = ["X", 'Y', 'Z']
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
    expression = ts.action_comment.substitute(dic)
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
        expression = expression + ts.temp_ini
        
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
            expression = expression + ts.position_loop
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
                    expression = expression + ts.action.substitute(dic)
            expression = expression + ts.position_end
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
                    expression = expression + ts.lap.substitute(dic)
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
                    expression = expression + ts.sym.substitute(dic) 
                expression = expression + ts.derive.substitute(dic)
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
              expression = expression + ts.action_final.substitute(dic)
            else:
              expression = expression + ts.action_final_pairing.substitute(dic)
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
