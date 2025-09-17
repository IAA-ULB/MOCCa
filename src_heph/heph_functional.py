#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Module governing the treatment of the functional in Hephaestos, for writing to
# Tantalus source files. 
#-------------------------------------------------------------------------------

import itertools
import numpy as np
import sys

from src_heph.heph_densities import Densities_needed, tab, sumindices, derstring
from src_heph.heph_densities import lapstring, OrderOfDen, ParseOperators
from src_heph.heph_densities import crossindices, Storage_Mapping, Multiplicity
from src_heph.heph_densities import deriv_needed, Isospinindices
from src_heph.heph_linechecker import *
from src_heph.heph_fields      import *

#-------------------------------------------------------------------------------
# Name of the functional_file
func_name = ''

#-------------------------------------------------------------------------------
# Array containing the expressions of all the functional terms. 
Functional_terms      = []
# Array containing all the coupling constants
coupling_constants    = []
# Array containing the isospin indices of the densities in the terms
isospin_indices       = []
# array containing the density dependence of the first density of the term
density_dependence    = []
# Array containing all parameters are used to compute coupling constants
paramparameters       = []
# ... and a corresponding array that tracks their types 
# Currently accepted: Real    =  'R '
#                     Integer =  'I ' 
paramtypes            = []
# Array indicating the indices of a term within a subgroup of terms with 
# identical structure. I.e. "x" in this list means this the (x+1)-th term of
# this structure in the Functional_terms list.
term_grouping         = []
term_number           = {}
# The additional function calls that need to be made for the calculation of
# coupling constants
extra_calls           = []
vmicro_found          = False # Whether or not this functional file will 
                              # place calls to vmicro
################################################################################
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Switch determining what order of derivatives is needed to be computed.
#
# Note that this is not directly i => ith order.
#  derivative_order <=> derivatives of spwfs calculated
#            1          1st order + trace of 2nd order (laplacian)
#            2          2nd order + trace of 3rd order
#            3          3rd order derivatives
#derivative_order = 1
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
#-------------------------------------------------------------------------------
# Functionality removed on 22/12/2020.
#   It is wrong and results in unexpected behavior.
#-------------------------------------------------------------------------------
#assume_locality = 0

def initfunctional(fname, so, density_spwf_summation, fam_active):
    """
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      Initialize everything relevant about this module
      
      1) read the functional from the fname file
      2) remove time-odd terms if necessary
      3) generate all the densities required to calculate this functional
      4) figure out what order of derivatives is going to be needed for this EDF
      5) Print some output to STDOUT
      
      Input: 
        fname                  : filename containing the functional description 
        so                     : set of symmetry-options, determining whether 
                                 time-odd terms get kept or not.
        density_spwf_summation : logical determining whether derivatives of 
                                 densities get calculated through summation 
                                 over spwfs or derivative calls.
      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    """
    global Functional_terms, Densities_needed, derivative_order, func_name
    global extra_calls, vmicro_found, pairing_action_derorder
    
    #---------------------------------------------------------------------------
    # Read the functional from a given file
    description = ReadFunctional(fname)
    #---------------------------------------------------------------------------
    # Check if time-reversal (or time-parity) is conserved
    if(so.timelike):
      RemoveTimeOddTerms(fam_active)
      
    #---------------------------------------------------------------------------
    # Reorder terms, such that everything which is "grouped" by structure 
    # is put together
    regroup_terms()

    #---------------------------------------------------------------------------
    # Set the name of the functional file, without the directory structure
    func_name = '"%s"'%fname.split('/')[-1].upper()

    #---------------------------------------------------------------------------
    # Generating the list of all densities and the total number of derivatives
    tempden     = []
    minder      = []
    temp_intermediate = []
    for term in Functional_terms:
       (densities,coup) = ParseDensities(term)
       totalder = 0
       totallap = 0
       for den in densities:
           (der,lap,left,right, coupling, cross) = \
                                                 ParseOperators(den,so.timelike)
           # Count the total number of derivatives present in this term
           totalder = totalder + der
           totallap = totallap + lap

       # add density and set minimum number of derivatives for this term
       for den in densities:
           tempden.append(den)
           temp_intermediate.append(False) # these are "true" densities
           if(len(densities)<=2):
             minder.append((totallap, totalder))
           else:
             # For trilinear and quadrilinear terms, we will need 
             # more derivatives, as the laplacians can "uncouple"
             # for the calculation of the fields
             minder.append((0, totalder+2*totallap))

    #---------------------------------------------------------------------------
    # If density_spwf_summation is true, we have to add a bunch of densities
    # to the list of needed ones.
    # 
    # We use the identities:
    # 
    #       N D^L,R  = D^NL, R  + D^L, NR 
    #       NN D^L,R = D^NNL, R + D^L, NNR + D^NL,NR 
    # 
    # and similar for C-like objects.
    #---------------------------------------------------------------------------
    if(density_spwf_summation):
      for j,den in enumerate(tempden):
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        # Note: because we add to the end of this list while traversing it, we
        #       automatically do the entire process recursively.
        # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        (der, lap, left, right, coup, cross) = \
                        ParseOperators(den,so.timelike)
        if(minder[j][1] > 0):
          lleft = 'N' + left.replace('I','')
          tempden.append(ReconstructDensity(der-1, lap, lleft, right))
          temp_intermediate.append(True) # these are intermediate objects
          minder.append((0,0))

          rright = 'N' + right.replace('I','')
          tempden.append(ReconstructDensity(der-1, lap,      left, rright))
          minder.append((0,0))
          temp_intermediate.append(True) # these are intermediate objects

        if(minder[j][0] != 0 or minder[j][1] == 2):
          # Add in all densities needed for the Laplacian ....
          lleft  =  'NN' + left.replace('I','')
          tempden.append(ReconstructDensity(der, lap-1, lleft, right))
          temp_intermediate.append(True) # these are intermediate objects
          minder.append((0,0))

          rright =  'NN' + right.replace('I','')
          tempden.append(ReconstructDensity(der, lap-1, left, rright))
          temp_intermediate.append(True) # these are intermediate objects
          minder.append((0,0))

          lleft  =  'N' + left.replace('I','')
          rright =  'N' + right.replace('I','')
          tempden.append(ReconstructDensity(der, lap-1, lleft, rright))
          temp_intermediate.append(True) # these are intermediate objects
          minder.append((0,0))

          # ... but also the external order one derivative
          lleft = 'N' + left.replace('I','')
          tempden.append(ReconstructDensity(der, lap-1, lleft, right))
          temp_intermediate.append(True) # these are intermediate objects
          minder.append((0,0))

          rright = 'N' + right.replace('I','')
          tempden.append(ReconstructDensity(der, lap-1,  left, rright))
          temp_intermediate.append(True) # these are intermediate objects
          minder.append((0,0))

        elif(minder[j][1]>2):
          print ("Hephaestos cannot combine DENSUM=1 with high order derivatives yet.")
          sys.exit(1)
          
    # Complete the needed derivatives from the "maximal" number of derivatives
    temp_deriv_needed = PopulateDeriv(minder)
    #---------------------------------------------------------------------------
    # Pruning the list
    # A) removing duplicates
    # B) removing contractions when the full density will be calculated
    Densities_needed.append(tempden[0])
    deriv_needed.append([])
    intermediate_status.append(False)  # the first density is a "real" one

    for i in range(len(tempden)):
        (deri, lapi, lefti, righti, coupi, crossi) = \
                        ParseOperators(tempden[i],so.timelike)

        Found = False
        for j in range(len(Densities_needed)):
            (derj, lapj, leftj, rightj, coupj, crossj, foundind) = \
               ParseOperators(Densities_needed[j],so.timelike, findindices=True)
            #-------------------------------------------------------------------
            # Two densities are identical if the left- and right-operators are 
            # the same.
            # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
            # TODO: expand this detection to also check if two densities are the
            #       same up to symmetry transformations. 
            #       Example: D_I_NN = D_NN_I 
            #-------------------------------------------------------------------
            if(leftj == lefti and rightj == righti):
                # Signal that the density is already present
                Found = True
                # Add the possibility of new derivatives
                deriv_needed[j] = deriv_needed[j] + temp_deriv_needed[i]
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
                            l = foundind[coupj.index(cj)]
                        except ValueError:  
                            l = foundind[len(coupj) + crossj.index(cj)]

                        Densities_needed[j] = Densities_needed[j].replace(l, '')

            #-------------------------------------------------------------------
        if(not Found):
            add = tempden[i] 
            # Getting the duplicates out of the derivatives by adding these
            # explicitly to deriv_needed 
            deriv_needed.append(temp_deriv_needed[i])
            # Remove all of the derivatives from the top
            add = add.replace(derstring + '_','').replace(lapstring+'_', '')
            for l in sumindices:
                add = add.replace(derstring + l + '_','')
            Densities_needed.append(add)
            intermediate_status.append(temp_intermediate[i])

    # Clean up the deriv_needed array and add all combinations that might be 
    # necessary for calculating fields etc. This call is not strictly needed 
    # here, but it will become so later for the setting up the fields (where 
    # it will be called again).
    src_heph.heph_functional.PruneDeriv(deriv_needed)
    #---------------------------------------------------------------------------
    # Finding out how many derivatives we need to take of the spwfs
    derivative_order = 1    
    for den in Densities_needed:
        (der, lap, left, right, cpl, cross) = ParseOperators(den,so.timelike)
        ders = right.count('N')
        derivative_order = max(derivative_order, ders)
    # ... but also how much figure in the delta_action routine(s)
    pairing_action_derorder = 0
    for den in Densities_needed:
      if('P' not in den):
        continue
      (der, lap, left, right, cpl, cross) = ParseOperators(den,so.timelike)
      pairing_action_derorder = max(pairing_action_derorder, right.count('N'))
    #---------------------------------------------------------------------------
    # Finding out whether this functional file places calls to subroutine vmicro
    vmicro_found = False
    for cc in extra_calls:
      if('vmicro' in cc):
        vmicro_found = True
    #---------------------------------------------------------------------------
    # Print some information to STDOUT
    print (' Functional form taken from file %s'%fname)
    print (' Description from file:')
    print ( description.replace('#', tab))
    print (' Number of terms     : %d'%len(Functional_terms))
    if(so.timelike):
      print (' ! Attention: terms with time-odd densities dropped. ' )
    print (' Order of derivatives: %d'%derivative_order)
    print (' Microscopic pairing treatment: ', vmicro_found)
    print (' # Parameters        : %d'%len(paramparameters))
    #print   paramparameters
    for i in range(int(len(paramparameters)/3)):
        print ('  ', paramparameters[3*i:3*i+3])
    if(len(paramparameters)%3 != 0):
        print ('  ', paramparameters[3*(i+1):3*(i+1)+len(paramparameters)%3])

    return (description)
    
def regroup_terms():
  """
    Reorder terms such that all terms with the same structure, but different 
    isospin couplings are grouped together.
  """
  
  global Functional_terms, coupling_constants, isospin_indices, extra_calls
  global density_dependence, term_grouping, term_number
  
  # Copy the lists into temporary lists
  tempterms = Functional_terms 
  tempcoupl = coupling_constants
  tempiso   = isospin_indices
  tempddep  = density_dependence
  tempextra = extra_calls
  
  term_grouping      = []
  Functional_terms   = []
  isospin_indices    = []
  coupling_constants = []
  density_dependence = []
  extra_calls        = []
  
  # First: find all terms with a unique structure  
  unique_terms = []
  for term in tempterms: 
    term_grouping.append(0)
    if(term not in unique_terms):
      unique_terms.append(term)
      term_number[term] = -1
      
  # Then, repopulate the Functional_terms with all terms "grouped"
  for l,uterm in enumerate(unique_terms):
    for k,term in enumerate(tempterms):
      if(term == uterm):
        Functional_terms.append(term)
        coupling_constants.append(tempcoupl[k])
        density_dependence.append(tempddep[k])
        isospin_indices.append(tempiso[k])
        extra_calls.append(tempextra[k])
        
        term_number[uterm] = term_number[uterm] + 1

  for k,aterm in enumerate(Functional_terms):
    for l,bterm in enumerate(Functional_terms[:k]):
      if(bterm == aterm):
        term_grouping[k] = term_grouping[k] +1

def PopulateDeriv(minder):
    """
      Complete the set of derivatives needed by populating an array with all
      lower order derivative combinations.
    """
    
    derivs = []
    for j,combo in enumerate(minder):
      derivs.append([])
      options = itertools.product(range(combo[0]+1), range(combo[1]+1))
      for opt in options:
          derivs[j].append(opt)
    
    return derivs
    
def PruneDeriv(derivs):
    """
      Add all of the possible combinations with less derivatives and laplacians,
      so that we can build the eventually needed combinations.
    """
    # a) remove duplicates in the list of needed derivative combinations
    # b) Order the list in increasing level of operations
    for i in range(len(derivs)):
        derivs[i] = list(set(derivs[i]))
        derivs[i] = sorted(derivs[i])

def ReadFunctional(fname):
  """
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Read the details of the EDF terms, their structure and coupling constants, 
     from the file named fname. 

     Such a file should be composed of 
      
      Part 1: description of the functional (in words)
      Part 2: enumeration of the parameters
      Part 3: !TERMS => signalling the end of Part 2 and the start of Part 1
      Part 4: term specification
    
     All of these parts can be interspersed with lines starting with '!'; these
     are comment lines and do not influence the code generation in any way. 
      
      Part 1: description of the type of functional, which will be included in
              the Hephaestos output.  
              Lines need to start with '#'    

      Part 2: enumeration of all parameters that should be read from a .param
              file. Entries should be separated by ';'.

      Part 3: '!TERMS' (no modification, EVER)           

      Part 4: line-by-line specification of all the terms in the functional. 
              These should have the form 
              
              E_[D1]_[D2]_[D3]_[D4] ; C ; alpha ; iso_1 ; iso_2 ; iso_3 ; iso_4 
                 (1)                 (2)  (3)      (4)
              
              (1)    enumeration of the densities in the term, including the way
                     they are coupled. Example: 
                     
                      E_D_I_Sm_D_I_Sm =  \sum_{mu=x/y/z} s_mu(r) s_mu(r)
                      
                     this version of the code allows for
                      (a) bilinear (two densities)
                      (b) trilinear
                      (c) quadrilinear terms

               (2)   coupling constant of the term. Can be given in terms of the 
                     Cc function coded in the Fortran templates. 
                     
               (3)   density dependence of the FIRST density, D1. 
                     If alpha != 1, the code will enforce iso_1 to be zero
               
               (4)   isospin indices of the densities; 0, 1, 'p' or 'n'.
                     normal densities should have isospin indices (0,1) and 
                     pairing densities should have p/n indices ('p', 'n')
                     
     Special exception: for pairing terms, the code can read ONE ADDITIONAL 
     entry into the row, useable to make the code call an additional extra
     routine to calculate energies and fields.
    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  """
  global Functional_terms, extra_calls

  def clean(a):
    # Quick'n'dirty string cleaning routine, strips spaces and newlines
    return a.replace(' ', '').replace('\n', '')

  description      = ''
  termsstart       = 0
      
  with open(fname, 'r') as f:
    for line in f:
      try:
      
        if(len(line.split()) == 0):
          # Forget about empty lines
          continue
        elif(line[0] == '#'):
          #-signs indicate PART 1, description of the EDF             
          description = description + line
          continue
        elif(line[0:6] == '!TERMS'):
          # Signal that the parameter specification, PART 2 is over.
          termsstart  = 1   
          continue            
        elif(line[0] == '!'):
          # comment line, don't do anything with it
          continue
          
        if(termsstart == 0):
          # Parse the parameters of the functional in PART 2
          split = line.split(';')
          for s in split:
             paramtype, param = identify_param(s)
             paramparameters.append(param)
             paramtypes.append(paramtype)
        elif(termsstart == 1):
          #  Start the actual terms of the functional in PART 4
          split = line.split(';')
          Functional_terms.append  (clean(split[0]))
          coupling_constants.append(clean(split[1]))
          density_dependence.append(clean(split[2]))                      
          (densities,trash) = ParseDensities(clean(split[0]))      

          iso_list = []
          for k in range(3,3+len(densities)):
            iso = clean(split[k])
            if(iso != '0' and iso != '1' and iso != 'p' and iso != 'n'):
              raise IndexError
            iso_list.append(iso)
          isospin_indices.append(iso_list)
          
          if(len(split)>3+len(densities)):
            extra = clean(split[-1])
            extra_calls.append(extra)            
          else:
            extra_calls.append('')            
      except IndexError:
          print ('Problem reading the following line in the func file.')
          print (line)
          sys.exit(1)      


  return description  

def RemoveTimeOddTerms(fam_active):
    """
      We remove all terms from the EDF specification that contain
      time-odd densities IF we are targetting a mean-field code.

      IF targetting a FAM executable, we keep ALL terms.

      Input:
      ------
        fam_active : if True, keep all terms.
    """
    global Functional_terms, coupling_constants, density_dependence
    global density_dependence, isospin_indices, extra_calls
    

    toremove = []
    if(not fam_active):
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
        sys.exit(1)
      
      if(timeodd):
        toremove.append(i)

    tempterms   = Functional_terms
    tempcc      = coupling_constants
    tempdd      = density_dependence
    tempiso     = isospin_indices
    temp_extra  = extra_calls

    Functional_terms      = []
    coupling_constants    = []
    density_dependence    = []
    isospin_indices       = [] 
    extra_calls           = []
            
    for j in range(len(tempterms)):
       if (j not in toremove):
          Functional_terms.append(tempterms[j])
          coupling_constants.append(tempcc[j])
          density_dependence.append(tempdd[j])
          isospin_indices.append(tempiso[j])
          extra_calls.append(temp_extra[j])


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

    for l in sumindices:
      for x in crossindices:
        for i in range(len(densities)): 
            if( derstring + x + l in densities[i]) :
                # Remove the coupling if it involves derivatives
                densities[i] = densities[i].replace(l, '')
                for j in range(len(densities)):
                    densities[j] = densities[j].replace(l, '')

    # Remove vector coupling indices that might remain
    for l in crossindices:
        for i in range(len(densities)): 
            densities[i] = densities[i].replace(l, '')
    #---------------------------------------------------------------------------
    return (densities, coupling)

def ProcessParameterization(fname, src, target):
    """
     Process the parameterization.f90 file to include the different parameters.
    """

    # Go and get the relevant snippets of FORTRAN CODE to combine. 
    import src_heph.fortran_templates.ProcessParameterization_templates as ts

    decl      = ''
    readparam = ''
    printparam= ''
    checkparam= ''
    resetparam= ''
    bcastparam= ''
    
    for k,s in enumerate(paramparameters):
        dic= {}
        dic['PARAM']     = s
        
        if(paramtypes[k] == 'real'):
          decl      = decl      + ts.decl_real.substitute(dic)
          printparam= printparam+ ts.print_real.substitute(dic)
          bcastparam= bcastparam+ ts.bcast_param_real.substitute(dic)
        else:
          decl      = decl      + ts.decl_int.substitute(dic)
          printparam= printparam+ ts.print_int.substitute(dic)
          bcastparam= bcastparam+ ts.bcast_param_int.substitute(dic)
          
        readparam = readparam + ts.read.substitute(dic)
        checkparam= checkparam+ ts.check_a.substitute(dic)
        checkparam= checkparam+ ts.check_b.substitute(dic)
        checkparam= checkparam+ ts.check_c.substitute(dic)
        checkparam= checkparam+ ts.check_d.substitute(dic)
        checkparam= checkparam+ '\n'
        resetparam= resetparam+ ts.reset.substitute(dic)
        
    # Remove the trailing comma and add line-end
    readparam = readparam[:-5] + '\n'    
    
    dic= {}
    dic['PARAMDECL']   = decl
    dic['READPARAMS']  = readparam
    dic['PRINTPARAMS'] = printparam
    dic['CHECKPARAMS'] = checkparam
    dic['RESETPARAMS'] = resetparam
    dic['BCASTPARAMS'] = bcastparam

    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))  

def ProcessFunctional(fname, src, target, so, oldso, ph_pp_decoupl, 
                      fam_active, density_spwf_summation):
    """
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Master routine calling the other ones to generate a functional.
     
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
     Input:
      - fname: name of the functional file, without the directory structure
      - src:  source directory where the templates are stored
      - target: target directory where the generated code will be placed
      - so:   set of symmetry options
      - oldso: set of symmetry options for input wf files
      - ph_pp_decoupl: boolean determining whether pairing terms are decoupled
                       from the particle-hole part of the functional.
      - fam_active:  boolean determining whether we are building a mean-field 
                     or a finite-amplitude linear response code.
      - density_spwf_summation: boolean determining whether derivatives of 
                               densities are calculated through summation over 
                               spwfs or derivative calls.
     Output:
      - pot_declaration: a (large) string containing the fortran code for the
                         declaration of mean-field potentials in vectors.f90

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    """

    declaration     = ''
    pot_declaration = ''
    calculation     = ''
    form            = ''
    printing        = ''
    calccoef        = ''
    printcoef_ph    = ''
    printcoef_pair  = ''
    sumtotal_even   = ''
    sumtotal_odd    = ''
    sumtotal_bi     = ''
    sumtotal_tri    = ''
    sumtotal_quad   = ''
    sumtotal_dd     = ''
    
    precond         = ''
    erear           = ''
    writing         = ''
    writing_hdf5    = ''
    reading         = ''
    reading_hdf5    = ''
    cleaning        = ''

    init            = ''
    add             = ''
    multiply        = ''
    
    inproduct       = ''
    
    pairtotal_neutron = ''
    pairtotal_proton  = ''
    #---------------------------------------------------------------------------
    # Generate the terms in the functional
    Quadri = False
    
    for i in range(len(Functional_terms)): 
        # Check if we have a quadrilinear term
        (tempden, coupling) = ParseDensities(Functional_terms[i])
        if(len(tempden) == 4):
          Quadri = True
   
        # Generate a bunch of strings to insert into the FORTRAN code for this
        # particular term 
        (d,c,p,cc, pc_ph, pc_pair, st,pt,er, T) = \
          GenTermExpression(Functional_terms[i], i, term_grouping[i],    
                        term_number[Functional_terms[i]], coupling_constants[i], 
                 isospin_indices[i], density_dependence[i], extra_calls[i], so, fam_active)
                 
        if( d != ''):
          declaration = declaration + d + '\n'
        calculation = calculation + c + '\n'
        printing    = printing    + p + '\n'
        calccoef    = calccoef    + cc+ '\n'
        printcoef_ph  = printcoef_ph    + pc_ph + '\n'
        printcoef_pair= printcoef_pair  + pc_pair+ '\n'

        # Decide in which of the subtotals of the energy the term belongs
        if(T):
          # Term to be added to the time-even subtotal
          sumtotal_even = sumtotal_even    + st+ '&\n'
        else:
          # Or to the time-odd subtotal
          sumtotal_odd  = sumtotal_odd     + st+ '&\n'

        (tempden, coupling) = ParseDensities(Functional_terms[i])
        #-----------------------------------------------------------------------
        # Figure out if there is a pairing density in this term
        # Attention: these calculations are all intended for the calculation 
        #            of the SPWF-energy. They are not intended as a calculation
        #            of ALL terms that are bi/tri/quadrilinear or density-dependent.
        pairing_term = False
        for k in range(len(tempden)):
          if ('P' in tempden[k]):
             pairing_term = True
        if((not pairing_term) or (not ph_pp_decoupl)):
          if(density_dependence[i] == '1'):
            if(len(tempden) == 2):
              sumtotal_bi   = sumtotal_bi   + st + '&\n'
            elif(len(tempden) == 3):
              sumtotal_tri  = sumtotal_tri  + st + '&\n'
            elif(len(tempden) == 4):
              sumtotal_quad = sumtotal_quad + st + '&\n'
          else:
            sumtotal_dd = sumtotal_dd + st + '&\n'
          
        if('P' in Functional_terms[i]):    
          if('p' in isospin_indices[i]):
            pairtotal_proton = pairtotal_proton   + pt+ '&\n'
          else:
            pairtotal_neutron= pairtotal_neutron  + pt+ '&\n'
                  
        if('P' in Functional_terms[i] and ph_pp_decoupl):
          # Drop the rearrangement terms due to density-dependent pairing terms
          pass
        else:
          erear       = erear       + er
    
    sumtotal_even  = rreplace( sumtotal_even, '&\n', '', 1)
    sumtotal_odd   = rreplace( sumtotal_odd , '&\n', '', 1)
    
    sumtotal_bi    = rreplace( sumtotal_bi  , '&\n', '', 1)
    sumtotal_tri   = rreplace( sumtotal_tri , '&\n', '', 1)
    sumtotal_quad  = rreplace( sumtotal_quad, '&\n', '', 1)
    sumtotal_dd    = rreplace( sumtotal_dd  , '&\n', '', 1)

    pairtotal_neutron = rreplace( pairtotal_neutron   , '&\n', '', 1)
    pairtotal_proton  = rreplace( pairtotal_proton    , '&\n', '', 1)

    if(len(pairtotal_proton) == 0):    
      pairtotal_proton = '0'
    if(len(pairtotal_neutron) == 0):    
      pairtotal_neutron = '0'

    #---------------------------------------------------------------------------
    # Generate the fields of the single-particle hamiltonian
    (fielddec,fieldini,fieldcalc,fieldcalc_perturbed, fieldprecon,fieldwrite, \
    fieldwrite_hdf5,fieldread, fieldread_hdf5,fieldadd,fieldmultiply,          \
    fieldinproduct,fieldINMk2, fieldINMk4) \
                                 =  GenerateFields(so, oldso,ph_pp_decoupl, fam_active)
    pot_declaration = pot_declaration + fielddec + '\n'
    writing     = writing     + fieldwrite 
    writing_hdf5     = writing_hdf5     + fieldwrite_hdf5
    reading     = reading     + fieldread
    reading_hdf5     = reading_hdf5     + fieldread_hdf5
    precond     = precond     + fieldprecon
    init        = init        + fieldini
    add         = add         + fieldadd
    multiply    = multiply    + fieldmultiply
    inproduct   = inproduct   + fieldinproduct
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
          SkyrmeAction = SkyrmeAction +                          \
                                       GenerateAction(field, 1, so)
          SkyrmeAction = SkyrmeAction +                          \
                                       GenerateAction(field,-1, so)
      else:
          SkyrmeAction = SkyrmeAction + GenerateAction(field, 0, so)

    #---------------------------------------------------------------------------
    # Generate the expressions for the actions of the pairing fields.
    PairingAction           = ''
    for field in  Pairing_Fields_needed:
      (left,right,coupling,cross) = ParseOperatorsField(field, so.timelike)
      PairingAction = PairingAction + GenerateAction(field, 0, so)

    #---------------------------------------------------------------------------
    # Now make sure all of the lines are not too long for compilation.
    declaration       = LineFormat(declaration)
    pot_declaration   = LineFormat(pot_declaration)
    calculation       = LineFormat(calculation)
    printing          = LineFormat(printing)
    calccoef          = LineFormat(calccoef)
    printcoef_ph      = LineFormat(printcoef_ph)
    sumtotal_even     = LineFormat(sumtotal_even)
    sumtotal_odd      = LineFormat(sumtotal_odd)

    fieldcalc           = LineFormat(fieldcalc)
    fieldcalc_perturbed = LineFormat(fieldcalc_perturbed)
    
    precond           = LineFormat(precond)
    SkyrmeAction      = LineFormat(SkyrmeAction)
    PairingAction     = LineFormat(PairingAction)
    erear             = LineFormat(erear)
    reading           = LineFormat(reading)
    writing           = LineFormat(writing)
    init              = LineFormat(init)
    add               = LineFormat(add)
    multiply          = LineFormat(multiply)
    inproduct         = LineFormat(inproduct)
    fieldINMk2        = LineFormat(fieldINMk2)
    fieldINMk4        = LineFormat(fieldINMk4)
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Substitute into the functional.f90 file.  
    dic={}

    dic['NTERMS']                 = len(Functional_terms)
    dic['DECLARATION']            = declaration
    dic['CALCULATION']            = calculation
    dic['INIPOTENTIALS']          = init
    dic['MULTIPLY_POTENTIALS']    = multiply
    dic['ADD_POTENTIALS']         = add
    dic['PRINT']                  = printing
    dic['CALCCOEF']               = calccoef
    dic['PRINTCOEF_PH']           = printcoef_ph
    dic['PRINTCOEF_PAIR']         = printcoef_pair

    dic['TOTAL_EVEN']     = sumtotal_even
    if(len(sumtotal_odd)>1):
      dic['TOTAL_ODD']      = sumtotal_odd
    else:
      dic['TOTAL_ODD']      = '0.0d0'

    dic['TOTAL_BI']     = sumtotal_bi
    if(len(sumtotal_tri)>1):
      dic['TOTAL_TRI']      = sumtotal_tri
    else:
      dic['TOTAL_TRI']      = '0.0d0'
    if(len(sumtotal_quad)>1):
      dic['TOTAL_QUAD']      = sumtotal_quad
    else:
      dic['TOTAL_QUAD']      = '0.0d0'

    if(len(sumtotal_dd)>1):
      dic['TOTAL_DD']      = sumtotal_dd
    else:
      dic['TOTAL_DD']      = '0.0d0'

    dic['TOTALPAIR_NEUTRON']= pairtotal_neutron
    dic['TOTALPAIR_PROTON'] = pairtotal_proton

    dic['K2POT'] = fieldINMk2
    dic['K4POT'] = fieldINMk4
    
    dic['CALCPOTENTIALS']           = fieldcalc
    dic['CALCPOTENTIALS_PERTURBED'] = fieldcalc_perturbed
    dic['POTENTIALPRECON']          = precond
    if('update' not in precond):
        dic['PRECON_ACTIVE'] = '!'
    else:
        dic['PRECON_ACTIVE'] = ' '

    if('update' not in precond):
        dic['PRECON_ACTIVE'] = '!'
    else:
        dic['PRECON_ACTIVE'] = ' '

    dic['SKYRMEACTION']   = SkyrmeAction
    dic['PAIRINGACTION']  = PairingAction
    dic['EREAR']          = erear
    dic['FUNC_NAME']      = func_name
    dic['POTENTIALNUMBER']= len(Densities_needed)
    dic['WRITEPOTENTIALS']= writing
    dic['WRITEPOTENTIALS_HDF5']= writing_hdf5
    dic['READPOTENTIALS'] = reading
    dic['READPOTENTIALS_HDF5'] = reading_hdf5
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    # Making sure to (un)comment the parts of the interfaces of the routines of
    #  sphamil, delta_action and the derivative routines. 
    if('dtemp' in SkyrmeAction):
      dic['DTEMPSPH'] = ' '
    else:
      dic['DTEMPSPH'] = '!'
    if('ddtemp' in SkyrmeAction):
      dic['D2TEMPSPH'] = ' '
    else:
      dic['D2TEMPSPH'] = '!'
    if('dddtemp' in SkyrmeAction):
      dic['D3TEMPSPH'] = ' '
    else:
      dic['D3TEMPSPH'] = '!'
    if('laptemp' in SkyrmeAction):
      dic['LAPTEMPSPH'] = ' '
    else:
      dic['LAPTEMPSPH'] = '!'

    if('dtemp' in PairingAction):
      dic['D1TEMPDELTA'] = ' '
    else:
      dic['D1TEMPDELTA'] = '!'
    if('ddtemp' in PairingAction):
      dic['D2TEMPDELTA'] = ' '
    else:
      dic['D2TEMPDELTA'] = '!'
    if('dddtemp' in PairingAction):
      dic['D3TEMPDELTA'] = ' '
    else:
      dic['D3TEMPDELTA'] = '!'
    if('laptemp' in PairingAction):
      dic['LAPTEMPDELTA'] = ' '
    else:
      dic['LAPTEMPDELTA'] = '!'
    
    if(Quadri):
      dic['QUADRI'] = ' '
    else:
      dic['QUADRI'] = '!'

    if(derivative_order == 1):
      dic['N2'] = ' '    
      dic['N3'] = '!'
    elif(derivative_order == 2): 
      dic['N2'] = ' '
      dic['N3'] = '!'
    elif(derivative_order == 3):
      dic['N2'] = '!'
      dic['N3'] = ' '

    if(derivative_order == 1):
      dic['N1'] = '!'
    else:
      dic['N1'] = ' '

    if(pairing_action_derorder == 0):
      dic['N1DELTA']  = '!'
      dic['N2DELTA']  = '!'
      dic['N3DELTA']  = '!'
      dic['SYMDELTA'] = '!'
    elif(pairing_action_derorder == 1):
      dic['N1DELTA']  = ' '
      dic['N2DELTA']  = '!'
      dic['N3DELTA']  = '!'
      dic['SYMDELTA'] = ''
    elif(pairing_action_derorder == 2):
      dic['N1DELTA']  = ' '
      dic['N2DELTA']  = ' '
      dic['N3DELTA']  = '!'
      dic['SYMDELTA'] = ' '
    elif(pairing_action_derorder == 3):
      dic['N1DELTA']  = ' '
      dic['N2DELTA']  = ' '
      dic['N3DELTA']  = ' '
      dic['SYMDELTA'] = ' '

    if('D_Nm_Nm' not in Densities_needed):
      if('D_N_N' not in Densities_needed):
        dic['TAUSCALAR'] = '!'
        dic['TAUTENSOR'] = '!'
        dic['NOTAU']     = '!'
      else:
        dic['TAUSCALAR'] = '!'
        dic['TAUTENSOR'] = ' '
        dic['NOTAU']     = ' '
      if('D_N_N' not in Densities_needed):
        dic['TAUSCALAR'] = '!'
        dic['TAUTENSOR'] = '!'
        dic['NOTAU']     = '!'
      else:
        dic['TAUSCALAR'] = '!'
        dic['TAUTENSOR'] = ' '
        dic['NOTAU']     = ' '
    else:
      dic['TAUSCALAR'] = ' '
      dic['TAUTENSOR'] = '!'
      dic['NOTAU']     = '!'

    if(so.timelike):
      dic['NTR'] = '!'
      dic['TR']  = ''
    else:
      dic['NTR'] = ''
      dic['TR']  = '!'

    if(fam_active):
      dic['FAM'] = '1'
    else:
      dic['FAM'] = '0'
      
    dic['PVECTORINPRODUCT'] = inproduct
    
    with open(src+fname, 'r') as template:
      with open(target+fname, 'w') as generated:
        for line in template:
          generated.write(Template(line).substitute(dic))

    return pot_declaration

def GenTermExpression(term,index,un_index,tnumber, ccoef, isoc, ddep, extra,so, fam_active):
    """
     Generate the FORTRAN expressions to calculate the terms in the functional.
     
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

     Input:
      term : string containing the actual expression
             example: E_D_Nm_Nm_D_No_No
      index: number of the term in the complete list of terms
      un_index : index of the term in the grouping of terms with same 
                  structure. If negative, this term is unique.
      t_number : total number of terms with this structure
      ccoef: string containg the coupling constant of the term
      isoc : isospin coupling of the terms
      ddep : density dependence of the FIRST density in the term
      extra: extra function call to perform
      so   : symmetry options
      fam_active :  generate code for mean-field (False) or FAM (True)
                    calculations

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      
     Output: a ton of strings that should be inserted in the FORTRAN templates
     
      declaration   : declaration of the (1) coupling constants
                                         (2) terms in the energy
      calculation   : FORTRAN statements to calculate the terms
      printing      : FORTRAN statements to print the content of the terms
      calccoef      : FORTRAN statement that calculates the coupling constant

      printcoef_ph  :
      printcoef_pair:

      sumtotal      : 
      pairtotal     :
      erear         : 
      
      timerev       :

      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

    """
    global sumindices
    
    # This module contains all of the little snippets of FORTRAN code that 
    # need to be combined in an intelligent way by this function. 
    import src_heph.fortran_templates.GenTermExpression_templates as ts

    declaration = ''
    calculation = ''
    printing    = ''
    
    #---------------------------------------------------------------------------
    # Parse the term of the functional.
    # The result is 
    #   tempden : a list of densities that make up the term
    #   coupling: the set of scalar couplings of Cartesian indices in the term
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
    # symmetries and terms, this is achieved somewhere else in Hephaestos.
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
        for s1 in sumindices:
          if(densities[i].count(s1) == 2):
            altterm = altterm.replace(s1,'')
            
            # WR 14/04/22: I'm unsure why I made things this 
            #              complicated. For future reference, 
            #              the below code fails when multiple
            #              couplings are present in one density
            #              Example : D_I_I_D_NmNm_NkNk term in
            #              N2LO functionals.
            # Replace internal couplings
#            for s2 in sumindices:
#              tryout = densities[i].replace(s1, s2)
#              nosum  = densities[i].replace(s1, '')
#              altterm = altterm.replace(tryout,nosum)           
#              print (s1, s2, densities[i].count(s1), altterm)    
    (rubbish, true_coupling) = ParseDensities(altterm)
       
    dic = {}
    index_encountered=0
    name = ''
    for i in range(len(densities)):
            name = name + '_' + densities[i]

    dic ['TERM' ] = term
    dic ['CPCTE'] = 'coupl_constant(%d)'%(index+1)  #'B' + dic ['TERM'][1:]   
    dic ['GROUPINDEX']  = ''
    dic ['GROUPNUMBER'] = '(%d)'%(tnumber+1)
    if( tnumber >= 0):
      dic['GROUPINDEX'] = '(%d)'%(un_index + 1) 
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
      true_args = []
      for a in args:
        for va in vec_args:
          true_args.append(a+va)
           
    if(un_index == 0):  
      # only construct a declaration for the first term in a set
      declaration = ts.decl.substitute(dic) 
    calculation = ts.comment.substitute(dic)
    calculation = calculation + ts.calc_z.substitute(dic)
    
    for arg in true_args: 
        dic['EDENT'] = ''

        sign      = +1
        prevorder =  0 
        for i in range(len(densities)):
            isodic = {}
            isodic['DEN'] = densities[i]
            isodic['ISO'] = Isospinindices(isoc[i])   
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
           
            if(i == 0 and ddep != '1'):
              # The first density for the first density in the term
              isodic['EXP']= ddep
              if(fam_active):
                dic['EDENT'] = dic['EDENT'] + ts.edent_DD_fam.substitute(isodic)+ '*'
              else:
                dic['EDENT'] = dic['EDENT'] + ts.edent_DD.substitute(isodic)+ '*'
            else:
              # No density dependence
              if(fam_active):
                dic['EDENT'] = dic['EDENT'] + ts.edent_fam.substitute(isodic)+ '*'
              else:
                dic['EDENT'] = dic['EDENT'] + ts.edent.substitute(isodic)+ '*'

            # Take out the final '*' which should not be necessary
            prevorder = prevorder + orders[i]
      
        if(sign == +1):
          dic['SIGN'] = '+'
        else:
          dic['SIGN'] = '-'
            
        dic['EDENT'] = dic['EDENT'][:-1]
        calculation = calculation + tab + '! indices = ' + str(arg) + '\n'

#        if('P' not in term ):
#        # Ordinary mean-field densities only in the term
        calculation = calculation + ts.calc_a.substitute(dic)
#        else:
#          # Pairing mean-field densities in the term.
#          calculation = calculation + ts.calc_a.substitute(dic)
            
    if(extra != ''):
      dic['EXTRA'] = extra
      calculation = calculation + ts.calc_extra.substitute(dic) 
    calculation = calculation + ts.calc_b.substitute(dic)
    calculation = calculation + '\n'

    dic['ISO_FULL'] = ''
    for k in range(4):
      try:
        dic['ISO_FULL'] = dic['ISO_FULL'] + "'%s',"%isoc[k]
      except IndexError:
        continue
    dic['ISO_FULL'] =     dic['ISO_FULL'][:-1] 

    # Format statements for printing in the FORTRAN code depend on whether
    # this term is bilinear, trilinear or quadrilinear
    if(len(densities) == 2):
      dic['FMT'] = 97
    if(len(densities) == 3):
      dic['FMT'] = 98
    if(len(densities) == 4):
      dic['FMT'] = 99
    
    if(len(densities) == 1):
      print ('This term was not parsed correctly.')
      print (term)
      print (densities)
      sys.exit(1)
    printing = ts.print.substitute(dic) 
    calculation = calculation + ts.end_comment 
    
    # Temporary
    erear         = ''
    printcoef_ph  = ''
    printcoef_pair= ''
    sumtotal      = ''
    pairtotal     = ''
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Put the coupling constant in the right spot
    dic['EXP']       = ccoef
    dic['ISO']       = str(isoc)
    dic['NCONSTANT'] = index + 1
    calccoef  = ts.calc_coef.substitute(dic)
    
    if('P' not in term):
        # Ordinary mean-field densities
        printcoef_ph  = ts.print_cpl_ph.substitute(dic)
        printcoef_pair= ''
    else:
        printcoef_ph  = ''
        printcoef_pair= ts.print_cpl_pair.substitute(dic) 
    
    sumtotal  = ts.sumtotal.substitute(dic)
    pairtotal = ts.pairtotal.substitute(dic)
    
    # Getting the contribution to the rearrangement energy if this is a 
    # density dependent term
    if( ddep == '1'):
        erear = '' 
    else:
        rearcoef = '+' + '0.5d0 * (%d - (%s))'%(2-len(densities)+1, ddep) 
        dic['REARCOEF'] = rearcoef
        erear = ts.rear.substitute(dic)
        
    return (declaration, calculation, printing, calccoef, printcoef_ph, 
                            printcoef_pair, sumtotal, pairtotal, erear, timerev)    

def identify_param(paramstring):
  """
     Identify a parameter from a string read from file.
  """    
  # First, strip sole "R" and "I"
  if('R ' == paramstring[0:2]):
    paramtype = 'real'
  elif('I ' == paramstring[0:2]):
    paramtype = 'integer'
  else:
    print ('Unrecognized parameter type.')
    print ('Offending entry: ', paramstring)
    sys.exit(1)

  # cleaning routine, strips spaces and newlines
  param = paramstring[2:].replace(' ', '').replace('\n', '')

  return paramtype, param
  
def rreplace(s, old, new, occurrence):
     li = s.rsplit(old, occurrence)  
     return new.join(li)
