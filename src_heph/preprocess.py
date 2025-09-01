#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Director of preprocessing of the FORTRAN files: knows which file gets which
# make-over. 
# 
# Currently linked: 
#
#       densities.f90   <=> heph_densities.py
#       derivatives.f90 <=> heph_derivatives.py
#       functional.f90  <=> heph_functional.py
#===============================================================================

import os
from string                      import Template
from src_heph.heph_densities     import ProcessDensities
from src_heph.heph_derivatives   import ProcessDerivatives
from src_heph.heph_functional    import ProcessFunctional, ProcessParameterization
from src_heph.heph_wavefunctions import ProcessWavefunctions
from src_heph.heph_pairing       import ProcessHFB, ProcessHartreeFock
from src_heph.heph_pairing       import ProcessPairing
from src_heph.heph_transform     import ProcessTransform
from src_heph.heph_IO            import ProcessIO
from src_heph.heph_cranking      import ProcessCranking
from src_heph.heph_multipoles    import ProcessMoments, ProcessFission_MOI
from src_heph.heph_coulomb       import ProcessCoulomb

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Almost all files are processed by a matching Hephaestos routine. The exception
# is the vectors.f90 file, which is processed with results from the
# heph_densities and heph_potentials modules. The vectors.f90 file contains
# the definition of the potential and density vectors, and thus needs the
# list of densities and potentials. To avoid having to call the entire machinery
# twice, I simply save the relevant strings here. This means it is important
# that densities.f90 and potentials.f90 get processed BEFORE vectors.f90.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
vectors_densities  = ''
vectors_potentials = ''
memory_densities   = ''

def preprocess(fname, src, target, so , oldso, ph_pp_decoupl,
               fam_active, density_spwf_summation):
    """
      Dispatching routine that selects the right preprocessing routine and
      additional info for every source code file.
    
      fname : filename of the Fortran template
      src   : directory of the templates
      target: directory to put the finished source code
      so    : SymmetryOption object, containing all the details on       
              the symmetries conserved during the calculation
      oldso : SymmetryOption object, containing all the details on
              the symmetries conserved ON THE INPUT WF FILE.
      ph_pp_decouple: Boolean. If True, do not include contributions of 
                      density-dependent pairing interactions to the potentials
                      associated with normal densities
      fam_active: Boolean. If True, the FAM is active and the code needs to be
                        prepared for complex densities - among other things.
      density_spwf_summation: Boolean. If .True., calculate the derivatives of
                              densities by summing 
    """
    
    global vectors_densities
    global vectors_potentials
    global memory_densities

    if(fname=='compilation.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='sphericalharmonics.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='geninfo.f90'):
        ProcessGeninfo(fname, src, target, so)
    if(fname=='derivatives.f90'):
        ProcessDerivatives(fname, src, target, so)
    if(fname=='wavefunctions.f90'):
        ProcessWavefunctions(fname, src, target, so)
    if(fname=='precondition.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='tantalus.f90'):
        ProcessGeneric(fname,src,target,so,fam_active)
    if(fname=='version.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='run_single.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='run_mpi.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='multirun_example.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='parameterization.f90'):
        ProcessParameterization(fname, src, target)
    if(fname=='functional.f90'):
        vectors_potentials = ProcessFunctional(fname, src, target,so, oldso, \
                                               ph_pp_decoupl,fam_active, density_spwf_summation)
    if(fname=='vectors.f90'):
        ProcessVectors(src,target,so,vectors_densities,vectors_potentials, \
                                     memory_densities,fam_active)
    if(fname=='nil8.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='scfiteration.f90'):
        ProcessGeneric(fname,src,target,so,fam_active)
    if(fname=='basis_transform.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='moments.f90'):
        ProcessMoments(fname, src, target, so)
    if(fname=='constants.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='diag.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='evolution.f90'):
        ProcessGeneric(fname, src, target, so, fam_active)
    if(fname=='IO.f90'):
        ProcessIO(fname, src, target, so, oldso, fam_active)
    if(fname=='coulomb.f90'):
        ProcessCoulomb(fname, src, target, so, fam_active)
    if(fname=='pairing.f90'):
        ProcessPairing(fname, src, target, so)
    if(fname=='pairing_strengths.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='pairingcutoffs.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='printing.f90'):
        ProcessGeneric(fname, src, target, so, fam_active)
    if(fname=='HFB_gradient.f90'):
        ProcessGeneric(fname, src, target, so, fam_active)
    if(fname=='HFB_direct.f90'):
        ProcessGeneric(fname, src, target, so, fam_active)
    if(fname=='HFB.f90' or fname == 'BCS.f90'): 
        # BCS.f90 and HFB.f90 have exactly the same needs in terms of 
        # preprocessing by Hephaestos
        ProcessHFB(fname, src, target, so)
    if(fname=='hartree-fock.f90'):
        ProcessHartreeFock(fname, src, target, so)
    if(fname=='folding.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='momentsofinertia.f90'):
        ProcessGeneric(fname, src, target, so, fam_active)
    if(fname=='fission_MOI.f90'):
        ProcessFission_MOI(fname, src, target, so)
    if(fname=='timing.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='transform.f90'):
        ProcessTransform(fname, src, target, so, oldso)
    if(fname=='densities.f90'):
        vectors_densities, memory_densities = ProcessDensities(fname, src, target, so, fam_active,  density_spwf_summation)
    if(fname=='cranking.f90'):
        ProcessCranking(fname, src, target, so)
    if(fname=='convergence.f90'):
        ProcessCranking(fname, src, target, so)
    if(fname=='fam.f90'):
        print('cp ' + src + fname + ' ' + target + fname)
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='fam_run.f90'):
        print('cp ' + src + fname + ' ' + target + fname)
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='fam_gmres.f90'):
        print('cp ' + src + fname + ' ' + target + fname)
        os.system('cp ' + src + fname + ' ' + target + fname)
    if(fname=='fam_testing.f90'):
        print('cp ' + src + fname + ' ' + target + fname)
        os.system('cp ' + src + fname + ' ' + target + fname)

def ProcessGeninfo(fname, src, target, so):
    """
     This one is already more complicated. 
    """    
    dic={}
    dic['NUMSYM'] = sum(so.ReduceAxes)
    dic['SYMSTRING'] = so.desc
    dic['SYMLEN'] = str(len(so.desc))

    dic['REDUX'] = str(so.ReduceAxes[0])
    dic['REDUY'] = str(so.ReduceAxes[1])
    dic['REDUZ'] = str(so.ReduceAxes[2])
    
    #NS: provides symmetry information to
    #    geninfo to calculate k_sh
    if(so.ReduceAxes[0] == 1):
       dic['LINESIZEX']    = '2*nx'
    else:
       dic['LINESIZEX']    = 'nx'

    if(so.ReduceAxes[1] == 1):
       dic['LINESIZEY']    = '2*ny'
    else:
       dic['LINESIZEY']    = 'ny'

    if(so.ReduceAxes[2] == 1):
       dic['LINESIZEZ']    = '2*nz'
    else:
       dic['LINESIZEZ']    = 'nz'

    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))   


def ProcessGeneric(fname, src, target, so, fam_active):
    """

    """
    from src_heph.heph_functional import derivative_order 
    
    global derivative_order

    dic = {}

    if(so.timelike):
      dic['TR'] = ''
      dic['NTR']= '!'
    else:
      dic['TR'] = '!'
      dic['NTR']= ''

    if(derivative_order == 1):
      dic['N2'] = ' '    
      dic['N3'] = '!'
    elif(derivative_order == 2): 
      dic['N2'] = ' '
      dic['N3'] = '!'
    elif(derivative_order == 3):
      dic['N2'] = '!'
      dic['N3'] = ' '

    if(fam_active):
      dic['FAM'] = 1
    else:
      dic['FAM'] = 0

    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))


def ProcessVectors(src, target, so, densities, potentials, memory_densities, fam_active):
  """
    Process the vectors.f90 template Fortran file to filled versions

    - vectors.f90     -> for mean-field calculations; real potentials and densities
    - vectors_FAM.f90 -> for FAM calculations; complex densities and potentials

    Note: densities, potentials and memory_densities are outputs from the ProcessDensities function.

    Input:
        src              : location of the original template
        target           : location to store the filled template
        so               : symmetry options of the calculations  [not useful right now]
        densities        : string containing the Fortran declaration of all individual local densities
        potentials       : string containing the Fortran declaration of all individual local potentials
        memory_densities : string containing the code for updating densities
        fam_active       : Boolean, if True, the FAM is active and vectors.f90 should declare the Coulomb 
                           quantities as complex.
  """

  dic = {}
  dic['DECLARATION']            = densities
  dic['DECLARATION_POTENTIALS'] = potentials
  dic['MEMORY_DENSITIES']       = memory_densities

  if(fam_active):
    # If the FAM is active, we need to declare the Coulomb quantities as complex
    dic['COULOMB_REAL']           = '!'
    dic['COULOMB_COMPLEX']        = ''
  else:   
    # If twe are constructing a mean-field code, the Coulomb quantities are real
    dic['COULOMB_REAL']           = ''
    dic['COULOMB_COMPLEX']        = '!'

  print ('Processing vectors.f90 template...')
  print (dic['DECLARATION'])

  with open(src+'vectors.f90', 'r') as template:
    with open(target+'vectors.f90', 'w') as generated:
        for line in template:
            generated.write(Template(line).substitute(dic))
