#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Director of preprocessing of the FORTRAN files: knows which file gets which
# make-over. 
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
from src_heph.heph_IO            import ProcessIO, ProcessIO_wf
from src_heph.heph_cranking      import ProcessCranking
from src_heph.heph_multipoles    import ProcessMoments, ProcessFission_MOI
from src_heph.heph_coulomb       import ProcessCoulomb

# String variables that should be passed from one preprocessing function to another
vectors_densities  = ''
vectors_potentials = ''
memory_densities   = ''

def preprocess(fname, src, target, so , oldso, ph_pp_decoupl,
               fam_active, density_spwf_summation, dry_run=False):
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
      dry_run: Boolean. If True, do not actually process anything, but do print output.
                        If False, actually process the files but do not print output.
    """
    
    global vectors_densities
    global vectors_potentials
    global memory_densities

    # List of files that require no preprocessing at all
    no_process_files = [ 'compilation.f90', 'sphericalharmonics.f90', \
                         'precondition.f90', 'constants.f90', 'version.f90', \
                         'run_single.f90', 'run_mpi.f90', \
                         'multirun_example.f90', 'nil8.f90', \
                         'pairing_strengths.f90', 'pairingcutoffs.f90', \
                         'folding.f90', 'timing.f90', \
                         'fam_run.f90', 'fam_gmres.f90', \
                         'fam_testing.f90', 'hdf5_auxiliary.f90', \
                         'IO_aux.f90', 'basis_transform.f90' ]
    # List of files that can be processed by a generic preprocessor
    generic_process_files = [ 'tantalus.f90', 'scfiteration.f90', \
                              'printing.f90', 'HFB_gradient.f90', \
                              'HFB_direct.f90', 'momentsofinertia.f90', \
                              'evolution.f90', 'fam.f90' ]

    if(fname in no_process_files):
        # Preprocessing = is a simple copy operation for these files
        if(not dry_run):
            os.system('cp ' + src + fname + ' ' + target + fname)
        return
    elif(fname in generic_process_files):
        # A generic preprocessor that makes a few simple substitutions                        
        ProcessGeneric(fname,src,target,so,fam_active, dry_run)
        return
    else: 
        # For these files, there is more work to do!
        if(fname=='geninfo.f90'):
            ProcessGeninfo(fname, src, target, so)
        if(fname=='derivatives.f90'):
            ProcessDerivatives(fname, src, target, so)
        if(fname=='wavefunctions.f90'):
            ProcessWavefunctions(fname, src, target, so)             
        if(fname=='parameterization.f90'):
            ProcessParameterization(fname, src, target)                    
        if(fname=='functional.f90'):
            vectors_potentials = ProcessFunctional(fname, src, target,so, oldso, ph_pp_decoupl,fam_active, density_spwf_summation)
        if(fname=='vectors.f90'):
            vectors_densities, memory_densities = ProcessDensities('densities.f90', src, target, so, fam_active,  density_spwf_summation)
            vectors_potentials = ProcessFunctional('functional.f90', src, target,so, oldso, ph_pp_decoupl,fam_active, density_spwf_summation)
            ProcessVectors(src,target,so,vectors_densities,vectors_potentials, memory_densities,fam_active, dry_run) 
        if(fname=='moments.f90'):
            ProcessMoments(fname, src, target, so, fam_active)
        if(fname=='evolution.f90'):
            ProcessGeneric(fname, src, target, so, fam_active)
        if(fname=='IO.f90'):
            ProcessIO(fname, src, target, so, oldso, fam_active)  
        if(fname=='coulomb.f90'):
            ProcessCoulomb(fname, src, target, so, fam_active)             
        if(fname=='pairing.f90'):
            ProcessPairing(fname, src, target, so)                        
        if(fname=='HFB.f90' or fname == 'BCS.f90'): 
            # BCS.f90 and HFB.f90 have exactly the same needs in terms of 
            # preprocessing by Hephaestos
            ProcessHFB(fname, src, target, so)
        if(fname=='hartree-fock.f90'):
            ProcessHartreeFock(fname, src, target, so)
        if(fname=='momentsofinertia.f90'):
            ProcessGeneric(fname, src, target, so, fam_active)
        if(fname=='fission_MOI.f90'):
            ProcessFission_MOI(fname, src, target, so)
        if(fname=='transform.f90'):
            ProcessTransform(fname, src, target, so, oldso) 
        if(fname=='densities.f90'): 
            vectors_densities, memory_densities = ProcessDensities(fname, src, target, so, fam_active,  density_spwf_summation)
        if(fname=='cranking.f90'):
            ProcessCranking(fname, src, target, so)
        if(fname=='convergence.f90'):
            ProcessCranking(fname, src, target, so)
        if(fname=='IO_wf.f90'):
            ProcessIO_wf(fname, src, target, so, oldso, fam_active) 
        if(fname=='fam.f90'):
            ProcessGeneric(fname, src, target, so, fam_active)

def ProcessGeninfo(fname, src, target, so):
    """
     This one is already more complicated. 
    """
    from src_heph.heph_substitute import substitute

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

    substitute(src+fname, target+fname, dic)
    # with open(src+fname, 'r') as template:
    #     with open(target+fname, 'w') as generated:
    #         for line in template:
    #             generated.write(Template(line).substitute(dic))   


def ProcessGeneric(fname, src, target, so, fam_active, dry_run=False):
    """

    """
    from src_heph.heph_functional import derivative_order
    from src_heph.heph_substitute import substitute
    
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

    if(not dry_run):    
        substitute(src+fname, target+fname, dic)


def ProcessVectors(src, target, so, densities, potentials, memory_densities, fam_active, dry_run=False):
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
  from src_heph.heph_substitute import substitute


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

  if(not dry_run):
    substitute(src+'vectors.f90', target+'vectors.f90', dic)