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

def preprocess(fname, src, target, so , oldso, ph_pp_decoupl,
               density_spwf_summation):
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
      density_spwf_summation: Boolean. If .True., calculate the derivatives of
                              densities by summing 
    """
      
    if(fname=='compilation.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='sphericalharmonics.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='geninfo.f90'):
        ProcessGeninfo(fname, src, target, so)
        return
    if(fname=='derivatives.f90'):
        ProcessDerivatives(fname, src, target, so)
        return
    if(fname=='wavefunctions.f90'):
        ProcessWavefunctions(fname, src, target, so)
        return
    if(fname=='precondition.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='tantalus.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='run_single.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='run_mpi.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='multirun_example.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='parameterization.f90'):
        ProcessParameterization(fname, src, target)
        return
    if(fname=='functional.f90'):
        ProcessFunctional(fname, src, target,so, oldso, ph_pp_decoupl, 
                          density_spwf_summation)
        return
    if(fname=='nil8.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='scfiteration.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='basis_transform.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='moments.f90'):
        ProcessMoments(fname, src, target, so)
        return
    if(fname=='constants.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='diag.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
        return
    if(fname=='evolution.f90'):
        ProcessGeneric(fname, src, target, so)
        return
    if(fname=='IO.f90'):
        ProcessIO(fname, src, target, so, oldso)
        return
    if(fname=='coulomb.f90'):
        ProcessCoulomb(fname, src, target, so)
        return
    if(fname=='pairing.f90'):
        ProcessPairing(fname, src, target, so)
        return
    if(fname=='pairing_strengths.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
#    if(fname=='BCS.f90'):
#        os.system('cp ' + src + fname + ' ' + target + fname)
#        return
    if(fname=='pairingcutoffs.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='printing.f90'):
        ProcessGeneric(fname, src, target, so)
        return
    if(fname=='HFB_gradient.f90'):
        ProcessGeneric(fname, src, target, so)
        return
    if(fname=='HFB_direct.f90'):
        ProcessGeneric(fname, src, target, so)
        return
    if(fname=='HFB.f90' or fname == 'BCS.f90'): 
        #BCS.f90 and HFB.f90 have exactly the same needs in terms of 
        # preprocessing by Hephaestos
        ProcessHFB(fname, src, target, so)
        return
    if(fname=='hartree-fock.f90'):
        ProcessHartreeFock(fname, src, target, so)
        return
    if(fname=='folding.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='temperature_projection.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='momentsofinertia.f90'):
        ProcessGeneric(fname, src, target, so)
        return
    if(fname=='fission_MOI.f90'):
        ProcessFission_MOI(fname, src, target, so)
        return
    if(fname=='particleinabox.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='timing.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='transform.f90'):
        ProcessTransform(fname, src, target, so, oldso)
        return
    if(fname=='densities.f90'):
        ProcessDensities(fname, src, target, so, density_spwf_summation)
        return
    if(fname=='cranking.f90'):
        ProcessCranking(fname, src, target, so)
        return
    if(fname=='convergence.f90'):
        ProcessCranking(fname, src, target, so)
        return
        

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


def ProcessGeneric(fname, src, target, so):  
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

    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))   
 
        
