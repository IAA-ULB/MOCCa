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
from string           import Template
from src_heph.heph_densities     import ProcessDensities
from src_heph.heph_derivatives   import ProcessDerivatives
from src_heph.heph_functional    import ProcessFunctional, ProcessParameterization
from src_heph.heph_wavefunctions import ProcessWavefunctions
from src_heph.heph_pairing       import ProcessHFB, ProcessHartreeFock
from src_heph.heph_transform     import ProcessTransform


def preprocess(fname, src, target, so , oldso):
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
        ProcessFunctional(fname, src, target)
        return
    if(fname=='nil8.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='scfiteration.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='moments.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='constants.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='diag.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
        return
    if(fname=='evolution.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='IO.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='coulomb.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='pairing.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='BCS.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='pairingcutoffs.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='printing.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='HFB.f90'):
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
        os.system('cp ' + src + fname + ' ' + target + fname)
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
        ProcessDensities(fname, src, target)
        return

def ProcessGeninfo(fname, src, target, so):
    """
     This one is already more complicated. 
    """    
    dic={}
    dic['NUMSYM'] = sum(so.ReduceAxes)
        
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))          
        
