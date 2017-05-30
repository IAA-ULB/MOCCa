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
# Director of preprocessing of the FORTRAN files: knows which file gets which make-over. 
# 
# Currently linked: 
#
#       densities.f90   <=> heph_densities.py
#       derivatives.f90 <=> heph_derivatives.py
#
#============================================================================================

import os
from string           import Template
from heph_densities   import ProcessDensities
from heph_derivatives import ProcessDerivatives
from heph_symmetries  import ReduceAxes
from heph_functional  import ProcessFunctional

def preprocess(fname, src, target):
      # Decides which routine to call on which file.
      
    if(fname=='compilation.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='geninfo.f90'):
        ProcessGeninfo(fname, src, target)
        return
    if(fname=='derivatives.f90'):
        ProcessDerivatives(fname, src, target)
        return
    if(fname=='wavefunctions.f90'):
        ProcessWavefunctions(fname, src, target)
        return
    if(fname=='tantalus.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='functional.f90'):
        ProcessFunctional(fname, src, target)
        return
    if(fname=='nil8.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='constants.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='diag.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='hartree-fock.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='IO.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='densities.f90'):
        ProcessDensities(fname, src, target)
        return

def ProcessGeninfo(fname, src, target):
    #=======================================
    # This one is already more complicated. 
    #
    #
    
    dic={}
    dic['NUMSYM'] = sum(ReduceAxes)
        
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))          

def ProcessWavefunctions(fname, src, target):
    
    dic={}
    dic['BLOCKS'] = 8
    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic)) 

        
