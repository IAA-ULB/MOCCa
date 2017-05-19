#-------------------------------------------------------------------------------
# Director of preprocessing of the FORTRAN files: knows which file gets which  
# make-over. 
#
import os
from string         import Template
from heph_densities import processdensities

def preprocess(fname, src, target):
      # Decides which routine to call on which file.
      
    if(fname=='compilation.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='geninfo.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='derivatives.f90'):
        processderivatives(fname, src, target)
        return
    if(fname=='wavefunctions.f90'):
        processwavefunctions(fname, src, target)
        return
    if(fname=='tantalus.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='functional.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
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
        processdensities(fname, src, target)
        return
          

def processwavefunctions(fname, src, target):
    
    dic={}
    dic['BLOCKS'] = 8
    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic)) 

def processderivatives(fname, src, target):
    #=======================================
    # This one is already more complicated. 
    #
    #
    
    dic={}
    dic['LINESIZEX'] = 20
    dic['LINESIZEY'] = 20
    dic['LINESIZEZ'] = 20
    dic['DX'] = 'D'
    dic['DY'] = 'D'
    dic['DZ'] = 'D'
        
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))
        
