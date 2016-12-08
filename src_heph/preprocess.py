#-------------------------------------------------------------------------------
# Director of preprocessing of the FORTRAN files: knows which file gets which  
# make-over. 
#
import os
from string  import Template

def preprocess(fname, src, target):
      # Decides which routine to call on which file.
      
    if(fname=='compilation.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='geninfo.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
    if(fname=='wavefunctions.f90'):
        processwavefunctions(fname, src, target)
        return
    if(fname=='tantalus.f90'):
        os.system('cp ' + src + fname + ' ' + target + fname)
        return
          

def processwavefunctions(fname, src, target):
    
    dic={}
    dic['BLOCKS'] = 8
    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic)) 
                
