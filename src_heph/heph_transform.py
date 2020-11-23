#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
from string                    import Template
from src_heph.heph_symmetries  import *
from math                      import copysign

def ProcessTransform(fname, src, target, so, oldso):
    """
      Generate the required Fortran code for the transformation of input 
      single-particle wavefunctions. 

      NONSPATIAL: .true. if we are breaking an antilinear, antihermitian 
                  conserved symmetry.

      TRANSFO_NONSPATIAL_1-4: transformation rules for the spinors under the 
                              antilinear, antihermitian symmetry  

    """

    temp = " %+d * temp(i,%d, si + wave)"

    dic = {}

    # Checking which symmetries we need to break
    tobreak = []
    for og in oldso.generators:
      found = False
      for ng in so.generators:
        if(og == ng):
          found = True
          break
      if(not found):
        tobreak.append(og)
    
    # First, check if we can find an antilinear, antihermitian symmetry
    nonspatial = False
    for g in tobreak:
      if( (not g.linear) and (not g.hermitian)):
        nonspatial = True
        s = g

    # Implement the correct transformation
    if(nonspatial):
      dic['NONSPATIAL'] = ".true."
      for i in range(4):
        ind = abs(s.permutation[i])
        sign = copysign(1, s.permutation[i])
        dic['TRANSFO_NONSPATIAL_%d'%(i+1)] = temp%(sign,ind)

    else:
      dic['NONSPATIAL'] = ".false."
      for i in range(4):
        dic['TRANSFO_NONSPATIAL_%d'%(i+1)] = "0.0d0"              

    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))


