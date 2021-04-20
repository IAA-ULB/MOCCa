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

      TR                    : logical for time-reversal conservation

      NONSPATIAL            : .true. if we are breaking an antilinear, 
                               antihermitian conserved symmetry.

      TRANSFO_NONSPATIAL_1-4: transformation rules for the spinors under the 
                              antilinear, antihermitian symmetry  

      SPATIAL               : .true. if we are breaking a spatial symmetry
      EXPANDX               : .true. if we want to obtain the full X-axis
      EXPANDY               : .true. if we want to obtain the full Y-axis
      EXPANDZ               : .true. if we want to obtain the full Z-axis
      
    """

    #---------------------------------------------------------------------------
    # One simple Fortran template is needed
    temp = " %+d * wftarget(i,%d)"
    #---------------------------------------------------------------------------
    dic = {}

    if(so.timelike):
        dic['TR'] = ''
    else:
        dic['TR'] = '!'

    #---------------------------------------------------------------------------
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
    #---------------------------------------------------------------------------
    # First, check if we can find an antilinear, antihermitian symmetry
    # to be broken
    nonspatial = False
    for g in tobreak:
      if( (not g.linear) and (not g.hermitian)):
        nonspatial = True
        s = g

    # Implement the correct information for breaking the antilinear, 
    # antihermitian symmetry
    dic['NONSPATIAL'] = PythLogicalToFortranString(nonspatial) 
    if(nonspatial):
      for i in range(4):
        ind = abs(s.permutation[i])
        sign = copysign(1, s.permutation[i])
        dic['TRANSFO_NONSPATIAL_%d'%(i+1)] = temp%(sign,ind)
    else:
      for i in range(4):
        dic['TRANSFO_NONSPATIAL_%d'%(i+1)] = "0.0d0"              
    #---------------------------------------------------------------------------
    # Then we look into spatial symmetries
    expandX = False
    expandY = False
    expandZ = False
    if(so.ReduceAxes == oldso.ReduceAxes):
       spatial = False
    else:
      spatial = True
      if(so.ReduceAxes[0] != oldso.ReduceAxes[0]):
        expandX = True
      if(so.ReduceAxes[1] != oldso.ReduceAxes[1]):
        expandY = True
      if(so.ReduceAxes[2] != oldso.ReduceAxes[2]):
        expandZ = True
         
      if((expandX and expandY) or (expandX and expandZ) or (expandY and expandZ)):
        # Sanity check
        print ('Hephaestos cannot create executables that are able to break more than one spatial symmetry at the same time.')    
        quit()
          
    # Tell the Fortran code along which axes to expand the s.p.w.fs.
    dic['EXPANDX'] = PythLogicalToFortranString(expandX)         
    dic['EXPANDY'] = PythLogicalToFortranString(expandY)          
    dic['EXPANDZ'] = PythLogicalToFortranString(expandZ)          
    dic['SPATIAL'] = PythLogicalToFortranString(spatial)
    
    if(spatial):
      #-------------------------------------------------------------------------
      # Find out which symmetry relation we can utilise, i.e. which one was
      # conserved in old-symmetry-options, but is broken in new-symmetry-options
      s = symmetry([+1,+2,+3,+4],[+1,+1,+1], True, True)
      c = []
      for i,sym in enumerate(oldso.syms):
        Found = True
        for newsym in so.syms:
          if(not(newsym is None)):
            if(newsym == sym):
              Found = False
        if(Found):
          s = sym
          c = oldso.combs[i]
      if(len(c) == 0):
        print ("Big problem in ProcessTransform.")
        quit()
            

      # Now to figure out how to do the actual expansion in the requested axis
      for i in range(4):
        ind = abs(s.permutation[i])
        sign = copysign(1, s.permutation[i])
        dic['TRANSFO_Z_%d_B1'%(i+1)] = temp%(sign,ind)
    else:
      # Now to figure out how to the actual expansion
      dic['TRANSFO_SPATIAL_Z_1'] = '0.0d0'    
      dic['TRANSFO_SPATIAL_Z_2'] = '0.0d0'    
      dic['TRANSFO_SPATIAL_Z_3'] = '0.0d0'    
      dic['TRANSFO_SPATIAL_Z_4'] = '0.0d0'    


    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))


def PythLogicalToFortranString(l):
  if(l):
    return '.True.'  
  else:
    return '.False.'  
  
