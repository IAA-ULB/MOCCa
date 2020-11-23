#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Module governing the wavefunctions module of the FORTRAN code.
from src_heph.heph_functional  import derivative_order
from string                    import Template


def ProcessWavefunctions(fname, src, target, so):
    """    
      Process the wavefunctions.f90 file of the Fortran code.       
        
      A 3D calculation, respecting some subgroup of the D^T(D)_2h groups, 
      will always have single-particle wavefunctions that can be divided into 
      at most eight blocks according to their behaviour under symmetries. 

      2 for protons versus neutrons
      2 for a conserved hermitian, linear symmetry (parity in EV8)
      2 for a conserved antihermitian, linear symmetry (z-signature in EV8) 

      A conserved antilinear, antihermitian symmetry will then allow for the 
      elimination of one of these sets in a practical calculation.

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      N2/N3 :  decides which derivative routines to comment/uncomment 
               depending on the order of derivatives in the functional
  
      ININX/NY/NZ: number of points in the box for the nilsson initialization

    """

    dic={}
    dic['N2'] = '!'
    dic['N3'] = ' '

    if(derivative_order == 1):
        dic['N2'] = ' '    
        dic['N3'] = '!'
    elif(derivative_order == 2): 
        dic['N2'] = ' '
        dic['N3'] = '!'
    elif(derivative_order == 3):
        dic['N2'] = '!'
        dic['N3'] = ' '

    dic['ININX'] = "nx/2"
    dic['ININY'] = "ny/2"
    dic['ININZ'] = "nz/2"
    if( so.ReduceAxes[0] == 1):
      dic['ININX'] = "nx"
    if( so.ReduceAxes[1] == 1):
      dic['ININY'] = "ny"
    if( so.ReduceAxes[2] == 1):
      dic['ININZ'] = "nz"

    nonspatial = False
    for g in so.generators:
        if(not g.linear and not g.hermitian):
          nonspatial = True
    if(nonspatial):
      dic['ININWT']  ='nwt'
      dic['ININWN']  ='nwn'
      dic['ININWP']  ='nwp'
    else:
      dic['ININWT']  ='nwt/2'
      dic['ININWN']  ='nwn/2'
      dic['ININWP']  ='nwp/2'
    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic)) 


    
    
