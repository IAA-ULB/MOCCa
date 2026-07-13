#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Module governing the wavefunctions module of the FORTRAN code.
import src_heph.heph_functional 
from src_heph.heph_symmetries  import *
from string                    import Template


def ProcessWavefunctions(fname, src, target, so, fam_active, dry_run=False):
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
            
  SXab  :  integers determining the behaviour under symmetry of the spwfs
            in different symmetry blocks.
    where 
      X  = X;Y;Z depending on the direction
      a  = block index, i.e. 1-4 depending on the symmetries of the spwf
      b  = component index, i.e. 1-4 depending on the spinor component
            we are dealing with

  PCON   | : Variables activating the relevant piece of code for the 
  PBROKEN|   calculation of the expectation values of Parity. This is 
              all "placeholder", as Hephaestos should in the future be able
              to generate this type of code itself.
  """

  from src_heph.heph_substitute         import substitute

  dic={}
  dic['N2'] = '!'
  dic['N3'] = ' '

  if(src_heph.heph_functional.derivative_order == 1):
      dic['N2'] = ' '    
      dic['N3'] = '!'
  elif(src_heph.heph_functional.derivative_order == 2): 
      dic['N2'] = ' '
      dic['N3'] = '!'
  elif(src_heph.heph_functional.derivative_order == 3):
      dic['N2'] = '!'
      dic['N3'] = ' '

  axes = ['X', 'Y', 'Z']
  for k in range(3):
    if(so.ReduceAxes[k] == 1):
      # This axes is not represented because of symmetry, hence we need to
      # generate a sign
      for B in range(1,5):

        #    A symmetry imposes a relation on an spwf
        #    
        #       [U psi] (x,y,z,c) = u psi(x,y,z,c)
        #                         = psi(ex,ey,ez,c')
        #      where  
        #       * U is a symmetry operator
        #       * u is the associated quantum number
        #       * ex, ey, ez are x,y,z with perhaps a sign
        #       * c and c', are components, i.e. 1,2,3,4

        mult = multiply_quantum_numbers(so.syms[k].permutation, so.combs[k],B)
        # => calling this routine lets us get the behaviour under the 
        #    symmetry operator for every component, i.e. it gets us the 
        #    number "u" for each block of spwfs (= a fixed set of quantum
        #    numbers.)

        for c in range(1,5):
          key = 'S' + axes[k] + '%d'%B + '%d'%c
          dic[key] = '%+2d'%(mult[c-1]/c)
    else:
      # This axis is completely represented in the calculation, we put 0
      for B in range(1,5):
        for c in range(1,5):
          key = 'S' + axes[k] + '%d'%B + '%d'%c
          dic[key] = ' 0'

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Ugly manual checking if parity is part of the generator set and 
  #  signalling this to the FORTRAN code
  symdic  = populatesymmetries()
  dic['PCON']    = '!'
  dic['PBROKEN'] = ' '

  for sym in so.generators:
    if (sym == symdic['P']):
      dic['PCON']    = ' '
      dic['PBROKEN'] = '!'

  if( fam_active ):
    dic['FAM']    = '1'
  else:
    dic['FAM']    = '0'
      
  if(not dry_run):
    substitute(src+fname, target+fname, dic)  


    
    
