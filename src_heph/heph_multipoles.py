#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# This module processes the 
#  (a) moments.f90 
#  (b) fission_MOI.f90 
# modules of MOCCa, tailoring several aspects of the calculation of multipole
# moments to the symmetry choices. 
#
# Important remarks:
# - - - - - - - - - - 
#   *) right now the code only figures out the symmetry restrictions on the 
#      multipole moment of the scalar density rho. This can be done without
#      reference to the symmetry properties of the density (since they are
#      trivial), but an extension to multipole moments of a different density
#      would need information on that.
#
#   *) The routines below check the behaviour of spherical harmonics under
#      reversal of individual coordinates (x -> -x, y -> -y, z -> -z). 
#      This is for sure sufficient in the immediate use cases of the code, but
#      might be reimagined in the future. 
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# We take spherical harmonics in the convention of Messiah, i.e. in terms
# of spherical coordinates theta and phi: 
#
#        Y^{l}_{m} = (-1)^{m}[(2*l+1)/4\pi (l-m)!/(l+m)!]^{1/2}
#                  *   P^{m}_{l}(cos(\theta)) e^{i m \phi}
#
# The P^{m}_{l} are the associated Legendre Polynomials and -l <= m <= l. 
# The multipole moments calculated by MOCCa are then
#
#        Q_{l m} = r^l Y^l_m(theta, phi).
#
# The nuclear multipole moments are then simply integrals of the nuclear density
#     
#     < Q_{l m} > = int d^3 r  Q_{lm}(r) rho(r)
# 
# These expectation values are not necessarily real. 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#
# This module sorts out which multipole moments are restricted by the symmetry
# choices. For example:
#
#         < Q_10 >  ~  int d^3 r  z rho(r) = 0 
# 
# if there is a reflection symmetry along the z-axis. 
#
# Which of all possible Q_{lm} are restricted by symmetry depends on 
# the choice of orientation of the nucleus with respect to the Cartesian axes.
# This is exactly the degree of freedom related to the correspondence between
# the spherical coordinates theta, phi and the Cartesian axes. 
# 
# There are six possibilities, and the distinction is made by two parameters:
#
#     QuantisationAxis = X/Y/Z (1,2,3)
#
# which determines which of the Cartesian coordinates determines the 
# polar angle theta
#    
#           X/Y/Z = r cos(theta)
#
# The second parameter
#
#     SecondaryAxis = 1/2
#
# simply determines the reference Cartesian axis for the angle \phi. 
# SecondaryAxis = "1" corresponds to lexographical order (x,y) while "2"
# corresponds to non-lexicographical order (y,x).
#
# The angle phi is determined by
#     \phi = atan2(x_2, x_1)
#
# Examples:
#    (QuantisationAxis, SecondaryAxis,Permutation)
#    (3,1, (x,y,z))
#    (3,2, (y,x,z))
#    (2,1, (x,z,y))
#    (2,2, (z,x,y))
#    (1,1, (y,z,x))
#    (1,2, (z,y,x))
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

from string                    import Template
from src_heph.heph_symmetries  import *
from scipy.special             import sph_harm
import numpy                   as np
from src_heph.heph_substitute import substitute

def ProcessMoments(fname, src, target, so, fam_active=False, dry_run=False):
  """
    Generate the required Fortran code to process the moments.f90 file. 
    
    Input: 
      fname : filename of the Fortran template
      src   : directory of the templates
      target: directory to put the finished source code
      so    : SymmetryOption object, containing all the details on the 
              symmetries conserved during the calculation.
      fam_active : Boolean, whether we are compiling a mf or fam executable.
      dry_run    : Boolean, if True, do not write any files (default: False)
  """
  from src_heph.heph_functional import derivative_order
  global derivative_order

  tab           = '    '
  template_list = Template( tab+'moment_list($ELL,$EMM,$IND) = $ON   ')
  template_comm = Template( '! $REIM Q_{ $ELL $EMM} \n')
  REIM  = ['Re', 'Im']
  
  # This routine is currently hardcoded to work up to and including l = 20 
  maxl = 20

  # We pick a point (1,2,3), and check what changing the sign on all coordinates
  # do to the spherical coordinates (which is determined by the quantisation 
  # and secondary axis, hidden in the symmetryoption object)
  r, theta_orig, phi_orig = CartToSpher(+1,+2,+3,so)
  r, theta_x, phi_x       = CartToSpher(-1,+2,+3,so)
  r, theta_y, phi_y       = CartToSpher(+1,-2,+3,so)
  r, theta_z, phi_z       = CartToSpher(+1,+2,-3,so)

  filling = ''
  dic = {}
  for l in range(1,maxl+1):
    filling = filling + tab + '! l = %d \n'%l    
    for m in range(0, l+1):
      # Attention: sph_harm in scipy takes arguments in what I would call
      #            a nonintuitive order.
      orig = sph_harm(m,l, phi_orig, theta_orig )
      
      # Then we calculate the values of the multipole moment at points 
      # obtained through coordinate reflection
      new  = [ sph_harm(m,l, phi_x   , theta_x    ), \
               sph_harm(m,l, phi_y   , theta_y    ), \
               sph_harm(m,l, phi_z   , theta_z    )]

      dic['ELL'] = '%d'%(l) # the Fortran arrays are 1-indexed
      dic['EMM'] = '%d'%(m) 
      
      restricted = [False, False]
      if(m == 0):
        restricted[1] = True
            
      for k in range(3):
        if(so.ReduceAxes[k] == 1):
          realpart = np.real(new[k])/np.real(orig)
          if(m != 0):
            impart   = np.imag(new[k])/np.imag(orig)
          else:
            impart   = 1    
          # If the multipole moment changes sign under coordinate reflection
          # and said coordinate reflection is part of the symmetries, then 
          # this multipole moment is constrained to zero and need not be
          # calculated. 
          if(abs(realpart + 1) < 0.1 ):
            restricted[0] = True
          if(abs(impart   + 1) < 0.1 ):
            restricted[1] = True
 
      for k in range(2):
        dic['IND']  = '%d'%k
        dic['REIM'] = REIM[k]
        if(restricted[k]):
          dic['ON']  = 0
        else:
          dic['ON']  = 1  
        # The code here is perfectly capable of writing ALL possible l,m 
        # combinations into the FORTRAN code. To keep things more readable, 
        # we only actually write multipole moments that should be 'switched on'.
        if(not restricted[k]):
          filling = filling + template_list.substitute(dic)
          filling = filling + template_comm.substitute(dic)
          
    filling = filling + '\n'
    
  dic={}
  dic['MAX_ELL']   = maxl 
  dic['FILL_LIST'] = filling
  dic['QUANT_AX']  = '%d'%(['X','Y','Z'].index(so.quant_axis) + 1)
  dic['SECOND_AX'] = so.second_axis
  
  if(so.timelike):
    dic['TR']  = ''
    dic['NTR'] = '!'
  else:
    dic['TR']  = '!'
    dic['NTR'] = ''

  if(fam_active):
    dic['FAM']  = 1
  else:
    dic['FAM']  = 0

  if (derivative_order >= 2 ):
    dic['CAN_DO_MAGNETIC'] = 1
  else:
    dic['CAN_DO_MAGNETIC'] = 0

  if(not dry_run):        
    substitute(src+fname, target+fname, dic)
        
def CartToSpher(x,y,z, so):
      """
        Translate Cartesian coordinates (x,y,z) into spherical coordinates
                    
                     (r, theta, phi)
                     
        There are six different ways to do this, depending on how we measure
        the angles with respect to a fixed set of Cartesian axes.
        
        In the "ordinary" way of doing things, 
              
                    r     = x^2 + y^2 + z^2
                    theta = arccos(z/r) 
                    phi   = atan2(y,x)

        The quantisation axis determines which Cartesian coordinates appears
        in the definition of theta, and the secondary axis 
        
        Input:
          x,y,z  : Cartesian coordinates 
          so     : Symmetry option object, which contains a choice for µ
                   quantisationaxis and secondary axis.
      """
  
      r = np.sqrt(x**2 + y**2 + z**2)
      
      if(so.quant_axis == 'X'):
        theta = np.arccos(x/r)
        if(so.second_axis == 1):
          phi = np.arctan2(z,y)
        else:
          phi = np.arctan2(y,z)
        
      elif(so.quant_axis == 'Y'):
        theta = np.arccos(y/r)
        if(so.second_axis == 1):
          phi = np.arctan2(z,x)
        else:
          phi = np.arctan2(x,z)

      elif(so.quant_axis == 'Z'):
        theta = np.arccos(z/r)
        if(so.second_axis == 1):
          phi = np.arctan2(y,x)
        else:
          phi = np.arctan2(x,y)
                  
      return r, theta, phi


def ProcessFission_MOI(fname, src, target, so, dry_run=False):
  """
    Generate the required Fortran code to process the fission_MOI.f90 file. 
    
    Input: 
      fname : filename of the Fortran template
      src   : directory of the templates
      target: directory to put the finished source code
      so    : SymmetryOption object, containing all the details on the 
              symmetries conserved during the calculation.
      dry_run    : Boolean, if True, do not write any files (default: False)
  """
  from src_heph.heph_substitute  import substitute

  dic = {}

  # Ugly Hack to check for parity conservation  
  symdic  = populatesymmetries()
  dic['PBROKEN']    = ' '
  dic['PCONSERVED'] = '!'
  for sym in so.generators:
    if (sym == symdic['P']):
      dic['PBROKEN'] = '!'
      dic['PCONSERVED'] = ''
  
  if(so.timelike):
    dic['TR']  = ''
    dic['NTR'] = '!'
  else:
    dic['TR']  = '!'
    dic['NTR'] = ''

  if(not dry_run):
    substitute(src+fname, target+fname, dic)
