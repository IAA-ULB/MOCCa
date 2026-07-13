#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
from string                    import Template
from src_heph.heph_symmetries  import *
from src_heph.heph_functional  import *

def ProcessPairing(fname, src, target, so, dry_run=False):
  """
   We process pairing.f90, depending on the symmetries imposed.
  """
  from src_heph.heph_substitute  import substitute

  dic = {}

  if(so.timelike):
      forbidBCS = '!'
      TR        = ''
      NTR       = '!'
  else:
      forbidBCS = ''
      TR        = '!'
      NTR       = ''

  dic['FORBIDBCS'] = forbidBCS
  dic['TR']        = TR
  dic['NTR']       = NTR
 
  # Checking if we are dealing with calls to subroutine vmicro or not
  # This reliance on a global definition on vmicro_found from the functional
  # module is ugly, but there is no obvious way around it due to the hierarchy
  # of the Fortran modules
  if(src_heph.heph_functional.vmicro_found):          
     dic['VMICRO']      = ' '
  else:
     dic['VMICRO']      = '!'
    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Ugly manual checking if parity is part of the generator set and 
  #  signalling this to the FORTRAN code
  symdic  = populatesymmetries()
  dic['PBROKEN'] = ' '

  for sym in so.generators:
    if (sym == symdic['P']):
     dic['PBROKEN'] = '!'

  if(not dry_run):
    substitute(src+fname, target+fname, dic)    
  
def ProcessHFB(fname, src, target, so, dry_run=False):
  """  
   We process HFB.f90, depending on the symmetries imposed and some choices
   made for the pairing functional.
  """

  from src_heph.heph_functional  import pairing_action_derorder
  from src_heph.heph_substitute  import substitute

  dic = {}

  if(so.timelike):
    dic["TR"]  = ' '
    dic["NTR"] = '!'
  else:
    dic["TR"]  = '!'
    dic["NTR"] = ' '

  if(pairing_action_derorder == 0):
    dic['N1DELTA']  = '!'
    dic['N2DELTA']  = '!'
    dic['N3DELTA']  = '!'
    dic['SYMDELTA'] = '!'
  elif(pairing_action_derorder == 1):
    dic['N1DELTA']  = ' '
    dic['N2DELTA']  = '!'
    dic['N3DELTA']  = '!'
    dic['SYMDELTA'] = ''
  elif(pairing_action_derorder == 2):
    dic['N1DELTA']  = ' '
    dic['N2DELTA']  = ' '
    dic['N3DELTA']  = '!'
    dic['SYMDELTA'] = ' '
  elif(pairing_action_derorder == 3):
    dic['N1DELTA']  = ' '
    dic['N2DELTA']  = ' '
    dic['N3DELTA']  = ' '
    dic['SYMDELTA'] = ' '

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Ugly manual checking if parity is part of the generator set and 
  #  signalling this to the FORTRAN code
  symdic  = populatesymmetries()
  dic['PCONSERVED'] = '!'
  dic['PBROKEN']    = ' '

  for sym in so.generators:
    if (sym == symdic['P']):
     dic['PBROKEN']    = '!'
     dic['PCONSERVED'] = ' '

  if(not dry_run):
    substitute(src+fname, target+fname, dic)    

def ProcessHartreeFock(fname, src, target, so, dry_run=False):
  """  
   We process HartreeFock.f90, depending on the symmetries imposed.
  """
  from src_heph.heph_substitute  import substitute
  
  dic = {}

  timelike = False  
  for g in so.generators:
    if(not g.linear and not g.hermitian):
      timelike = True

  if(timelike == 1):
    dic["TR"]  = ' '
    dic["NTR"] = '!'
  else:
    dic["TR"]  = '!'
    dic["NTR"] = ' '

  if(not dry_run):
    substitute(src+fname, target+fname, dic)
