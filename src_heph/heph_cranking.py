#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
#
#===============================================================================

from string          import Template
from src_heph.heph_functional import Densities_needed
#===============================================================================
def ProcessCranking(fname, src, target, so, dry_run=False):
  """
    Process the cranking.f90 file. 

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  """
  from src_heph.heph_substitute  import substitute

  directions = []

  dic = {}

  dic['NTR'] = '!'
  dic['TR']  = ''
  if( not so.timelike):
    directions = [3]
    dic['NTR'] = ''
    dic['TR']  = '!'

  dic['CRANKDIR'] = ""
  for d in directions:
    dic['CRANKDIR'] = dic['CRANKDIR'] + "%d,"%d
    
  dic['CRANKLEN'] = len(directions)

  if('D_Nm_Nm' not in Densities_needed and 'D_N_N' not in Densities_needed):
    dic['TAUPRESENT'] = 0
  else:
    dic['TAUPRESENT'] = 1

  if(not dry_run):
    substitute(src+fname, target+fname, dic)