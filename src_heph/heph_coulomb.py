#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Module governing the Coulomb module of the FORTRAN code.
import src_heph.heph_functional 
from src_heph.heph_symmetries  import *
from string                    import Template


def ProcessCoulomb(fname, src, target, so, fam_active, dry_run=False):
  """
  Preprocess the Coulomb.f90 file.

  Input:
    fname:   name of the file to be processed
    src:     source directory of the file
    target:  target directory of the file
    so:      symmetry object for the current system
  
  Output:
    None, but writes the processed file to the target directory.
  """

  from src_heph.heph_substitute import substitute
  
  dic = {}
  
  axes = ['X', 'Y', 'Z']
  
  for k in range(3):
    if(so.ReduceAxes[k] == 1):
      dic['REDU%s'%axes[k]] = ' '
      dic['FULL%s'%axes[k]] = '!'

    else:
      dic['REDU%s'%axes[k]] = '!'
      dic['FULL%s'%axes[k]] = ' '

  if(fam_active):
   dic['FAM'] = '1'
  else:   
   dic['FAM'] = '0'       

  if(not dry_run):
    substitute(src + fname, target+fname, dic)  


