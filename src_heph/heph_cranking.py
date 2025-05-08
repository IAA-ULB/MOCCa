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

def ProcessCranking(fname, src, target, so):
  """
    Process the cranking.f90 file. 

    - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

  """
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

  with open(src+fname, 'r') as template:
    with open(target+fname, 'w') as generated:
        for line in template:
            generated.write(Template(line).substitute(dic))
