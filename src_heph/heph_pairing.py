#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
from string          import Template

def ProcessPairing(fname, src, target, so):
  """
   We process pairing.f90, depending on the symmetries imposed.
  """

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

  with open(src+fname, 'r') as template:
    with open(target+fname, 'w') as generated:
        for line in template:
            generated.write(Template(line).substitute(dic))  


def ProcessHFB(fname, src, target, so):
  """  
   We process HFB.f90, depending on the symmetries imposed.
  """

  dic = {}

  if(so.timelike):
    dic["TR"]  = ' '
    dic["NTR"] = '!'
  else:
    dic["TR"]  = '!'
    dic["NTR"] = ' '

  with open(src+fname, 'r') as template:
    with open(target+fname, 'w') as generated:
        for line in template:
            generated.write(Template(line).substitute(dic))  


def ProcessHartreeFock(fname, src, target, so):
  """  
   We process HartreeFock.f90, depending on the symmetries imposed.
  """

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

  with open(src+fname, 'r') as template:
    with open(target+fname, 'w') as generated:
        for line in template:
            generated.write(Template(line).substitute(dic))  


# Note: the BCS module needs no processing at the moment.
# def ProcessBCS(fname, src, target, so): 
