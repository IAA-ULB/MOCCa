#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
from string          import Template

def ProcessHFB(fname, src, target, generators):
  """  
   We process HFB.f90, depending on the symmetries imposed.
  """

  dic = {}

  TimeReversal = 1
  if(TimeReversal == 1):
    dic["TR"]  = ' '
    dic["NTR"] = '!'
  else:
    dic["TR"]  = '!'
    dic["NTR"] = ' '

  with open(src+fname, 'r') as template:
    with open(target+fname, 'w') as generated:
        for line in template:
            generated.write(Template(line).substitute(dic))  
