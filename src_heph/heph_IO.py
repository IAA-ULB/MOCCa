#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
from string          import Template
from src_heph.heph_symmetries import symmetryencoding

def ProcessIO(fname, src, target, so, oldso):
    """
      
    """
    
    dic = {}

    SYM_CODE   = symmetryencoding(so)
    TRANS_CODE = symmetryencoding(oldso)

    dic['SYM_CODE'] = SYM_CODE
    dic['TRANS_CODE'] = TRANS_CODE

    with open(src+fname, 'r') as template:
      with open(target+fname, 'w') as generated:
        for line in template:
            generated.write(Template(line).substitute(dic))  

