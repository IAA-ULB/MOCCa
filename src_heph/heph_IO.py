#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
from string          import Template
from src_heph.heph_symmetries import symmetryencoding
from src_heph.heph_functional import Densities_needed

def ProcessIO(fname, src, target, so, oldso, fam_active):
    """
    Preprocess the IO.f90 file.
    
    Input:
      fname:       name of the file to be processed
      src:         source directory of the file
      target:      target directory of the file
      so:          symmetry object for the current system
      oldso:       symmetry object characterising a previous mf run
      fam_active:  whether we are compiling a mean-field code (False)
                   or a finite amplitude linear response code (True)
    """
    
    dic = {}

    SYM_CODE   = symmetryencoding(so)
    TRANS_CODE = symmetryencoding(oldso)

    dic['SYM_CODE'] = SYM_CODE
    dic['TRANS_CODE'] = TRANS_CODE
    
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
          
    if(so.timelike):
      dic['TR'] = ''
      dic['NTR']= '!'
    else:
      dic['TR'] = '!'
      dic['NTR']= ''

    if('D_Nm_Nm' not in Densities_needed):
      dic['TAUSCALAR'] = '!'
      dic['TAUTENSOR'] = ' '
    else:
      dic['TAUSCALAR'] = ' '
      dic['TAUTENSOR'] = '!'
      
    if(fam_active):
      dic['FAM'] = '!'
      dic['MF']  = ''
    else:
      dic['FAM'] = ''
      dic['MF']  = '!'

    with open(src+fname, 'r') as template:
      with open(target+fname, 'w') as generated:
        for line in template:
            generated.write(Template(line).substitute(dic))  

