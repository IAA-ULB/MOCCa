#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Module governing the wavefunctions module of the FORTRAN code.
import heph_functional
from string           import Template


def ProcessWavefunctions(fname, src, target):
    
    dic={}
    dic['BLOCKS'] = 8
    dic['N2'] = '!'
    dic['N3'] = ' '
    

    if(heph_functional.derivative_order == 1):
        dic['N2'] = ' '    
        dic['N3'] = '!'
    elif(heph_functional.derivative_order == 2): 
        dic['N2'] = ' '
        dic['N3'] = '!'
    elif(heph_functional.derivative_order == 3):
        dic['N2'] = '!'
        dic['N3'] = ' '
    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic)) 


    
    
