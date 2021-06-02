#-------------------------------------------------------------------------------
# | | | |  ___  _ __  | |__    __ _   ___  ___ | |_  ___   ___ 
# | |_| | / _ \| '_ \ | '_ \  / _` | / _ \/ __|| __|/ _ \ / __|
# |  _  ||  __/| |_) || | | || (_| ||  __/\__ \| |_| (_) |\__ \
# |_| |_| \___|| .__/ |_| |_| \__,_| \___||___/ \__|\___/ |___/
#              |_|                                             
#-------------------------------------------------------------------------------
# Module governing the wavefunctions module of the FORTRAN code.
import src_heph.heph_functional 
from src_heph.heph_symmetries  import *
from string                    import Template


def ProcessWavefunctions(fname, src, target, so):
    """    
      Process the wavefunctions.f90 file of the Fortran code.       
        
      A 3D calculation, respecting some subgroup of the D^T(D)_2h groups, 
      will always have single-particle wavefunctions that can be divided into 
      at most eight blocks according to their behaviour under symmetries. 

      2 for protons versus neutrons
      2 for a conserved hermitian, linear symmetry (parity in EV8)
      2 for a conserved antihermitian, linear symmetry (z-signature in EV8) 

      A conserved antilinear, antihermitian symmetry will then allow for the 
      elimination of one of these sets in a practical calculation.

     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

      N2/N3 :  decides which derivative routines to comment/uncomment 
               depending on the order of derivatives in the functional

    """

#    self.syms       = syms
#    self.combs      = combs

#    print (multiply_quantum_numbers(so.syms[0].permutation, so.combs[0], 1 ))
#    print (multiply_quantum_numbers(so.syms[0].permutation, so.combs[0], 2 ))
#    print (multiply_quantum_numbers(so.syms[0].permutation, so.combs[0], 3 ))
#    print (multiply_quantum_numbers(so.syms[0].permutation, so.combs[0], 4 ))
#    print ()
#    print (multiply_quantum_numbers(so.syms[1].permutation, so.combs[1], 1 ))
#    print (multiply_quantum_numbers(so.syms[1].permutation, so.combs[1], 2 ))
#    print (multiply_quantum_numbers(so.syms[1].permutation, so.combs[1], 3 ))
#    print (multiply_quantum_numbers(so.syms[1].permutation, so.combs[1], 4 ))
#    exit()

    dic={}
    dic['N2'] = '!'
    dic['N3'] = ' '

    if(src_heph.heph_functional.derivative_order == 1):
        dic['N2'] = ' '    
        dic['N3'] = '!'
    elif(src_heph.heph_functional.derivative_order == 2): 
        dic['N2'] = ' '
        dic['N3'] = '!'
    elif(src_heph.heph_functional.derivative_order == 3):
        dic['N2'] = '!'
        dic['N3'] = ' '
    
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic)) 


    
    
