"""
    This module implements functions to shortcut the typical linear 
    combinations of Skyrme parameters that occur in coupling constants 
    of EDF models.
"""

def skyrme_cc(t,x,p,s,it):
    """
    Computes the quantity C^+/-_cST, i.e. a useful intermediate object 
    that occurs repeatedly in the construction of the coupling constants 
    of *CENTRAL* Skyrme terms.
    
    These follow the conventions of Equation B1 in appendix B of 
        W. Ryssens and M. Bender, Physical Review C 104, 044308 (2021).

    In table format:

      C^+/-_cST       +              -
                  t      tx      t    tx
           
        c00     +3/8     0     +5/8  +1/2 
        c01     -1/8   -1/4    +1/8  +1/4
        c10     -1/8   +1/4    +1/8  +1/4
        c11     -1/8     0     +1/8   0

    Args: 
        t:  real, parameter 
        x:  real, parameter   
        p:  integer, parameter, +1/-1 for C^+/C^-
        s:  integer, spin     , 0/1 for S=0/1
        it: integer, isospin  , 0/1 for T=0/1
    Returns:
        c:  real, C^+/-_ST
    
    """
    c = 0
    match p:
        case 1: # C^+_c
            match s:
                case 0: # C^+_c0
                    match it:
                        case 0:
                            c = +3.0/8.0 * t                     # C^+_c00
                        case 1:
                            c = -1.0/8.0 * t  - 1.0/4.0 * t * x  # C^+_c01
                case 1: # C^+_c1
                    match it:
                        case 0:
                            c = -1.0/8.0 * t  + 1.0/4.0 * t * x  # C^+_c10
                        case 1:
                            c = -1.0/8.0 * t                     # C^+_c11

        case -1: # C^-_c
            match s:
                case 0: # C^-_c0
                    match it:
                        case 0:
                            c = +5.0/8.0 * t  + 1.0/2.0 * t * x  # C^-_c00
                        case 1:
                            c = +1.0/8.0 * t  + 1.0/4.0 * t * x  # C^-_c01
                case 1: # C^-_c1
                    match it:
                        case 0:
                            c = +1.0/8.0 * t  + 1.0/4.0 * t * x  # C^-_c10
                        case 1:
                            c = +1.0/8.0 * t                     # C^-_c11
    return c

def skyrme_ct(t,p,it):
    """
    Computes the quantity C^+/-_ST, i.e. a useful intermediate object
    that occurs repeatedly in the construction of the coupling constants
    of *TENSOR* and *spin-orbit* Skyrme terms.

    To the best of my knowledge, these have not been published so far; 
    they do correspond to the conventions of M. Bender.

    In table format:

      C^+/-_tST    +       -
                   t       t        
        c10     +1/4    +3/4
        c11     -1/4    +1/4

    Args:
        t:  real, parameter 
        p:  integer, parameter, +1/-1 for C^+/C^-
        it: integer, isospin  , 0/1 for T=0/1

    Returns:
        c:  real, C^+/-_ST
    """
    c = 0

    match p:
        case  1: # C^+_t
            match it:
                case 0:
                    c =  1.0/4.0 * t  # C^+_t10
                case 1:
                    c = -1.0/4.0 * t  # C^+_t11
        case -1: # C^-_t
            match it:
                case 0:
                    c =  3.0/4.0 * t  # C^-_t10
                case 1:
                    c =  1.0/4.0 * t  # C^-_t11
    return c
