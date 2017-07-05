#============================================================================================
#          _______  _______           _______  _______  _______ _________ _______  _______ 
#|\     /|(  ____ \(  ____ )|\     /|(  ___  )(  ____ \(  ____ \\__   __/(  ___  )(  ____ \
#| )   ( || (    \/| (    )|| )   ( || (   ) || (    \/| (    \/   ) (   | (   ) || (    \/
#| (___) || (__    | (____)|| (___) || (___) || (__    | (_____    | |   | |   | || (_____ 
#|  ___  ||  __)   |  _____)|  ___  ||  ___  ||  __)   (_____  )   | |   | |   | |(_____  )
#| (   ) || (      | (      | (   ) || (   ) || (            ) |   | |   | |   | |      ) |
#| )   ( || (____/\| )      | )   ( || )   ( || (____/\/\____) |   | |   | (___) |/\____) |
#|/     \|(_______/|/       |/     \||/     \|(_______/\_______)   )_(   (_______)\_______)
#============================================================================================
# Director of preprocessing of the FORTRAN files: knows which file gets which make-over. 
# 
# Currently linked: 
#
#       densities.f90 <=> heph_densities.py
#
#============================================================================================

from string           import Template
from heph_symmetries  import ReduceAxes
import heph_functional

def ProcessDerivatives(fname, src, target):
    #===========================================================================
    # LINESIZEX/Y/Z
    # these determine the size of the box in that direction
    # as a function of nx,ny and nz.
    #
    # DERX/Y/Z_ONE  
    # DERX/Y/Z_TWO
    #
    # Determines the derivative matrices as a function of variables C and D
    # in the FORTRAN file. 
    #      C = coefficients of the line-segment itself
    #      D = coefficients of the symmetry-partner of the line-segment
    #
    # For every Cartesian direction, the viable options are thus
    #    
    # 1) Symmetry reduced, symmetry partner is the line itself   
    #    DERX/Y/Z_ONE = C + D
    #    DERX/Y/Z_TWO = C - D
    # 2) Symmetry partner is not the line itself
    #    DERX/Y/Z_ONE = C
    #    DERX/Y/Z_TWO = D
    # 3) Not symmetry reduced
    #    DERX/Y/Z_ONE = C
    #    DERX/Y/Z_TWO = 0 ( The D-variable is no longer relevant) 
    #
    # 
    # DERSYMX/Y/Z
    #  Variables that allow the commenting out of the inclusion of 
    #  symmetry partners for every direction. There are two options:
    #
    # 1) Symmetry reduced, symmetry partner is the line itself 
    #    OR Not symmetry reduced
    #    $DERSYMX/Y/Z = '!'
    # 2) Symmetry reduced, symmetry partner is not the line itself
    #    $DERSYMX/Y/Z = ''
    #
    #
    # SYMPARTNERX/Y/Z
    #  Variable indicating the symmetric partner of the line-segment
    #  as a function of (i,j,k) which are the (x,y,z) indices on 
    #  the mesh. Only relevant when there is a symmetry present
    #  that does not directly point back to the same line. 
    # 
    #===========================================================================
    # Something that is rather determined by the functional is the maximum order
    # of derivatives that is needed. This is determined by the following strings
    #
    # $N2DIAG  => 1st order derivatives and the diagonal 2nd order ones
    # $N2ALL   => 1st order derivatives and all of the 2nd order ones
    # $N3ALL   => up to and including 3rd order derivatives
    #===========================================================================
    dic={}

    dic['DERSYMX'] = '!'
    dic['DERSYMY'] = '!' 
    dic['DERSYMZ'] = '!'

    if(ReduceAxes[0] == 1):
       dic['LINESIZEX']    = '2*nx'
       dic['DERX_ONE']     = 'C + D'
       dic['DERX_TWO']     = 'C - D'
       dic['SYMPARTNERX']  = ''  

    else:
       dic['LINESIZEX']    = 'nx'
       dic['DERX_ONE']     = 'C'
       dic['DERX_TWO']     = '0'  
       dic['SYMPARTNERX']  = ''          

    if(ReduceAxes[1] == 1):
       dic['LINESIZEY']    = '2*ny'
       dic['DERY_ONE']     = 'C + D'
       dic['DERY_TWO']     = 'C - D'
       dic['SYMPARTNERY']  = ''  
    else:
       dic['LINESIZEY']    = 'ny'
       dic['DERY_ONE']     = 'C'
       dic['DERY_TWO']     = '0'
       dic['SYMPARTNERY']  = ''  
    
    if(ReduceAxes[2] == 1):
       dic['LINESIZEZ']    = '2*nz'
       dic['DERZ_ONE']     = 'C + D'
       dic['DERZ_TWO']     = 'C - D'
       dic['SYMPARTNERZ']  = ''  
    else:
       dic['LINESIZEZ']    = 'nz'
       dic['DERZ_ONE']     = 'C'
       dic['DERZ_TWO']     = '0'
       dic['SYMPARTNERZ']  = ''  
    
    dic['DY'] = 'D'
    dic['DZ'] = 'D'
        
    if(heph_functional.derivative_order == 1):    
        dic['N2DIAG'] = ' '
        dic['N2ALL']  = '!'  
        dic['N3ALL']  = '!'   
    elif(heph_functional.derivative_order == 2):    
        dic['N2DIAG'] = '!'
        dic['N2ALL']  = ' '   
        dic['N3ALL']  = '!'   
    elif(heph_functional.derivative_order == 3):    
        dic['N2DIAG'] = '!'
        dic['N2ALL']  = '!'   
        dic['N3ALL']  = ' ' 
        
    with open(src+fname, 'r') as template:
        with open(target+fname, 'w') as generated:
            for line in template:
                generated.write(Template(line).substitute(dic))
