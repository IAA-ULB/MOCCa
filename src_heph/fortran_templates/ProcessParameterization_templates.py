#-------------------------------------------------------------------------------
#  This python module contains a set of Templates that are used by the
#              
#    ProcessParameterization function in the functional module of Hephaestos
#  
#  to generate code for the reading of all parameters of a parameterization
#  from file.
#-------------------------------------------------------------------------------
from string                  import Template as T
tab = '    '

decl_real    = T(  tab + 'real(KIND=dp) :: $PARAM = -123456789 \n')
decl_int     = T(  tab + 'integer       :: $PARAM = -123456789 \n')


read    = T(  tab + '& $PARAM, & \n')
print_real   = T(2*tab + 'print "(a6,2x, f10.3)", "$PARAM", $PARAM \n ')
print_int    = T(2*tab + 'print "(a6,2x,    i5)", "$PARAM", $PARAM \n ')

reset   = T(  tab + '$PARAM = -123456789 \n')

check_a = T(  tab + 'if($PARAM .eq. -123456789) then \n')
check_b = T(2*tab + '   print *, "$PARAM not read from .param file." \n')
check_c = T(2*tab + '   stop \n')        
check_d = T(  tab + 'endif \n')
