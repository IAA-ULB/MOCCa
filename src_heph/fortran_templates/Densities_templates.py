from string import Template as T
tab = '    '
#-------------------------------------------------------------------------------
# Density calculation template to fill in

Den_1    = T( 2*tab+'$NAME(i$IND,it) = $NAME(i$IND,it) + $WEIGHT * (')
Den_diag = T( \
                    tab+'$SIGN $LEFTWF(i$LIND,$LCOMP,$LEFTWAVE) * ' + \
                        '$RIGHTWF(i$RIND,$RCOMP,$RIGHTWAVE)')

Ini      = T(   tab+'if(.not.allocated($NAME)) then     \n' + \
                           2*tab+'allocate($NAME(mv$DIM,2)) \n'          + \
                           2*tab+'$NAME = 0.0d0 \n'                      + \
                             tab+'endif \n'                              ) 
Zero_template  = T(   tab+'$NAME = 0.0d0 \n')
Clean_template = T(   tab+'if(allocated($NAME)) deallocate($NAME)')

Dec   = T( \
                     tab + 'real*8, allocatable, target :: $NAME(:$TOTALIND,:)')
Der_indep = T( 2*tab + \
             'call Derive_$DIR($NAME(:$IND,it), $PS,der_$NAME(:$DERIND,it)) \n') 
Lap   = T( 2*tab + \
       'call Derive_lap ($NAME(:$IND,it), $PX,$PY,$PZ, lap_$NAME(:$IND,it)) \n')

#-------------------------------------------------------------------------------
# Some templates for comments to put into the densities file
Den_comment          = T(2*tab+'! Calculation of density $NAME \n')

Den_comment_deriv    = T(2*tab+ \
                                   '! Derivation of density $NAME(:$IND,it) \n')
Den_comment_deriv_b  =         (2*tab+'! LAP = %d, DER = %d \n')
Den_line             = T(2*tab+ \
  '! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  \n')

