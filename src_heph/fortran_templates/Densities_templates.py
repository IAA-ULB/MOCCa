from string import Template as T
tab = '    '
#-------------------------------------------------------------------------------
# Density calculation template to fill in

Den_1    = T( 2*tab+'$NAME(i$IND,it) = $NAME(i$IND,it) + $WEIGHT * (')
Den_diag = T( \
                    tab+'$SIGN $LEFTWF(i$LIND,$LCOMP,$LEFTWAVE) * ' + \
                        '$RIGHTWF(i$RIND,$RCOMP,$RIGHTWAVE)')

Ini      = T(   tab+'if(.not.allocated($NAME)) then     \n' + \
                           2*tab+'allocate($NAME(mv$DIM,$ISOSIZE)) \n'   + \
                           2*tab+'$NAME = 0.0d0 \n'                      + \
                             tab+'endif \n'                              ) 
Zero_template  = T(   tab+'$NAME = 0.0d0 \n')
Clean_template = T(   tab+'if(allocated($NAME)) deallocate($NAME)')

Dec       = T(   tab + 'real*8, allocatable, target :: $NAME(:$TOTALIND,:)')
Der_indep = T( 2*tab + \
             'call Derive_$DIR($NAME(:$IND,it), $PS,der_$NAME(:$DERIND,it)) \n') 
Lap   = T( 2*tab + \
       'call Derive_lap ($NAME(:$IND,it), $PX,$PY,$PZ, lap_$NAME(:$IND,it)) \n')

Der_sum = T( 2*tab + 'der_$NAME(:$DERIND,it) = $LEFTDEN(:$DERLIND,it) + $RIGHTDEN(:$DERRIND,it) \n') 
Der_der_sum_a = T( 2*tab + 'Der_$NAME(:$DERIND,it) = & \n') 
Der_der_sum_b = T( 3*tab + '&  $LEFTDEN(:$DERLIND,it) + 2*$CENTRALDEN(:$DERCIND,it) + $RIGHTDEN(:$DERRIND,it) \n') 

Lap_sum_a = T( 2*tab + 'lap_$NAME(:$IND,it) = & \n') 
Lap_sum_b = T( 3*tab + '& +  $LEFTDEN(:$DERLIND,it) + 2*$CENTRALDEN(:$DERCIND,it) + $RIGHTDEN(:$DERRIND,it) &\n') 


iso_normal   =  T( tab +'$NAME(:$IND,3) = $NAME(:$IND,1) + $NAME(:$IND,2) \n' \
                  +tab +'$NAME(:$IND,4) = $NAME(:$IND,1) - $NAME(:$IND,2) \n' )
iso_der      =  T( tab +'Der_$NAME(:$DERIND,3) = Der_$NAME(:$DERIND,1) + Der_$NAME(:$DERIND,2) \n' \
                  +tab +'Der_$NAME(:$DERIND,4) = Der_$NAME(:$DERIND,1) - Der_$NAME(:$DERIND,2) \n' )
iso_lap      =  T( tab +'Lap_$NAME(:$IND,3) = Lap_$NAME(:$IND,1) + Lap_$NAME(:$IND,2) \n' \
                  +tab +'Lap_$NAME(:$IND,4) = Lap_$NAME(:$IND,1) - Lap_$NAME(:$IND,2) \n' )

mpi          =  T( tab + 'call MPI_ALLREDUCE(MPI_IN_PLACE,$NAME(:$TOTALIND,1:2),$TRANS_SIZE, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpi_err)')

#-------------------------------------------------------------------------------
# Some templates for comments to put into the densities file
Den_comment          = T(2*tab+'! Calculation of density $NAME \n')
Den_iso_comment      = T( tab+'! Isospin representation of density $NAME(:$IND,:) and its derivatives \n')

Den_comment_deriv    = T(2*tab+ \
                                   '! Derivation of density $NAME(:$IND,it) \n')
Den_comment_deriv_b  = (2*tab+'! LAP = %d, DER = %d \n')
#                      Note that this is not a TEMPLATE string for historical 
#                      reasons that I'm too lazy to fix
Den_line             = T(2*tab+ \
  '! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  \n')

