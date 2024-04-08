from string import Template as T
tab = '    '
#-------------------------------------------------------------------------------
# Density calculation template to fill in

Den_1    = T( 2*tab+'R%$NAME(i$IND,it) = R%$NAME(i$IND,it) + $WEIGHT * (')
Den_diag = T( \
                    tab+'$SIGN $LEFTWF(i$LIND,$LCOMP,$LEFTWAVE) * ' + \
                        '$RIGHTWF(i$RIND,$RCOMP,$RIGHTWAVE)')

Ini      = T(   tab+'if(.not.allocated(R%$NAME)) then     \n' + \
                           2*tab+'allocate(R%$NAME(mv$DIM,$ISOSIZE)) \n'   + \
                           2*tab+'R%$NAME = 0.0d0 \n'                      + \
                             tab+'endif \n'                              ) 
Zero_template  = T(   tab+'R%$NAME = 0.0d0 \n')
Clean_template = T(   tab+'if(allocated(R%$NAME)) deallocate(R%$NAME)')

Dec       = T(   tab + 'real*8, allocatable :: $NAME(:$TOTALIND,:)')
Der_indep = T( 2*tab + \
             'call Derive_$DIR(R%$NAME(:$IND,it),$PS,R%der_$NAME(:$DERIND,it)) \n') 
Lap   = T( 2*tab + \
       'call Derive_lap(R%$NAME(:$IND,it), $PX,$PY,$PZ,R%lap_$NAME(:$IND,it)) \n')

Der_sum = T( 2*tab + 'R%der_$NAME(:$DERIND,it) = R%$LEFTDEN(:$DERLIND,it) + R%$RIGHTDEN(:$DERRIND,it) \n') 
Der_der_sum_a = T( 2*tab + 'R%Der_$NAME(:$DERIND,it) = & \n') 
Der_der_sum_b = T( 3*tab + '&  R%$LEFTDEN(:$DERLIND,it) + 2*R%$CENTRALDEN(:$DERCIND,it) + R%$RIGHTDEN(:$DERRIND,it) \n') 

Lap_sum_a = T( 2*tab + 'R%lap_$NAME(:$IND,it) = & \n') 
Lap_sum_b = T( 3*tab + '& +  R%$LEFTDEN(:$DERLIND,it) + 2*R%$CENTRALDEN(:$DERCIND,it) + R%$RIGHTDEN(:$DERRIND,it) &\n') 


iso_normal   =  T( tab +'R%$NAME(:$IND,3) = R%$NAME(:$IND,1) + R%$NAME(:$IND,2) \n' \
                  +tab +'R%$NAME(:$IND,4) = R%$NAME(:$IND,1) - R%$NAME(:$IND,2) \n' )
iso_der      =  T( tab +'R%Der_$NAME(:$DERIND,3) = R%Der_$NAME(:$DERIND,1) + R%Der_$NAME(:$DERIND,2) \n' \
                  +tab +'R%Der_$NAME(:$DERIND,4) = R%Der_$NAME(:$DERIND,1) - R%Der_$NAME(:$DERIND,2) \n' )
iso_lap      =  T( tab +'R%Lap_$NAME(:$IND,3) = R%Lap_$NAME(:$IND,1) + R%Lap_$NAME(:$IND,2) \n' \
                  +tab +'R%Lap_$NAME(:$IND,4) = R%Lap_$NAME(:$IND,1) - R%Lap_$NAME(:$IND,2) \n' )

mpi          =  T( tab + 'call MPI_ALLREDUCE(MPI_IN_PLACE,R%$NAME(:$TOTALIND,1:2),$TRANS_SIZE, MPI_REAL8, MPI_SUM, MPI_COMM_WORLD, mpi_err)')

Add_template        = T( tab + ' R%$NAME = R1%$NAME + R2%$NAME')
Multiply_template   = T( tab + ' R%$NAME = a * R1%$NAME')

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

