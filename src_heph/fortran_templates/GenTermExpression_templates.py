#-------------------------------------------------------------------------------
#  This python module contains a set of Templates that are used by the
#              
#    GenTermExpression function in the functional module of Hephaestos
#  
#  to generate code for every individual term of the EDF.
#
#-------------------------------------------------------------------------------
from string                  import Template as T
tab = '   '


decl   = T( tab + 'real(KIND=dp) :: $CPCTE(2,2), $TERM(2,2)')
edent  = T('sum($DEN(:$IND,:),2)')
edenq  = T('$DEN(:$IND,$IT)')

comment= T(   tab + '!' + 38 * '- ' + '\n' +               \
                            tab + '! Calculation of $TERM \n')
                            
end_comment     =             tab + '!' + 38 * '- ' + '\n'
                            
doloop    =    tab + 'do %s = 1, 3 \n'
enddoloop =    tab + 'enddo \n'

sumtotal  = T( 2*tab + ' & + $TERM(:,1)')
pairtotal  = T( 2*tab + ' & + $TERM(:,1)')

calc_z = T(   tab + 'Edensity = 0.0_dp \n')
calc_a = T(   tab + 'EDensity(:,3) = Edensity(:,3) $SIGN $EDENT\n')
calc_b = T(   tab + 'Edensity(:,1) = Edensity(:,1) $SIGN $EDENN\n' + \
                     tab + 'Edensity(:,2) = Edensity(:,2) $SIGN $EDENP\n'  )
                     
                     
calc_DD= T(   tab + 'do m=1,3 \n' +                                 \
                   2*tab + 'EDensity(:,m) = Edensity(:,m) * $DD \n' +      \
                     tab + 'enddo \n')
calc_c = T(   tab + '$TERM(1,2) = $CPCTE(1,2) * sum( Edensity(:,3)) * dv \n')
calc_d = T(   tab + '$TERM(2,2) = $CPCTE(2,2) * dv & \n'+           \
              tab + '&'+6*tab+' * sum(Edensity(:,1) + Edensity(:,2) ) \n')

# This new division of contributions is only correct for bilinear terms.
# If we want to include trilinear terms, we should probably explicitly 
# construct isoscalar and isovector densities.    
calc_e = T(   tab + '$TERM(1,1) = dv * ( & \n' + 
                            tab + '& ($CPCTE(1,2)+0.5*$CPCTE(2,2)) * sum(Edensity(:,3))) \n')
calc_f = T(   tab + '$TERM(2,1) = dv * & \n' +
                            tab + '& (         $CPCTE(2,2) * sum(Edensity(:,1) + Edensity(:,2)) & \n' + 
                            tab + '&   - 0.5 * $CPCTE(2,2) * sum(Edensity(:,3))               )   \n')

calc_pair_a= T(   tab + 'Edensity(:,1) = Edensity(:,1) + $EDENN\n' + \
                            tab + 'Edensity(:,2) = Edensity(:,2) + $EDENP\n'  )
calc_pair_c = T(  tab + '$TERM(1,1) = $CPCTE(1,1) * sum( Edensity(:,1)) * dv \n')
calc_pair_d = T(  tab + '$TERM(2,1) = $CPCTE(2,1) * sum( Edensity(:,2)) * dv \n')

calc_coef  = T(tab + '$CPCTE(1,1) = $EXP1 \n' + \
                      tab + '$CPCTE(2,1) = $EXP2 \n' + \
                      tab + '$CPCTE(1,2) = $CPCTE(1,1) - $CPCTE(2,1) \n' + \
                      tab + '$CPCTE(2,2) =             2*$CPCTE(2,1) \n')       

calc_coef_pair = T(tab + '$CPCTE(1,1) = $EXP1 \n' + \
                             tab + '$CPCTE(2,1) = $EXP2 \n')    
                          
write_edensity = T( tab + ' call output_Edensity(Edensity, "$FILENAME")' )
print          = T(tab +" print('(2x, a28, 3f15.6)'), rps('$TERM',28),    & \n"+ 
                        tab +"                           $TERM(:,1),& \n"+
                        tab +"                       sum($TERM(:,1))  \n")
print_P  = T(tab +" print('(2x, a28, 30x, f15.6)'), rps('$TERM',28),& \n"+
                        tab +"                       sum($TERM(:,1))  \n")


print_cpl_pn   = T(tab +" print('(a30 , 4f15.6)'), '$CPCTE', $CPCTE(:,1)")
print_cpl_iso  = T(tab +" print('(a30 , 4f15.6)'), '$CPCTE', $CPCTE(:,2)")
print_cpl_pair = T(tab +" print('(a30 , 4f15.6)'), '$CPCTE', $CPCTE(:,1)")

rear      = T(tab +" e_rear = e_rear $REARCOEF*sum($TERM(:,2))\n")
rear_p    = T(tab +" e_rear = e_rear $REARCOEF*sum($TERM(:,1))\n")
