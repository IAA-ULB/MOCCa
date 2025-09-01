#-------------------------------------------------------------------------------
#  This python module contains a set of Templates that are used by the
#  "GenTermExpression" function in the functional module of Hephaestos
#  to generate code for every individual term of the EDF.
#
#-------------------------------------------------------------------------------
from string                  import Template as T
tab = '    '

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Declaration statements
decl   = T( tab + 'real(KIND=dp) :: $TERM$GROUPNUMBER')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statements for the calculation of terms 
comment         = T(   tab + '!' + 38 * '- ' + '\n' +               \
                            tab + '! Calculation of $TERM \n')
end_comment     =             tab + '!' + 38 * '- ' + '\n' 

edent     = T('(R%$DEN(:$IND,$ISO))')
edent_DD  = T('(pow(R%$DEN(:$IND,$ISO), $EXP))')

edent_fam    = T('DBLE(R%$DEN(:$IND,$ISO))')
edent_DD_fam = T('DBLE(pow(R%$DEN(:$IND,$ISO), $EXP))')

calc_z    = T(   tab + 'Edensity = 0.0_dp \n')
calc_a    = T(   tab + 'EDensity = Edensity $SIGN $EDENT\n')
calc_extra= T(   tab + 'EDensity = Edensity * $EXTRA\n')
calc_b    = T(   tab + '$TERM$GROUPINDEX = $CPCTE * sum( Edensity) * dv \n')

doloop    =    tab + 'do %s = 1, 3 \n'
enddoloop =    tab + 'enddo \n'

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculation of the coupling constant
calc_coef  = T(tab + '! Coupling constant of $TERM, $ISO \n' + \
               tab + 'coupl_constant($NCONSTANT) = $EXP \n')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Printing of the value of the energy
print    = T(tab +" print $FMT,                           rps('$TERM',38), & \n"+
                        tab +"                           $ISO_FULL,       & \n"+ 
                        tab +"                           $TERM$GROUPINDEX \n")
#print_P  = T(tab +" print('(2x, a28, 30x, f15.6)'), rps('$TERM',28),& \n"+
#                        tab +"                       sum($TERM(:,1))  \n")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Small statement to have this term be added in different summations of energies
sumtotal   = T( 2*tab + ' & + $TERM$GROUPINDEX ')
pairtotal  = T( 2*tab + ' & + $TERM$GROUPINDEX ')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Contributions to the rearrangement energy
rear      = T(tab +"e_rear = e_rear $REARCOEF*$TERM$GROUPINDEX\n")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statements for printing the coupling constants
print_cpl_ph   = T(tab +" print $FMT, rps('$TERM',38), $ISO_FULL, $NCONSTANT, $GROUPINDEX, coupl_constant($NCONSTANT)")
print_cpl_pair = T(tab +" print $FMT, rps('$TERM',38), $ISO_FULL, $NCONSTANT, $GROUPINDEX, coupl_constant($NCONSTANT)")

#-------------------------------------------------------------------------------
