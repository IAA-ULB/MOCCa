from string import Template


decl = Template(   tab + 'real(KIND=dp), allocatable :: $NAME(:$DECLIND,:)')
Ini  = Template(   tab + 'if(.not.allocated($NAME)) then \n'      \
                        + 2*tab + 'allocate($NAME(mv$ALLOCIND,2)) \n'      \
                        +   tab + 'endif  \n'                              \
                        +   tab + '$NAME = 0.0_dp  \n')
calc_a   = Template(   tab + '$NAME(:$ARG,:) = & \n')
calc_b   = Template( '$SIGN $DEN(:$IND,:) ')
