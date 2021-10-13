#
#
#
#
#-------------------------------------------------------------------------------

from string import Template as T    
tab = '    '

field_decl = T(  tab + 'real(KIND=dp), allocatable, target :: $FIELD(:$DECLIND,:) \n')
fhist_decl = T(  tab + 'real(KIND=dp), allocatable :: ${FIELD}_hist(:$DECLIND,:) \n')

field_allo = T(2*tab + 'if(.not.allocated($FIELD)) then\n'+\
               3*tab + 'allocate($FIELD(mv$ALLOCIND,2)) \n'   +\
               3*tab + 'allocate(${FIELD}_hist(mv$ALLOCIND,2)) \n'+\
               3*tab + '$FIELD = 0.0 ; ${FIELD}_hist = 0.0 \n' + \
               2*tab + 'endif \n')

field_allo_b = T(3*tab + 'if(.not.allocated($FIELD)) then\n'+\
                     4*tab + 'allocate($FIELD(mv$ALLOCIND,2)) \n'   +\
                     4*tab + 'allocate(${FIELD}_hist(filemv$ALLOCIND,2)) \n'+\
                     4*tab + '$FIELD = 0.0 ; ${FIELD}_hist = 0.0 \n' + \
                     3*tab + 'endif \n')

field_hist   = T( 2*tab + 'if(calcall) then \n' + 
                          3*tab + '${FIELD}_hist = $FIELD \n' + 
                          3*tab + '$FIELD = 0.0 \n'           + 
                          2*tab + 'endif \n')

field_line             = T(2*tab+ \
'!----------------------------------------------------------------------\n')
field_calc_a_start = T( 2*tab + '! Calculation of $FIELD \n')
field_calc_b_start = T( 3*tab + '$FIELD(:$IND,it) = $FIELD(:$IND,it)  & \n')

field_calc_den_a     = T('* sum($DENSITY(:$DENIND,:),$SUMIND)  ')
field_calc_den_b     = T('* $DENSITY(:$DENIND,it)')
field_calc_den_c     = T('* $DENSITY(:$DENIND,3-it)')

# Note the convention for coupling constants for the pairing term is 
# different. We do not deal with isoscalar and isovector coupling constants
# but rather with the coupling constants for neutrons and protons.
field_calc_den_pair  = T('* $DENSITY(:$DENIND,it)')

isoloop = T(2*tab + 'maxit = 2 \n' + \
        2*tab + 'if((.not.calcall).and.(any(${FIELD}.ne.0.0))) maxit=0\n'+\
        2*tab + 'do it=1,maxit \n')

isoloop_end = 2*tab + 'enddo\n'

field_calc_b = T( 3*tab + '& $SIGN $DD $CPLCTE(1,2)  $EXPR1 & \n') 
field_calc_c = T( 3*tab + '& $SIGN $DD $CPLCTE(2,2)  $EXPR2 & \n') 
field_calc_d = T( 3*tab + '& $SIGN $DD $CPLCTE(2,2)  $EXPR3 & \n') 

field_pair_a = T( 3*tab + '& $SIGN $DD $CPLCTE(it,1)   $EXPR2  & \n') 
field_pair_b = T( 3*tab + '& $SIGN $DD $CPLCTE(3-it,1) $EXPR3  & \n') 

doloop    = 2*tab + 'do %s = 1, 3 \n'
enddoloop = 2*tab + 'enddo \n'

field_write_a = T(tab + ('write(chan, iostat=io) "$FIELDFILLED" \n'))
field_write_b = T(tab + ('write(chan, iostat=io)  $FIELD  \n'))

field_read_a  = T(2*tab + ('case("$FIELD") \n'))
field_read_b  = T(3*tab +  'read(chan, iostat=io)  ${FIELD}_hist \n')
field_read_c  = T(3*tab +  'if(symtransfo_needed) then    \n')
field_read_d  = T(4*tab + '$FIELD = ${FIELD}_hist \n'                  \
                            +   4*tab + '$UNDOREAD deallocate($FIELD, ${FIELD}_hist) \n')
field_read_e  = T(3*tab +  'else \n')
field_read_f  = T(3*tab +  'endif \n')

field_transfo = \
T( \
       + 4*tab + 'do it=1,2 \n'                                                                 \
       + 5*tab + '${FIELD}(:$IND,it) = & \n'                                                    \
       + 5*tab + '&  changeboxsize_function(${FIELD}_hist(:$IND,it), filenx, fileny, filenz) \n'\
       + 4*tab + 'enddo \n'                                                              
       )

field_transfo_end = \
T( \
       + 4*tab + 'deallocate(${FIELD}_hist)              \n' \
       + 4*tab + 'allocate(${FIELD}_hist(mv$ALLOCIND,2)) \n' \
       + 4*tab + '${FIELD}_hist = 0.0d0 \n' )


# T for cleaning fields
clean   = T(   tab+'if(allocated($FIELD)) deallocate($FIELD)')
clean_b = T(   tab+'if(allocated(${FIELD}_hist)) deallocate(${FIELD}_hist)')

