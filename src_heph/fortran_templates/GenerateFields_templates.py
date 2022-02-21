#
#
#-------------------------------------------------------------------------------

from string import Template as T    
tab = '    '

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Declaration and allocation statements
field_decl = T(tab + 'real(KIND=dp), allocatable, target :: $FIELD(:$DECLIND,:) \n')
fhist_decl = T(tab + 'real(KIND=dp), allocatable :: ${FIELD}_hist(:$DECLIND,:) \n')

field_allo = T( tab + 'if(.not.allocated($FIELD)) then\n'+\
               2*tab + 'allocate($FIELD(mv$ALLOCIND,$ISOSIZE)) \n'   +\
               2*tab + 'allocate(${FIELD}_hist(mv$ALLOCIND,$ISOSIZE)) \n'+\
               2*tab + '$FIELD = 0.0 ; ${FIELD}_hist = 0.0 \n' + \
                 tab + 'endif \n')

field_allo_b = T(3*tab + 'if(.not.allocated($FIELD)) then\n'+\
                     4*tab + 'allocate($FIELD(mv$ALLOCIND,$ISOSIZE)) \n'   +\
                     4*tab + 'allocate(${FIELD}_hist(filemv$ALLOCIND,$ISOSIZE)) \n'+\
                     4*tab + '$FIELD = 0.0 ; ${FIELD}_hist = 0.0 \n' + \
                     3*tab + 'endif \n')

field_hist   = T(   tab + 'if(calcall) then \n' + 
                  2*tab + '${FIELD}_hist = $FIELD \n' + 
                  2*tab + '$FIELD = 0.0 \n'           + 
                    tab + 'endif \n')

field_condition_start = T( tab + 'if(calcall .or. (all($FIELD .eq. 0.0))) then \n')
field_condition_end   = T( tab + 'endif \n')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculation statements
field_line           = T(tab+ '!' + 75 * '-' +  '\n')
field_calc_a_start   = T(tab + '! Calculation of $FIELD \n')
field_calc_iso_start = T(2*tab + '!' + 36 * '- '  + '\n' + 
                         2*tab + '! Isospin = $ISO \n')
field_calc_b_start   = T(2*tab + '$FIELD(:$IND,$ISOIND) = $FIELD(:$IND,$ISOIND)  & \n')
field_calc_den       = T(' * $DENSITY(:$DENIND,$ISOALT)')
field_calc_DD        = T(' * pow($DENSITY(:$DENIND,$ISOALT), $DD)')

field_calc_full      = T(2*tab + '& $SIGN $CPLCTE $EXPR1 $EXTRA & \n') 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Recombination statements
field_recombination = T(2*tab + '!' + 36 * '- '  + '\n' + 
                        2*tab + '! Recombination \n' +
                        2*tab + '$FIELD(:$DECLIND,1) = $FIELD(:$DECLIND,3) + $FIELD(:$DECLIND,4)\n' +
                        2*tab + '$FIELD(:$DECLIND,2) = $FIELD(:$DECLIND,3) - $FIELD(:$DECLIND,4) \n' )

#-------------------------------------------------------------------------------
# Note the convention for coupling constants for the pairing term is 
# different. We do not deal with isoscalar and isovector coupling constants
# but rather with the coupling constants for neutrons and protons.
field_calc_den_pair  = T('* $DENSITY(:$DENIND,it)')

#-------------------------------------------------------------------------------
# Statements for writing the field to file
field_write_a = T(tab + ('write(chan, iostat=io) "$FIELDFILLED" \n'))
field_write_b = T(tab + ('write(chan, iostat=io)  $FIELD  \n'))
# ... and to read it from file
field_read_a  = T(2*tab + ('case("$FIELD") \n'))
field_read_b  = T(3*tab +  'read(chan, iostat=io)  ${FIELD}_hist \n')
field_read_c  = T(3*tab +  'if(symtransfo_needed) then    \n')
field_read_d  = T(4*tab + '$FIELD = ${FIELD}_hist \n'                  \
                            +   4*tab + '$UNDOREAD deallocate($FIELD, ${FIELD}_hist) \n')
field_read_e  = T(3*tab +  'else \n')
field_read_f  = T(3*tab +  'endif \n')
# ..... and to transform fields with different symmetries
field_transfo = \
T( \
       + 4*tab + 'do it=1,2 \n'                                                                 \
       + 5*tab + '${FIELD}(:$IND,it) = & \n'                                                    \
       + 5*tab + '&  changeboxsize_function(${FIELD}_hist(:$IND,it), filenx, fileny, filenz) \n'\
       + 4*tab + 'enddo \n'
       )

field_transfo_recomb = \
T( \
       + 4*tab + '${FIELD}(:$DECLIND,3) = ${FIELD}(:$DECLIND,1) + ${FIELD}(:$DECLIND,2) \n' \
       + 4*tab + '${FIELD}(:$DECLIND,4) = ${FIELD}(:$DECLIND,1) - ${FIELD}(:$DECLIND,2) \n')
#field_transfo_end = \
#T( \
#       + 4*tab + 'deallocate(${FIELD}_hist)              \n' \
#       + 4*tab + 'allocate(${FIELD}_hist(mv$ALLOCIND,2)) \n' \
#       + 4*tab + '${FIELD}_hist = 0.0d0 \n' )

#-------------------------------------------------------------------------------
# Template for cleaning fields after a run
clean   = T(   tab+'if(allocated($FIELD)) deallocate($FIELD)')
clean_b = T(   tab+'if(allocated(${FIELD}_hist)) deallocate(${FIELD}_hist)')

