from string import Template as T    
tab = '    '

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Declaration and allocation statements
field_decl = T(tab + 'real(KIND=dp), allocatable :: $FIELD(:$DECLIND,:) \n')

field_allo = T( tab + 'if(.not.allocated(F%$FIELD)) then\n'+\
               2*tab + 'allocate(F%$FIELD(mv$ALLOCIND,$ISOSIZE)) \n'   +\
               2*tab + 'F%$FIELD = 0.0 \n' + \
                 tab + 'endif \n')

field_condition_start = T(  \
                       tab + 'if(.not. allocated(F%$FIELD)) then\n'+ \
                     2*tab + 'allocate(F%$FIELD(mv$ALLOCIND,$ISOSIZE)) \n'+ \
                     2*tab + 'F%$FIELD = 0.0 \n')
field_condition_end   = T( tab + 'endif \n')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Calculation statements
field_line           = T(tab+ '!' + 75 * '-' +  '\n')
field_calc_a_start   = T(tab + '! Calculation of $FIELD \n')
field_calc_iso_start = T(2*tab + '!' + 36 * '- '  + '\n' + 
                         2*tab + '! Isospin = $ISO \n')
field_calc_b_start   = T(2*tab + 'F%$FIELD(:$IND,$ISOIND) = F%$FIELD(:$IND,$ISOIND)  & \n')
field_calc_den       = T(' * R%$DENSITY(:$DENIND,$ISOALT)')
field_calc_den_pert  = T(' * R_pert%$DENSITY(:$DENIND,$ISOALT)')
field_calc_DD        = T(' * pow(R%$DENSITY(:$DENIND,$ISOALT), $DD)')
#field_calc_DD_pert   = T(' * $DD * pow(R%$DENSITY(:$DENIND,$ISOALT), $DD_pert) * R_pert%$DENSITY(:$DENIND,$ISOALT)')
field_calc_DD_pert   = T(' * lr_pow(R%$DENSITY(:$DENIND,$ISOALT), R_pert%$DENSITY(:$DENIND,$ISOALT), $DD_pert)')

field_calc_full      = T(2*tab + '& $SIGN $CPLCTE $EXPR1 $EXTRA & \n') 
field_calc_pert      = T(2*tab + '& $SIGN $CPLCTE $EXPR_PERT $EXTRA & \n') 
field_calc_INM = T(  tab + 'pot(:,$ISOIND) = pot(:,$ISOIND) $SIGN $CPLCTE $EXPR1 \n')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Dealing with potential vectors
Add       = T( tab + ' F%$FIELD = F1%$FIELD + F2%$FIELD')
Multiply  = T( tab + ' F%$FIELD = a * F1%$FIELD')

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Recombination statements
field_recombination = T(2*tab + '!' + 36 * '- '  + '\n' + 
                        2*tab + '! Recombination \n' +
                        2*tab + 'F%$FIELD(:$DECLIND,1) = F%$FIELD(:$DECLIND,3) +F%$FIELD(:$DECLIND,4)\n' +
                        2*tab + 'F%$FIELD(:$DECLIND,2) = F%$FIELD(:$DECLIND,3) -F%$FIELD(:$DECLIND,4) \n' )

#-------------------------------------------------------------------------------
# Note the convention for coupling constants for the pairing term is 
# different. We do not deal with isoscalar and isovector coupling constants
# but rather with the coupling constants for neutrons and protons.
field_calc_den_pair  = T('* $DENSITY(:$DENIND,it)')

#-------------------------------------------------------------------------------
# Statements for writing the field to file
field_write_a = T(tab + ('write(chan, iostat=io) "$FIELDFILLED" \n'))
field_write_b = T(tab + ('write(chan, iostat=io)  F%$FIELD  \n'))
field_write_hdf5 = T(tab + ('call hdf5_writepot(file_id,"$FIELDFILLED",F%$FIELD,size(F%$FIELD)) \n')) 
# ... and to read it from file
field_read_a  = T(2*tab + ('case("$FIELD") \n'))
field_read_b  = T(3*tab +  'allocate(F%${FIELD}(mv$ALLOCIND,$ISOSIZE),F_temp%${FIELD}(mv$ALLOCIND,$ISOSIZE))\n' +
                  3*tab +  'if(MPI_RANK .eq. 0) read(chan, iostat=io) F_temp%${FIELD}\n' + 
                  '#if (USE_MPI > 0) \n' +
                  3*tab +  'call MPI_BCAST(F_temp%${FIELD},size(F_temp%${FIELD}),MPI_REAL8,0,MPI_COMM_WORLD,mpi_err) \n' + 
                  '#endif \n')

field_read_c  = T(3*tab +  '!if(symtransfo_needed) then    \n')
field_read_d  = T(4*tab + '$POTREAD F%$FIELD = F_temp%${FIELD} \n'            \
                + 4*tab + '$UNDOREAD deallocate(F%$FIELD, F_temp%${FIELD}) \n')
field_read_e  = T(3*tab +  '!else \n')
field_read_f  = T(3*tab +  '!endif \n')
# hdf5 option
field_read_hdf5_a = T(tab +  'allocate(F%${FIELD}(mv$ALLOCIND,$ISOSIZE),F_temp%${FIELD}(mv$ALLOCIND,$ISOSIZE))\n')
field_read_hdf5_b = T(tab +  'if(MPI_RANK .eq. 0) then \n' +
                    2*tab +  'allocate(Ftmp(mv$ALLOCINDHDF5*$ISOSIZE)) \n' +
                    2*tab +  'call hdf5_readpot(group_id, "$FIELD", Ftmp, size(Ftmp))\n' +
                    2*tab +  'F_temp%${FIELD}=reshape(Ftmp,(/mv$ALLOCIND,$ISOSIZE/))\n' +
                    2*tab +  'deallocate(Ftmp)                                 \n' +
                    tab +  'endif                              \n'+
                    '#if (USE_MPI > 0) \n' +
                    tab +  'call MPI_BCAST(F_temp%${FIELD},size(F_temp%${FIELD}),MPI_REAL8,0,MPI_COMM_WORLD,mpi_err) \n' + 
                  '#endif \n')
field_read_hdf5_c = T(tab +  'if(symtransfo_needed) then    \n')
field_read_hdf5_d = T(2*tab + '$POTREAD F%$FIELD = F_temp%${FIELD} \n'            \
                + 2*tab + '$UNDOREAD deallocate(F%$FIELD, F_temp%${FIELD}) \n')
field_read_hdf5_e = T(tab +  'else \n')
field_read_hdf5_f = T(tab +  'endif \n')
# ..... and to transform fields with different symmetries
field_transfo = \
T( \
       + 4*tab + '!do it=1,2 \n'                                                                   \
       + 5*tab + '!F%${FIELD}(:$IND,it) = & \n'                                                    \
       + 5*tab + '!&  changeboxsize_function(F_temp%${FIELD}(:$IND,it), filenx, fileny, filenz) \n'\
       + 4*tab + '!enddo \n'
       )

field_transfo_recomb = \
T( \
       + 4*tab + 'F%${FIELD}(:$DECLIND,3) = F%${FIELD}(:$DECLIND,1) + F%${FIELD}(:$DECLIND,2) \n' \
       + 4*tab + 'F%${FIELD}(:$DECLIND,4) = F%${FIELD}(:$DECLIND,1) - F%${FIELD}(:$DECLIND,2) \n')
field_transfo_hdf5 = \
T( \
       + 2*tab + 'do it=1,2 \n'                                                                   \
       + 3*tab + '!F%${FIELD}(:$IND,it) = & \n'                                                    \
       + 3*tab + '!&  changeboxsize_function(F_temp%${FIELD}(:$IND,it), filenx, fileny, filenz) \n'\
       + 2*tab + 'enddo \n'
       )

field_transfo_recomb_hdf5 = \
T( \
       + 2*tab + 'F%${FIELD}(:$DECLIND,3) = F%${FIELD}(:$DECLIND,1) + F%${FIELD}(:$DECLIND,2) \n' \
       + 2*tab + 'F%${FIELD}(:$DECLIND,4) = F%${FIELD}(:$DECLIND,1) - F%${FIELD}(:$DECLIND,2) \n')
#-------------------------------------------------------------------------------
# Templates for preconditioning 
field_precon_start= \
T( 1*tab + '! Preconditioning of the field ${FIELD} \n')
field_precon_update = \
T( 1*tab +  'update= F_out%${FIELD}(:$IND,:) - F_in%${FIELD}(:$IND,:) \n')
field_precon_call = \
T( 2*tab +  'update=  PreconditionPotential(update,-preconfactor/$DIVISOR,1.0_dp, $PX,$PY,$PZ) \n')
field_precon_add  = \
T( 1*tab +  'F%${FIELD}(:$IND,:) = F_in%${FIELD}(:$IND,:) + update \n')

#-------------------------------------------------------------------------------
# Template for inproduct
field_inproduct = T(1*tab + 'x = x + sum(F1%$FIELD * F2%$FIELD)\n')
#-------------------------------------------------------------------------------
# Template for cleaning fields after a run
clean   = T(   tab+'if(allocated($FIELD)) deallocate($FIELD)')
clean_b = T(   tab+'!if(allocated(${FIELD}_hist)) deallocate(${FIELD}_hist)')

