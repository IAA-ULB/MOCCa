from string import Template as T
tab = '   '

action_final            = T(  tab + \
'hpsi(:,$IND) =  hpsi(:,$IND) $SIGN $LMULT $TEMP(:$LIND,$RCOMP)\n')
action_final_pairing    = T(  tab + \
'deltapsi(:,$IND) =  deltapsi(:,$IND) $SIGN $LMULT $TEMP(:$LIND,$RCOMP)\n')


temp_ini        = tab + 'temp = 0.0 \n'
action    = T(2*tab + \
'temp(i,$IND) =  temp(i,$IND) $SIGN $RMULT $FIELD(i$FIELDIND,it) * $WF(i$RIND,$RCOMP)\n')

#derive     = T(tab + \
#'call Derive_$DIR($DNUMBER(:,$RCOMP), $SYM, d$DNUMBER(:,$DIRIND,$RCOMP)) \n')
derive     = T(tab + \
'call Derive_$DIR($DNUMBER(:,:), sym, d$DNUMBER(:,$DIRIND,:)) \n')
lap        = T(tab + \
'call Derive_lap(temp(:,$RCOMP), $SYMX, $SYMY, $SYMZ, laptemp(:,$RCOMP)) \n')
sym        = T(tab + \
'sym($RCOMP) = $SYM \n')

position_loop   = tab + 'do i=1,mv\n'
position_end    = tab + 'enddo    \n'

action_comment  = T(tab + '!' + 75*'-' + '\n'\
                    +    tab + '! Action of $FIELD symmetrized: $SYM \n')
