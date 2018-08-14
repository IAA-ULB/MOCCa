module Printing
 !==============================================================================
 !  #######   ##   #    # #####   ##   #      #    #  ####
 !     #     #  #  ##   #   #    #  #  #      #    # #
 !     #    #    # # #  #   #   #    # #      #    #  ####
 !     #    ###### #  # #   #   ###### #      #    #      #
 !     #    #    # #   ##   #   #    # #      #    # #    #
 !     #    #    # #    #   #   #    # ######  ####   ####
 !
 !  Copyright W. Ryssens & M. Bender
 !
 !==============================================================================
 use pairing
 use wavefunctions

 implicit none
 
contains

  subroutine PrintSpwfs
    !---------------------------------------------------------------------------
    ! Print the info on the single-particle wave-functions in the HFBasis.
    !---------------------------------------------------------------------------
    
    10 format (21 ('-'), ' Sp wavefunctions ', 41('-'))
    20 format (90 ('-'))
    30 format (90 ('_'),/,3x , 'Neutron wavefunctions')
    40 format (90 ('_'),/,3x , 'Proton  wavefunctions')
    50 format (90 ('_'),/,3x , 'HF Basis')
    60 format (90 ('_'),/,3x , 'CAN Basis')
    
    11 format (i3, 3x, f5.2, 3x, f7.4, 3x, f10.3, 3x, e10.3, 3x, f7.4)

    integer       :: wave,k
    integer       :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp) :: p
  
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! HFBasis
  
    ! Order the spwfs according to growing energy
    ProtonOrder = OrderSpwfsISO(+1)
    NeutronOrder= OrderSpwfsISO(-1)
    
    print 10
    print 50
    print 30
    print 10
    do k=1,nwn 
        wave = NeutronOrder(k)
        
        if(wave .lt. HFBlocks(1)) p = +1
        if(wave .gt. HFBlocks(1)) p = -1
        
        if(pairingtype.eq.1) then
          print 11, wave, p, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), BCSgaps(wave)
        elseif(pairingtype.eq.2) then
          print 11, wave, p,2* rho_pairing(wave,wave), spenergies(wave), &
          &               dispersions(wave), maxval(abs(HFBgaps(wave,:)))
        else
          print 11, wave, p, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), 0.0
        endif
    enddo
    
    print 40  
    print 10
    do k=1,nwp
        wave = ProtonOrder(k)
        
      if(wave .lt. sum(HFBlocks(1:5))) p = +1
      if(wave .gt. sum(HFBlocks(1:5))) p = -1
        
        if(pairingtype.eq.1) then
          print 11, wave, p, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), BCSgaps(wave)
        elseif(pairingtype.eq.2) then
          print 11, wave, p, 2*rho_pairing(wave,wave), spenergies(wave), &
          &               dispersions(wave), maxval(abs(HFBgaps(wave,:)))
        else
          print 11, wave, p, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), 0.0
        endif
    enddo
    print 20
    if(PairingType.ne.2) return
  end subroutine PrintSpwfs
  
end module
