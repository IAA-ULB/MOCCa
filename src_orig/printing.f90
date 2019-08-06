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
    20 format (80 ('-'))
    30 format (80 ('_'),/,3x , 'Neutron wavefunctions')
    40 format (80 ('_'),/,3x , 'Proton  wavefunctions')
    60 format (2x,'i',4x,'P',3x,'occ',7x,'E',8x,'d2h',4x,'Delta',3x,'JxT',4x,&
    &             'JyT', 4x ,'Jz', 5x, 'J')    

    11 format (i3, 1x, f4.1, 2x, f6.4, 1x, f9.3, 1x, es8.1,1x,f6.2, 4(2x, f5.2))

    integer       :: wave,k
    integer       :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp) :: p, Jx, Jy, Jz, JJ
  
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! HFBasis
  
    ! Order the spwfs according to growing energy
    ProtonOrder = OrderSpwfsISO(+1)
    NeutronOrder= OrderSpwfsISO(-1)
    
    print 10
    print 30
    print 60
    print 10
    do k=1,nwn 
        wave = NeutronOrder(k)
        
        if(wave .lt. HFBlocks(1)) p = +1
        if(wave .gt. HFBlocks(1)) p = -1

        !-----------------------------------------------------------------------
        ! Depending on the symmetries, select different quantities to print 
        ! for Jx, Jy, Jz. Currently configured for EV8 mode!
        Jx = angmom_xt_real(HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))
        Jy = angmom_yt_imag(HFPsi(:,:,wave),HFPsi(:,:,wave),HFDPsi(:,:,:,wave))
        Jz = angmom_z_real (HFPsi(:,:,wave),HFPsi(:,:,wave),HFdPsi(:,:,:,wave))
        
        JJ = & 
        &   angmom_x_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) &
        & + angmom_y_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) &
        & + angmom_z_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) 
        JJ = (-1. + sqrt(1. + 4*JJ))/2.

        if(pairingtype.eq.1) then
          print 11, wave, p, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), BCSgaps(wave),Jx, Jy, Jz, JJ
        elseif(pairingtype.eq.2) then
          print 11, wave, p,2* rho_pairing(wave,wave), spenergies(wave),       &
          &               dispersions(wave), maxval(abs(HFBgaps(wave,:))),     &
          &               Jx, Jy, Jz, JJ
        else
          print 11, wave, p, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), 0.0, Jx, Jy, Jz, JJ                
        endif
    enddo
    
    print 40  
    print 60
    print 10
    do k=1,nwp
        wave = ProtonOrder(k)
        
        if(wave .lt. sum(HFBlocks(1:5))) p = +1
        if(wave .gt. sum(HFBlocks(1:5))) p = -1

        !-----------------------------------------------------------------------
        ! Depending on the symmetries, select different quantities to print 
        ! for Jx, Jy, Jz. Currently configured for EV8 mode!
        Jx = angmom_xt_real(HFPsi(:,:,wave),HFPsi(:,:,wave),HFDPsi(:,:,:,wave))
        Jy = angmom_yt_imag(HFPsi(:,:,wave),HFPsi(:,:,wave),HFDPsi(:,:,:,wave))
        Jz = angmom_z_real (HFPsi(:,:,wave),HFPsi(:,:,wave),HFDPsi(:,:,:,wave))

        JJ = & 
        &   angmom_x_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) &
        & + angmom_y_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) &
        & + angmom_z_quad(HFPsi(:,:,wave),HFdPsi(:,:,:,wave), &
        &                 HFPsi(:,:,wave),HFdPsi(:,:,:,wave)) 
        JJ = -0.5*(1. - sqrt(1. + 4*JJ))       

        if(pairingtype.eq.1) then
          print 11, wave, p, rho_can(wave), spenergies(wave),                  &
          &               dispersions(wave), BCSgaps(wave),Jx, Jy, Jz, JJ
        elseif(pairingtype.eq.2) then
          print 11, wave, p, 2*rho_pairing(wave,wave), spenergies(wave),       &
          &               dispersions(wave), maxval(abs(HFBgaps(wave,:))),     &
          &               Jx, Jy, Jz, JJ
        else
          print 11, wave, p, rho_can(wave), spenergies(wave),                  &
          &               dispersions(wave), 0.0, Jx, Jy, Jz, JJ
        endif
    enddo
    print 20
    if(PairingType.ne.2) return
  end subroutine PrintSpwfs

  subroutine printqps
    !---------------------------------------------------------------------------
    ! Print all relevant info on quasiparticles.
    ! Very bare-bones for the moment.
    !---------------------------------------------------------------------------
    integer :: i, N, B, si, sb
    
    1  format (33 ('-'), 'Quasiparticles',33('-'))
    2  format (80 ('_'))
    3  format ( i3, 1f7.2, 1es12.2 )
    4  format ('Block ', i1, /,  a8, ' parity ', a8)
    5  format ( '  N    Eqp     f_n')

    11  format(80 ('-'))

    if(PairingType.eq.0) return

    print 1
    
    si = 0
    sb = 0
    do B=1,8
        N = HFblocks(B)

        if(N.eq.0) cycle
        
        print 2
        select case (B)
        case(1)
            print 4,  B , 'positive', 'neutrons' 
        case(3)
            print 4,  B , 'negative' , 'neutrons'
        case(5)        
            print 4,  B , 'positive' , 'protons'
        case(7)        
            print 4,  B , 'negative' , 'protons'
        end select
        print 5
        print 2
        do i=1,N
            select case(pairingtype)
            case(2)
              print 3, i, QPenergies(si+i), configmatrix(sb+i)
            case(1)
              print 3, i, BCSqps(si+i), BCSf(si+i)
            end select
        enddo
        si = si + N
        sb = sb + 2*N
    enddo
    print 11


  end subroutine printqps
  
end module
