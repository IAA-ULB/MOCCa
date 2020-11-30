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
  
    ! Order the spwfs according to growing energy
    ProtonOrder = OrderSpwfsISO(+1)
    NeutronOrder= OrderSpwfsISO(-1)

    print 10
    print 30
    print 60
    print 10
    do k=1,nwn 
        wave = NeutronOrder(k)
        
        if(wave .le. HFBlocks(1)) p = +1
        if(wave .gt. HFBlocks(1)) p = -1

        Jx = spwf_JTR(1,wave)
        Jy = spwf_JTI(2,wave)
        Jz = spwf_J(3,wave)
        JJ = spwf_JJ(wave)

        if(pairingtype.eq.1) then
          print 11, wave, p, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), BCSgaps(wave),Jx, Jy, Jz, JJ
        elseif(pairingtype.eq.2) then
          print 11, wave, p,   rho_pairing(wave,wave), spenergies(wave),       &
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
        
        if(wave .le. sum(HFBlocks(1:5))) p = +1
        if(wave .gt. sum(HFBlocks(1:5))) p = -1

        Jx = spwf_JTR(1,wave)
        Jy = spwf_JTI(2,wave)
        Jz = spwf_J(3,wave)
        JJ = spwf_JJ(wave)

        if(pairingtype.eq.1) then
          print 11, wave, p, rho_can(wave), spenergies(wave),                  &
          &               dispersions(wave), BCSgaps(wave),Jx, Jy, Jz, JJ
        elseif(pairingtype.eq.2) then
          print 11, wave, p,     rho_pairing(wave,wave), spenergies(wave),     &
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
    integer :: i, N, B, si, sb, ind, N2
    integer, allocatable :: indices(:)
    
    1  format (33 ('-'), 'Quasiparticles',33('-'))
    2  format ( i3, 1f7.2, 1es12.2 )

    11  format(80 ('-'))

    if(PairingType.eq.0) return

    print 1
    
    si = 0
    sb = 0
    do B=1,8,2
        N = HFblocks(B) ;      if(N.eq.0) cycle
        N2 = HFBlocks(B+1)

        call print_qp_header(B)
        select case(pairingtype)
        case(2)
          do i=1,N
            print 2, i, QPenergies(si+i), configmatrix(sb+i)
          enddo
          if(N2.ne.0) then 
            print *
            call print_qp_header(B+1)
            do i=1,N2
              print 2, i, QPenergies(si+N+i), configmatrix(sb+N+i)
            enddo
          endif
        case(1)
          ! The BCS qp energies are not ordered by energy
          indices = order(BCSqps(si+1:si+N))
          do i=1, N
            ind  = indices(i)
            print 2, i, BCSqps(si+ind), BCSf(si+ind)
          enddo
        end select
        si = si +   N +  N2
        sb = sb + 2*N +2*N2
        print *
    enddo
    print 11
  end subroutine printqps

  subroutine print_qp_header(B)
    
    integer, intent(in) :: B  

    1  format ('Block ', i1, ':  P=',a1,'1', ',Rz=', a1,'i ',  a8)
    2  format ( '  N    Eqp     f_n')
    3  format (80 ('_'))
  
    select case (B)
    case(1)
        print 1,  B , '+', '+', 'neutrons' 
    case(2)
        print 1,  B , '+', '-', 'neutrons' 
    case(3)
        print 1,  B , '-' ,'+', 'neutrons'
    case(4)
        print 1,  B , '-' ,'-', 'neutrons'
    case(5)        
        print 1,  B , '+' ,'+',  'protons'
    case(6)        
        print 1,  B , '+' ,'-',  'protons'
    case(7)        
        print 1,  B , '-' ,'+',  'protons'
    case(8)        
        print 1,  B , '-' ,'-',  'protons'
    end select
    print 2
    print 3
  end subroutine print_qp_header

end module
