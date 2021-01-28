module Printing
 !==============================================================================
 !_________ _______  _       _________ _______  _                 _______ 
 !\__   __/(  ___  )( (    /|\__   __/(  ___  )( \      |\     /|(  ____ \
 !   ) (   | (   ) ||  \  ( |   ) (   | (   ) || (      | )   ( || (    \/
 !   | |   | (___) ||   \ | |   | |   | (___) || |      | |   | || (_____ 
 !   | |   |  ___  || (\ \) |   | |   |  ___  || |      | |   | |(_____  )
 !   | |   | (   ) || | \   |   | |   | (   ) || |      | |   | |      ) |
 !   | |   | )   ( || )  \  |   | |   | )   ( || (____/\| (___) |/\____) |
 !   )_(   |/     \||/    )_)   )_(   |/     \|(_______/(_______)\_______)
 !                                                                       
 !  Copyright W. Ryssens & M. Bender
 !
 !==============================================================================
 use pairing
 use wavefunctions
 use convergence

 implicit none
 
contains

  subroutine PrintSpwfs
    !---------------------------------------------------------------------------
    ! Print the info on the single-particle wave-functions in the HFBasis.
    !---------------------------------------------------------------------------
    
    10 format (21 ('-'), ' Sp wavefunctions ', 41('-'))
    12 format (21 ('-'), ' Canonical basis  ', 41('-'))
    20 format (80 ('-'))
    30 format (80 ('_'),/,3x , 'Neutron wavefunctions')
    40 format (80 ('_'),/,3x , 'Proton  wavefunctions')
    60 format (2x,'i',4x,'P',3x, 'Rz', 3x,'occ',7x,'E',8x,'d2h',4x,'Delta',3x, &
    &             'JxT',4x, 'JyT', 4x ,'Jz', 5x, 'J')    

    11 format (i3, 1x, f4.1, 1x, f4.1, 2x, f6.4, 1x, f9.3, 1x, es8.1,1x,f6.2, 4(2x, f5.2))

    integer       :: wave,k
    integer       :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp) :: p, Jx, Jy, Jz, JJ, s
  
    ! Order the spwfs according to growing energy
    ProtonOrder = OrderSpwfsISO(+1)
    NeutronOrder= OrderSpwfsISO(-1)

    print 10
    print 30
    print 60
    print 10
    do k=1,nwn 
        wave = NeutronOrder(k)
        
        if(wave .le. sum(HFBlocks(1:2))) p = +1
        if(wave .gt. sum(HFBlocks(1:2))) p = -1


        if(wave .le. sum(HFBlocks(1:2))) then
            if(wave .le. HFBlocks(1)) then
               s = +1
            else
               s = -1
            endif
        else
            if(wave .le. sum(HFBlocks(1:3))) then
               s = +1
            else
               s = -1
            endif
        endif

        Jx = spwf_JTR(1,wave)
        Jy = spwf_JTI(2,wave)
        Jz = spwf_J(3,wave)
        JJ = spwf_JJ(wave)

        if(pairingtype.eq.1) then
          print 11, wave, p, s, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), BCSgaps(wave),Jx, Jy, Jz, JJ
        elseif(pairingtype.eq.2) then
          print 11, wave, p, s, rho_pairing(wave,wave), spenergies(wave),      &
          &               dispersions(wave), maxval(abs(HFBgaps(wave,:))),     &
          &               Jx, Jy, Jz, JJ
        else
          print 11, wave, p, s, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), 0.0, Jx, Jy, Jz, JJ                
        endif
    enddo
    
    print 40  
    print 60
    print 10
    do k=1,nwp
        wave = ProtonOrder(k)
        
        if(wave .le. sum(HFBlocks(1:6))) p = +1
        if(wave .gt. sum(HFBlocks(1:6))) p = -1

        if(wave .le. sum(HFBlocks(1:6))) then
            if(wave .le. sum(HFBlocks(1:5))) then
               s = +1
            else
               s = -1
            endif
        else
            if(wave .le. sum(HFBlocks(1:7))) then
               s = +1
            else
               s = -1
            endif
        endif

        Jx = spwf_JTR(1,wave)
        Jy = spwf_JTI(2,wave)
        Jz = spwf_J(3,wave)
        JJ = spwf_JJ(wave)

        if(pairingtype.eq.1) then
          print 11, wave, p, s, rho_can(wave), spenergies(wave),               &
          &               dispersions(wave), BCSgaps(wave),Jx, Jy, Jz, JJ
        elseif(pairingtype.eq.2) then
          print 11, wave, p, s,  rho_pairing(wave,wave), spenergies(wave),     &
          &               dispersions(wave), maxval(abs(HFBgaps(wave,:))),     &
          &               Jx, Jy, Jz, JJ
        else
          print 11, wave, p, s, rho_can(wave), spenergies(wave),               &
          &               dispersions(wave), 0.0, Jx, Jy, Jz, JJ
        endif
    enddo
    print 20
  
    ! Return if we are not doing a HFB calculation
    if(PairingType.ne.2) return

    ! Otherwise, print the properties of the canonical basis
    print 12
    print 30
    print 60
  
    ! Order the canonical basis, not the HF one
    ProtonOrder = OrderSpwfsISO(+1, .true.) 
    NeutronOrder= OrderSpwfsISO(-1, .true.)

    do k=1,nwn 
      wave = NeutronOrder(k) 
      if(wave .le. sum(HFBlocks(1:2))) p = +1
      if(wave .gt. sum(HFBlocks(1:2))) p = -1

      if(wave .le. sum(HFBlocks(1:2))) then
          if(wave .le. HFBlocks(1)) then
             s = +1
          else
             s = -1
          endif
      else
          if(wave .le. sum(HFBlocks(1:3))) then
             s = +1
          else
             s = -1
          endif
      endif

      Jx = can_JTR(1,wave)
      Jy = can_JTI(2,wave)
      Jz = can_J(3,wave)
      JJ = can_JJ(wave)

      print 11, wave, p,  s,   rho_can(wave), canenergies(wave),             &
      &               0.0, 0.0, Jx, Jy, Jz, JJ
    enddo
    print 40  
    print 60
    print 10
    do k=1,nwp 

      wave = ProtonOrder(k) 

      if(wave .le. sum(HFBlocks(1:6))) p = +1
      if(wave .gt. sum(HFBlocks(1:6))) p = -1

      if(wave .le. sum(HFBlocks(1:6))) then
          if(wave .le. sum(HFBlocks(1:5))) then
             s = +1
          else
             s = -1
          endif
      else
          if(wave .le. sum(HFBlocks(1:7))) then
             s = +1
          else
             s = -1
          endif
      endif

      Jx = can_JTR(1,wave)
      Jy = can_JTI(2,wave)
      Jz = can_J(3,wave)
      JJ = can_JJ(wave)

      print 11, wave, p, s,    rho_can(wave), canenergies(wave),             &
      &               0.0, 0.0, Jx, Jy, Jz, JJ
    enddo
    print 20

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
          ! HFB QP energies
          do i=1,N
            print 2, i, QPenergies(sb+i), configmatrix(sb+N-i+1)
          enddo
          do i=N+N2+1,2*N + N2
            print 2, i, QPenergies(sb+i), configmatrix(sb+i)
          enddo

          if(N2.ne.0) then 
            print *
            call print_qp_header(B+1)
            do i=1,N2
              print 2, i, QPenergies(sb+N+i), configmatrix(sb+N+N2-i+1)
            enddo
            do i=2*N+N2+1,2*N+2*N2
              print 2, i, QPenergies(sb+i), configmatrix(sb+i)
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

  subroutine convergence_report()
    !---------------------------------------------------------------------------
    ! Print a report on the observed convergence.
    !---------------------------------------------------------------------------

    use functional

    1 format (80('-'))
    2 format ('Convergence report')
    3 format (19x,'     Change^(i)     Change^(i-1)   rate (approx)')
    4 format ('  Con.  Energy      :', 3es15.2)
    5 format ('  Con.  Routhian    :', 3es15.2)
    6 format ('  Con.  E_fu - E_sp :', 3es15.2)
    7 format ('  Con.  |delta rho| :', 3es15.2)
    
    print 1
    print 2
    print 3
    print 4, totalE   - Ehistory(1), Ehistory(1) - Ehistory(2), con_rates(1)
    print 5, Routhian - Rhistory(1), Rhistory(1) - Rhistory(2), con_rates(2)
    print 6, SpwfEnergy     - totalE      - SpwfHistory(1) + Ehistory(1), &
    &        SpwfHistory(1) - Ehistory(1) - SpwfHistory(2) + Ehistory(2), &
    &        con_rates(3)
    print 7,  sqrt(sum((D_I_I - D_I_I_hist(:,:,1))**2)*dv) , 0.0 , con_rates(4)
    print 1

  end subroutine convergence_report

end module
