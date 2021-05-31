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
    ! Print the info of the physical Hartree-Fock basis.
    !---------------------------------------------------------------------------
    
    10 format (21 ('-'), ' Sp wavefunctions ', 61('-'))
    12 format (21 ('-'), ' Canonical basis  ', 61('-'))
    20 format (100 ('-'))
    30 format (100 ('_'),/,3x , 'Neutron wavefunctions')
    40 format (100 ('_'),/,3x , 'Proton  wavefunctions')
    60 format (2x,'i',4x,'P',3x, 'Rz', 3x,'occ',7x,'E',8x,'d2h',4x,'Delta',  &
    &             ' | ', 2x, 'JxT',4x, 'JyT', 4x,'Jz', 6x, 'J', 2x,          &
    &             ' | ', 2x, 'SxT',4x, 'SyT', 4x,'Sz')    

    11 format (i3, 1x, f4.1, 1x, f4.1, 2x, f6.4, 1x, f9.3, 1x, es8.1,1x,f6.2,  &
    &          1x,'|', 4(2x, f5.2), 1x, '|', 3(2x, f5.2) )

    integer       :: wave,k, B, si, N, T, wavebar
    integer       :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp) :: p, Jx, Jy, Jz, JJ, s, Delta, Sx, Sy, Sz
    
    real(KIND=dp), allocatable :: HF_gaps(:,:), can_gaps(:,:)
    
    ! We transform the gaps to the Hartree-Fock basis for printing
    if(pairingtype.eq.2) then
      si = 0
      allocate(HF_gaps(nwt,nwt)) ; HF_gaps = 0.0d0
      do B=1,8,2
        N = HFblocks(B) ; if(N.eq.0) cycle
        T = HFBlocks(B+1) + N
        
        HF_gaps(si+1:si+T, si+1:si+T) = &
        & matmul(transpose(HFtransfo(si+1:si+T,si+1:si+T)),&
        &                                          HFBgaps(si+1:si+T,si+1:si+T))
        HF_gaps(si+1:si+T, si+1:si+T) = &
        &    matmul(HF_gaps(si+1:si+T,si+1:si+T),HFtransfo(si+1:si+T,si+1:si+T))

        si = si + T
      enddo
    endif
    !---------------------------------------------------------------------------
    ! Start of the actual printing.   
    ! Order the spwfs according to growing single-particle energy.
    ProtonOrder = OrderSpwfsISO(+1)
    NeutronOrder= OrderSpwfsISO(-1)

    print 10
    print 30
    print 60
    print 20
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

        Jx = HF_JTR(1,wave) ; SX = HF_STR (1,wave)
        Jy = HF_JTI(2,wave) ; SY = HF_STI (2,wave)
        Jz = HF_J(3,wave)   ; SZ = HF_spin(3,wave)
        JJ = HF_JJ(wave)

        if(pairingtype.eq.1) then
          print 11, wave, p, s, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), BCSgaps(wave),      &
          &               Jx, Jy, Jz, JJ, Sx, Sy, Sz
        elseif(pairingtype.eq.2) then
          print 11, wave, p, s, rho_HF(wave), spenergies(wave),                &
          &               dispersions(wave), maxval(abs(HF_gaps(wave,:))),     &
          &               Jx, Jy, Jz, JJ, Sx, Sy, Sz
        else
          print 11, wave, p, s, rho_can(wave), spenergies(wave), &
          &               dispersion(wave), 0.0, Jx, Jy, Jz, JJ, Sx, Sy, Sz              
        endif
    enddo
    
    print 40  
    print 60
    print 20
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

        Jx = HF_JTR(1,wave) ; SX = HF_STR (1,wave)
        Jy = HF_JTI(2,wave) ; SY = HF_STI (2,wave)
        Jz = HF_J(3,wave)   ; SZ = HF_spin(3,wave)
        JJ = HF_JJ(wave)

        if(pairingtype.eq.1) then
          print 11, wave, p, s, rho_can(wave), spenergies(wave),               &
          &               dispersions(wave), BCSgaps(wave),                    &
          &               Jx, Jy, Jz, JJ, Sx, Sy, Sz
        elseif(pairingtype.eq.2) then
          print 11, wave, p, s,  rho_HF(wave), spenergies(wave),               &
          &               dispersions(wave), maxval(abs(HF_gaps(wave,:))),     &
          &               Jx, Jy, Jz, JJ, Sx, Sy, Sz
        else
          print 11, wave, p, s, rho_can(wave), spenergies(wave),               &
          &               dispersions(wave), 0.0, Jx, Jy, Jz, JJ,              &
          &               Sx, Sy, Sz
        endif
    enddo
    print 20
    deallocate( HF_gaps)  
    !---------------------------------------------------------------------------
    ! Return if we are not doing a HFB calculation
    if(PairingType.ne.2) return
    !---------------------------------------------------------------------------
    ! Otherwise, print the properties of the canonical basis
    print 12
    print 30
    print 60
    print 20

    ! Prepare by calculating the gaps in the canonical basis  
    can_gaps = matmul(transpose(cantransfo), HFBgaps)
    can_gaps = matmul(can_gaps, cantransfo)
  
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

      Jx = can_JTR(1,wave) ; SX = can_STR (1,wave)
      Jy = can_JTI(2,wave) ; SY = can_STI (2,wave)
      Jz = can_J(3,wave)   ; SZ = can_spin(3,wave)
      JJ = can_JJ(wave)
    
     if(allocated(conjugp)) then
       wavebar  = conjugp(wave)
     else
       wavebar  = wave     
     endif
      if(wavebar .eq.0) then
          Delta = 0.0
      else
          Delta = can_gaps(wave, wavebar)
      endif    
      print 11, wave, p,  s,   rho_can(wave), canenergies(wave),             &
      &               0.0, Delta , Jx, Jy, Jz, JJ, Sx, Sy, Sz
    enddo
    print 40  
    print 60
    print 20
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


      Jx = can_JTR(1,wave) ; SX = can_STR (1,wave)
      Jy = can_JTI(2,wave) ; SY = can_STI (2,wave)
      Jz = can_J(3,wave)   ; SZ = can_spin(3,wave)
      JJ = can_JJ(wave)
      
      if(allocated(conjugp)) then
        wavebar  = conjugp(wave)
      else
        wavebar  = wave     
      endif
      if(wavebar .eq.0) then
          Delta = 0.0
      else
          Delta = can_gaps(wave, wavebar)
      endif   

      print 11, wave, p, s,    rho_can(wave), canenergies(wave),             &
      &               0.0, Delta, Jx, Jy, Jz, JJ, Sx, Sy, Sz
    enddo
    print 20

  end subroutine PrintSpwfs

  subroutine printqps
    !---------------------------------------------------------------------------
    ! Print all relevant info on quasiparticles.
    ! Rather bare-bones for the moment.
    !---------------------------------------------------------------------------
    integer              :: i, N, B, si, sb, ind, N2, k, U(1),V(1)
    integer, allocatable :: indices(:)
    real(KIND=dp)        :: ov
    character(len=1)     :: Bstr, Pstr
    
    1  format (33 ('-'), 'Quasiparticles',33('-'))
    2  format ( i3, 1f10.2, 2x, 1es12.2,' | ', 2i4, ' | ', 2x, a1, 2x, a1,     &
    &           2x, 1f5.3, ' | ',  3(2x,f5.2))

    11  format(80 ('-'))

    if(PairingType.eq.0) return
    
    call update_qp_angmom(Bogoliubov)

    print 1
    
    si = 0
    sb = 0
    do B=1,8,2
        N = HFblocks(B) ;      if(N.eq.0) cycle
        N2 = HFBlocks(B+1)

        call print_qp_header(B)
        select case(pairingtype)
        case(2)
          ! HFB case
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The unselected quasi-particles
          do i=1,N+N2
            ! What are the single-particles dominating these qps? 
            if(maxval(abs(Bogoliubov(sb     +1:sb+  N+  N2,sb+i))).gt. 0.1) then
              U = maxloc(Bogoliubov(sb     +1:sb+  N+  N2,sb+i)**2)+si
            else
              U = 0
            endif
            if(maxval(abs(Bogoliubov(sb+N+N2+1:sb+2*N+2*N2,sb+i))).gt. 0.1) then
              V = maxloc(Bogoliubov(sb+N+N2+1:sb+2*N+2*N2,sb+i)**2)+si
            else
              V = 0
            endif    
          
            print 2, i, QPenergies(sb+i), 1-configmatrix(sb+2*N+2*N2-i+1),     &
            &           U(1), V(1), '-', '-', 0.0d0,                           &
            &            qp_JTR(1,sb+i), qp_JTI(2,sb+i), QP_J(3,sb+i)
          enddo
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The selected quasi-particles
          do i=N+N2+1,2*N+2*N2
            ! What are the single-particles dominating these qps? 
            if(maxval(abs(Bogoliubov(sb     +1:sb+  N+  N2,sb+i))).gt. 0.1) then
              U = maxloc(Bogoliubov(sb     +1:sb+  N+  N2,sb+i)**2)+si
            else
              U = 0
            endif
            if(maxval(abs(Bogoliubov(sb+N+N2+1:sb+2*N+2*N2,sb+i))).gt. 0.1) then
              V = maxloc(Bogoliubov(sb+N+N2+1:sb+2*N+2*N2,sb+i)**2)+si
            else
              V = 0
            endif          
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Check if this particular qp is a blocked, partner or not special
            if(allocated(blocked_qps)) then
              Bstr = '-'
              Pstr = '-'
              ov   = 0.0d0 
              do k=1,size(blocked_qps)
                if(N2.eq.0) then
                  if(blocked_qps(k) .eq. si+i-N-N2) then
                    Bstr = 'B'
                    Pstr = 'P'
                    ov = 1.0
                  endif
                else
                  if(blocked_qps(k) .eq. si+i-N-N2) then
                     Bstr = 'B'
                     Pstr = '-'
                     ov = partner_overlaps(k)
                  endif
                  if(partner_qps(k) .eq. si+i-N-N2) then
                    Bstr = '-'
                    Pstr = 'P'
                    ov = partner_overlaps(k)
                  endif
                endif
              enddo
            else 
              Bstr = '-'
              Pstr = '-'
              ov = 0.0d0
            endif
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            print 2, i, QPenergies(sb+i), configmatrix(sb+i), U(1), V(1),      &
            &           Bstr, Pstr, ov,                                        &
            &           qp_JTR(1,sb+i), qp_JTI(2,sb+i), QP_J(3,sb+i)
          enddo

        case(1)
          ! The BCS qp energies are not ordered by energy
          indices = order(BCSqps(si+1:si+N))
          do i=1, N
            ind  = indices(i)
            print 2, i, BCSqps(si+ind), BCSf(si+ind), 0,0,'-', '-', 0.0d0
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

    1  format ('Block ', i1, ':  P=',a1,'1',2x,  a8)
    2  format ( '  N      Eqp       f_n      |   U   V  |   B  P  ov_TR |',4x, &
    &           'JxT',4x,'JyT',4x,'Jz')
    3  format (80 ('_'))
  
    select case (B)
    case(1)
        print 1,  B , '+', 'neutrons' 
!    case(2)
!        print 1,  B , '+', '-', 'neutrons' 
    case(3)
        print 1,  B , '-' , 'neutrons'
!    case(4)
!        print 1,  B , '-' ,'-', 'neutrons'
    case(5)        
        print 1,  B , '+' ,  'protons'
!    case(6)        
!        print 1,  B , '+' ,'-',  'protons'
    case(7)        
        print 1,  B , '-' ,  'protons'
!    case(8)        
!        print 1,  B , '-' ,'-',  'protons'
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
