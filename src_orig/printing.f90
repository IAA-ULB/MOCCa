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
    
    10 format (21 ('-'), ' Sp wavefunctions ', 41('-'))
    12 format (21 ('-'), ' Canonical basis  ', 41('-'))
    20 format (80 ('-'))
    30 format (80 ('_'),/,3x , 'Neutron wavefunctions')
    40 format (80 ('_'),/,3x , 'Proton  wavefunctions')
    60 format (2x,'i',4x,'P',3x, 'Rz', 3x,'occ',7x,'E',8x,'d2h',4x,'Delta',3x, &
    &             'JxT',4x, 'JyT', 4x ,'Jz', 5x, 'J')    

    11 format (i3, 1x, f4.1, 1x, f4.1, 2x, f6.4, 1x, f9.3, 1x, es8.1,1x,f6.2, 4(2x, f5.2))

    integer       :: wave,k, l, i, j, B, si, N, T, wavebar
    integer       :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp) :: p, Jx, Jy, Jz, JJ, s, Delta
    
    real(KIND=dp), allocatable :: HF_J(:,:),  HF_JTR(:,:), HF_JTI(:,:)
    real(KIND=dp), allocatable :: HF_J2(:,:), HF_JJ(:), HF_rho(:)
    real(KIND=dp), allocatable :: HF_gaps(:,:), can_gaps(:,:)
    !---------------------------------------------------------------------------
    ! The single-particle wavefunctions in storage are not necessarily 
    ! converging to the physical Hartree-Fock basis. Hence, we need to convert 
    ! the observables before printing.
    
    ! We start with angular momentum observables. 
    ! I have written out the matrix multiplications explicitly, as we need only
    ! the diagonal matrix elements. 
    allocate(HF_J(3,nwt), HF_JTR(3,nwt), HF_JTI(3,nwt), HF_rho(nwt))
    allocate(HF_J2(3,nwt), HF_JJ(nwt))
    
    si = 0
    do B=1,8
      N = HFBlocks(B); if(N.eq.0) cycle
      do k=1,3
        HF_J  (k,si+1:si+N) = 0.0
        HF_J2 (k,si+1:si+N) = 0.0
        HF_JTR(k,si+1:si+N) = 0.0
        HF_JTI(k,si+1:si+N) = 0.0
        do i=si+1,si+N
          do j=si+1,si+N
            do l=si+1,si+N
              HF_J  (k,i) = HF_J  (k,i) &
              &           + HFtransfo(l,i) * spwf_J  (k,l,j) * HFtransfo(j,i)                   
              HF_J2 (k,i) = HF_J2 (k,i) &
              &           + HFtransfo(l,i) * spwf_J2 (k,l,j) * HFtransfo(j,i)                   
              HF_JTR(k,i) = HF_JTR(k,i) &
              &           + HFtransfo(l,i) * spwf_JTR(k,l,j) * HFtransfo(j,i)                   
              HF_JTI(k,i) = HF_JTI(k,i) &
              &           + HFtransfo(l,i) * spwf_JTI(k,l,j) * HFtransfo(j,i)                   
            enddo
          enddo
        enddo
      enddo
      si = si + N
    enddo
    ! The JJ variable is still simple
    do wave=1,nwt
      HF_JJ(wave) = (-1. + sqrt(1. + 4*sum(HF_J2(:,wave))))/2.
    enddo

    ! In the case of HFB pairing, we would like the diagonal matrix elements
    ! of the density as well.
    if(pairingtype.eq.2) then
      HF_rho = 0.0d0
      do k=1,nwt  
        HF_rho(k) = 0.0d0
        do l=1,nwt
          do j=1,nwt
            HF_rho(k) = HF_rho(k) + &
            &                 HFtransfo(l,k) * rho_pairing(l,j) * HFtransfo(j,k)
          enddo
        enddo
      enddo
    endif
    
    ! finally, we transform the gaps to the Hartree-Fock basis for printing
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

        Jx = HF_JTR(1,wave)
        Jy = HF_JTI(2,wave)
        Jz = HF_J(3,wave)
        JJ = HF_JJ(wave)

        if(pairingtype.eq.1) then
          print 11, wave, p, s, rho_can(wave), spenergies(wave), &
          &               dispersions(wave), BCSgaps(wave),Jx, Jy, Jz, JJ
        elseif(pairingtype.eq.2) then
          print 11, wave, p, s, HF_rho(wave), spenergies(wave),                &
          &               dispersions(wave), maxval(abs(HF_gaps(wave,:))),     &
          &               Jx, Jy, Jz, JJ
        else
          print 11, wave, p, s, rho_can(wave), spenergies(wave), &
          &               dispersion(wave), 0.0, Jx, Jy, Jz, JJ                
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

        Jx = HF_JTR(1,wave)
        Jy = HF_JTI(2,wave)
        Jz = HF_J(3,wave)
        JJ = HF_JJ(wave)

        if(pairingtype.eq.1) then
          print 11, wave, p, s, rho_can(wave), spenergies(wave),               &
          &               dispersions(wave), BCSgaps(wave),Jx, Jy, Jz, JJ
        elseif(pairingtype.eq.2) then
          print 11, wave, p, s,  HF_rho(wave), spenergies(wave),               &
          &               dispersions(wave), maxval(abs(HF_gaps(wave,:))),     &
          &               Jx, Jy, Jz, JJ
        else
          print 11, wave, p, s, rho_can(wave), spenergies(wave),               &
          &               dispersions(wave), 0.0, Jx, Jy, Jz, JJ
        endif
    enddo
    print 20
    deallocate(HF_J, HF_JTR, HF_JTI, HF_rho, HF_J2, HF_JJ, HF_gaps)  
    !---------------------------------------------------------------------------
    ! Return if we are not doing a HFB calculation
    if(PairingType.ne.2) return
    !---------------------------------------------------------------------------
    ! Otherwise, print the properties of the canonical basis
    print 12
    print 30
    print 60
    print 10

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

      Jx = can_JTR(1,wave)
      Jy = can_JTI(2,wave)
      Jz = can_J(3,wave)
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
      &               0.0, Delta , Jx, Jy, Jz, JJ
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
      &               0.0, Delta, Jx, Jy, Jz, JJ
    enddo
    print 20

  end subroutine PrintSpwfs

  subroutine printqps
    !---------------------------------------------------------------------------
    ! Print all relevant info on quasiparticles.
    ! Very bare-bones for the moment.
    !---------------------------------------------------------------------------
    integer              :: i, N, B, si, sb, ind, N2
    integer, allocatable :: indices(:)
    
    1  format (33 ('-'), 'Quasiparticles',33('-'))
    2  format ( i3, 1f10.2, 2x, 1es12.2)

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
          do i=1,N+N2
            print 2, i, QPenergies(sb+i), 1-configmatrix(sb+2*N+2*N2-i+1)
          enddo
          do i=N+N2+1,2*N+2*N2
            print 2, i, QPenergies(sb+i), configmatrix(sb+i)
          enddo

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

    1  format ('Block ', i1, ':  P=',a1,'1',2x,  a8)
    2  format ( '  N      Eqp     f_n')
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
