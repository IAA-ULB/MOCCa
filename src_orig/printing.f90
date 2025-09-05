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
 ! Hephaestos keywords
 !
 ! TR        : $TR
 ! NTR       : $NTR
 !==============================================================================
 use geninfo
 use pairing
 use wavefunctions
 use convergence

 implicit none
 
contains

  subroutine PrintSpwfs(print_advanced,print_last)
    !---------------------------------------------------------------------------
    ! Print the info of the (physical) Hartree-Fock basis and the canonical
    ! basis in the case of HFB calculations.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !   print_advanced : logical, if .true. print ALL details of the spwfs
    !                             if .false., skip some properties
    !
    !   printlast      : integer, if non-zero, print only information on this
    !                             amount of highest-energy spwfs in each
    !                             isospin block
    !---------------------------------------------------------------------------

    10 format (42 ('-'), ' Hartree-Fock basis', 80('-'))
    13 format (42 ('-'), ' Canonical    basis', 80('-'))
    20 format (133 ('-'))
    30 format (133 ('_'),/,3x , 'Neutron wavefunctions')
    40 format (133 ('_'),/,3x , 'Proton  wavefunctions')
    60 format (3x,' n ', 4x, 'i', 4x,'P',4x, 'Rz', 3x,'occ',10x,'E',7x,       &
    &             'd2h',4x,'Delta', 1x,                                      &
    &             ' | ', 2x, 'JxT',4x, 'JyT', 4x,'Jz', 6x, 'J', 2x,          &
    &             ' | ', 2x, 'SxT',4x, 'SyT', 4x,'Sz', '   | r_rms ',        &
    &             ' | MPI_RANK ' )    

    11 format (1x, i5, 1x, i5, 1x, f5.2, 1x, f4.1, 2x, f6.4, 1x, a1, &
    &          1x, f9.3, 1x,es8.1,1x, f6.2,  1x,'|', 4(2x, f5.2), 1x, '|',   &
    &          3(2x, f5.2), ' | ', f6.2 , ' | ', i4)

    12 format (1x, i5, 1x, i5, 1x, f5.2, 1x, f4.1, 2x, f6.4, 1x, a1, &
    &          1x, f9.3, 1x,es8.1,1x, f6.2,  1x,'|', 4(3x, '*', 3x), 1x, '|',  &
    &          3(3x, '*', 3x), ' | ', 3x, '*', 2x , ' | ', i4)

    logical, intent(in) :: print_advanced
    integer, intent(in) :: print_last

    integer          :: wave,k, B, si, N, T, wavebar, l
    integer          :: ProtonOrder(nwp), NeutronOrder(nwn), sumocc
    real(KIND=dp)    :: p, Jx, Jy, Jz, JJ, s, Delta, Sx, Sy, Sz, r2
    character(len=1) :: blo
    real(KIND=dp), allocatable :: HF_gaps(:,:), can_gaps(:,:)

    ! We transform the gaps to the Hartree-Fock basis for printing
    if(pairingtype.eq.2) then
      si = 0
      allocate(HF_gaps(nwt,nwt)) ; HF_gaps = 0.0d0
      do B=1,8,2
        N = HFblocks_global(B) ; if(N.eq.0) cycle
        T = HFBlocks_global(B+1) + N

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

        if(print_last .ne. 0 .and. k .lt. (nwn - print_last) ) cycle

$NTR    sumocc = k
$TR     sumocc = 2*k

#if(PASTA == 0)
        P = P_hf(wave)
#else
        ! Temporary hack
        P = 0
#endif
        if(wave .le. sum(HFBlocks_global(1:2))) then
            if(wave .le. HFBlocks_global(1)) then
               s = +1
            else
               s = -1
            endif
        else
            if(wave .le. sum(HFBlocks_global(1:3))) then
               s = +1
            else
               s = -1
            endif
        endif

        if(print_advanced) then
          Jx = HF_JTR(1,wave) ; SX = HF_STR (1,wave)
          Jy = HF_JTI(2,wave) ; SY = HF_STI (2,wave)
          Jz = HF_J(3,wave)   ; SZ = HF_spin(3,wave)
          JJ = HF_JJ(wave)

          r2 = sqrt(spwf_r2_hf(wave))
          if(pairingtype.eq.1) then
            print 11, sumocc, wave, p, s, rho_can(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), BCSgaps(wave),                    &
            &               Jx, Jy, Jz, JJ, Sx, Sy, Sz, r2, rank_map(wave)
          elseif(pairingtype.eq.2) then
            print 11, sumocc, wave, p, s, rho_HF(wave), ' ', spenergies(wave),   &
            &               dispersions(wave), maxval(abs(HF_gaps(wave,:))),     &
            &               Jx, Jy, Jz, JJ, Sx, Sy, Sz, r2, rank_map(wave)
          else
            print 11, sumocc, wave, p, s, rho_can(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), 0.0, Jx, Jy, Jz, JJ, Sx, Sy, Sz,  &              
            &               r2, rank_map(wave)
          endif
        else
          if(pairingtype.eq.1) then
            print 12, sumocc, wave, p, s, rho_can(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), BCSgaps(wave), rank_map(wave)
          elseif(pairingtype.eq.2) then
            print 12, sumocc, wave, p, s, rho_HF(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), maxval(abs(HF_gaps(wave,:))), rank_map(wave)
          else
            print 12, sumocc, wave, p, s, rho_can(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), 0.0, rank_map(wave)
          endif
        endif
    enddo

    print 40  
    print 60
    print 20
    do k=1,nwp
        wave = ProtonOrder(k)

        if(print_last .ne. 0 .and. k .lt. (nwp - print_last) ) cycle

$NTR    sumocc = k
$TR     sumocc = 2*k


#if(PASTA == 0)
        P = P_hf(wave)
#else
        ! Temporary hack
        P = 0
#endif

        if(wave .le. sum(HFBlocks_global(1:6))) then
            if(wave .le. sum(HFBlocks_global(1:5))) then
               s = +1
            else
               s = -1
            endif
        else
            if(wave .le. sum(HFBlocks_global(1:7))) then
               s = +1
            else
               s = -1
            endif
        endif

        if(print_advanced) then
          Jx = HF_JTR(1,wave) ; SX = HF_STR (1,wave)
          Jy = HF_JTI(2,wave) ; SY = HF_STI (2,wave)
          Jz = HF_J(3,wave)   ; SZ = HF_spin(3,wave)
          JJ = HF_JJ(wave)

          r2 = sqrt(spwf_r2_hf( wave))
          if(pairingtype.eq.1) then
            print 11, sumocc, wave, p, s, rho_can(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), BCSgaps(wave),                    &
            &               Jx, Jy, Jz, JJ, Sx, Sy, Sz, r2, rank_map(wave)
          elseif(pairingtype.eq.2) then
            print 11, sumocc, wave, p, s, rho_HF(wave), ' ', spenergies(wave),   &
            &               dispersions(wave), maxval(abs(HF_gaps(wave,:))),     &
            &               Jx, Jy, Jz, JJ, Sx, Sy, Sz, r2, rank_map(wave)
          else
            print 11, sumocc, wave, p, s, rho_can(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), 0.0, Jx, Jy, Jz, JJ, Sx, Sy, Sz,  &              
            &               r2, rank_map(wave)
          endif
        else
          if(pairingtype.eq.1) then
            print 12, sumocc, wave, p, s, rho_can(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), BCSgaps(wave), rank_map(wave)
          elseif(pairingtype.eq.2) then
            print 12, sumocc, wave, p, s, rho_HF(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), maxval(abs(HF_gaps(wave,:))),rank_map(wave)
          else
            print 12, sumocc, wave, p, s, rho_can(wave), ' ', spenergies(wave),  &
            &               dispersions(wave), 0.0, rank_map(wave)
          endif
        endif
    enddo
    print 20
    if(allocated(HF_gaps)) deallocate( HF_gaps)  
    !---------------------------------------------------------------------------
    ! Return if we are not doing a HFB calculation
    if(PairingType.ne.2) return
    !---------------------------------------------------------------------------
    ! Otherwise, print the properties of the canonical basis
    print 13
    print 30
    print 60
    print 20

    ! Prepare by calculating the gaps in the canonical basis
    allocate(can_gaps(nwt,nwt))
    
    can_gaps = matmul(transpose(cantransfo), HFBgaps)
    can_gaps = matmul(can_gaps, cantransfo)

    ! Order the canonical basis, not the HF one
    ProtonOrder = OrderSpwfsISO(+1, .true.) 
    NeutronOrder= OrderSpwfsISO(-1, .true.)

    do k=1,nwn 
      wave = NeutronOrder(k) 

      if(print_last .ne. 0 .and. k .lt. (nwn - print_last) ) cycle

$NTR    sumocc = k
$TR     sumocc = 2*k

      !P = P_can(wave)

      if(wave .le. sum(HFBlocks_global(1:2))) then
          if(wave .le. HFBlocks_global(1)) then
             s = +1
          else
             s = -1
          endif
      else
          if(wave .le. sum(HFBlocks_global(1:3))) then
             s = +1
          else
             s = -1
          endif
      endif

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

      blo = ' ' 
      if(allocated(blocked_sps)) then
        do l = 1, blocknumber
          if(wave.eq.blocked_sps(l)) blo = '*'
        enddo
      endif    

      if(print_advanced) then
          if(allocated(canpsi)) then
            Jx = can_JTR(1,wave) ; SX = can_STR (1,wave)
            Jy = can_JTI(2,wave) ; SY = can_STI (2,wave)
            Jz = can_J(3,wave)   ; SZ = can_spin(3,wave)
            JJ = can_JJ(wave)
          !else
          !  Jx = spwf_JTR(1,wave) ; SX = spwf_STR (1,wave)
          !  Jy = spwf_JTI(2,wave) ; SY = spwf_STI (2,wave)
          !  Jz = spwf_J(3,wave)   ; SZ = spwf_spin(3,wave)
          !  JJ = spwf_JJ(wave)
          endif    

          r2 = sqrt(spwf_r2_can(wave))

          print 11, sumocc, wave, p,  s,   rho_can(wave), blo , canenergies(wave), &
          &              0.0, Delta , Jx, Jy, Jz, JJ, Sx, Sy, Sz, r2, rank_map(wave)
      else
          print 12, sumocc, wave, p,  s,   rho_can(wave), blo , canenergies(wave), &
          &              0.0, Delta, rank_map(wave)
      endif
    enddo
    print 40  
    print 60
    print 20

    do k=1,nwp 
      wave = ProtonOrder(k) 

      if(print_last .ne. 0 .and. k .lt. (nwp - print_last) ) cycle

$NTR    sumocc = k
$TR     sumocc = 2*k

      !P = P_can(wave)

      if(wave .le. sum(HFBlocks_global(1:6))) then
          if(wave .le. sum(HFBlocks_global(1:5))) then
             s = +1
          else
             s = -1
          endif
      else
          if(wave .le. sum(HFBlocks_global(1:7))) then
             s = +1
          else
             s = -1
          endif
      endif

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
    
      blo = ' ' 
      if(allocated(blocked_sps)) then
        do l = 1, blocknumber
          if(wave.eq.blocked_sps(l)) blo = '*'
        enddo
      endif   

      if(print_advanced) then
        if(allocated(canpsi)) then
          Jx = can_JTR(1,wave) ; SX = can_STR (1,wave)
          Jy = can_JTI(2,wave) ; SY = can_STI (2,wave)
          Jz = can_J(3,wave)   ; SZ = can_spin(3,wave)
          JJ = can_JJ(wave)
        !else
        !  Jx = spwf_JTR(1,wave) ; SX = spwf_STR (1,wave)
        !  Jy = spwf_JTI(2,wave) ; SY = spwf_STI (2,wave)
        !  Jz = spwf_J(3,wave)   ; SZ = spwf_spin(3,wave)
        !  JJ = spwf_JJ(wave)
        endif    
        r2 = sqrt(spwf_r2_can(wave))

        print 11, sumocc, wave, p, s,   rho_can(wave),  blo, canenergies(wave),  &
        &               0.0, Delta, Jx, Jy, Jz, JJ, Sx, Sy, Sz, r2, rank_map(wave)
      else
        print 12, sumocc, wave, p, s,   rho_can(wave),  blo, canenergies(wave),  &
        &               0.0, Delta, rank_map(wave)
      endif
    enddo
    print 20

    deallocate(can_gaps)
  end subroutine PrintSpwfs

  subroutine printqps(print_advanced)
    !---------------------------------------------------------------------------
    ! Print all relevant info on quasiparticles.
    ! - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Input:
    !   print_advanced : logical, if .false. does not print angular momenta.
    !                    This option exists because we cannot be sure that the 
    !                    underlying single-particle angular momentum expectation
    !                    values are up-to-date.
    !---------------------------------------------------------------------------
    logical, intent(in)  :: print_advanced
    integer              :: i, N, B, si, sb, ind, N2, k, U(1),V(1), T
    integer, allocatable :: indices(:)
    real(KIND=dp)        :: ov, v2, u2
    character(len=1)     :: Bstr, Pstr

    1  format (48 ('-'), 'Quasiparticles',50('-'))
    2  format ( i5, 1f10.2, 2x, 1es12.2, 1es12.2, ' | ', 2i4,  2x, 2f8.5,  &
    &           ' | ', 2x, a1,   &
    &           2x, a1,  2x, 1f5.3, ' | ',  3(2x,f5.2))
    3  format ( i5, 1f10.2, 2x, 1es12.2, 1es12.2, ' | ', 2i4,  2x, 2f8.5,  &
    &           ' | ', 2x, a1,   &
    &           2x, a1,  2x, 1f5.3, ' | ',  3(2x,'*', 2x))

    11  format(110('-'))
 
    ! Trash statement to stop the cray compiler complaining about non-allocated
    ! arrays because of the (possible) early return below.
    allocate(indices(1)) ; deallocate(indices)

    if(PairingType.eq.0) return
    if(PairingType.eq.2) call update_qp_angmom(Bogoliubov)

    print 1

    si = 0
    sb = 0
    do B=1,8,2
        N = HFblocks_global(B) ;      if(N.eq.0) cycle
        N2 = HFBlocks_global(B+1)
        T = N + N2
        ! ^------ these are all global indices, i.e. spanning all MPI ranks

        call print_qp_header(B)
        select case(pairingtype)
        case(2)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          ! HFB case
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The first half of the Bogoliubov matrix
          do i=1,T
            ! What are the single-particles dominating these qps? 
            if(maxval(abs(Bogoliubov(sb+1:sb+T,sb+i))).gt. 0.1) then
              U =  maxloc(Bogoliubov(sb+1:sb+T,sb+i)**2)+si
            else
              U = 0
            endif
            if(maxval(abs(Bogoliubov(sb+T+1:sb+2*T,sb+i))).gt. 0.1) then
              V =  maxloc(Bogoliubov(sb+T+1:sb+2*T,sb+i)**2)+si
            else
              V = 0
            endif  

            u2 = sum(Bogoliubov(sb  +1:sb+  T,sb+i)**2) 
            v2 = sum(Bogoliubov(sb+T+1:sb+2*T,sb+i)**2) 

            if(print_advanced) then
              print 2, i, QPenergies(sb+i), 1-configmatrix(sb+2*T-i+1),        &
              &           qpdispersions(sb+i),                                 &
              &           U(1), V(1), u2, v2, '-', '-', 0.0d0,                 &
              &            qp_JTR(1,sb+i), qp_JTI(2,sb+i), QP_J(3,sb+i)
            else
              print 3, i, QPenergies(sb+i), 1-configmatrix(sb+2*T-i+1),        &
              &           qpdispersions(sb+i),                                 &
              &           U(1), V(1), u2, v2, '-', '-', 0.0d0
            endif
          enddo
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          ! The second half of the Bogoliubov matrix
          do i=T+1,2*T
            ! What are the single-particles dominating these qps? 
            if(maxval(abs(Bogoliubov(sb     +1:sb+T,sb+i))).gt. 0.1) then
              U =  maxloc(Bogoliubov(sb     +1:sb+T,sb+i)**2)+si
            else
              U = 0
            endif
            if(maxval(abs(Bogoliubov(sb+T+1:sb+2*T,sb+i))).gt. 0.1) then
              V =  maxloc(Bogoliubov(sb+T+1:sb+2*T,sb+i)**2)+si
            else
              V = 0
            endif          

            u2 = sum(Bogoliubov(sb  +1:sb+  T,sb+i)**2) 
            v2 = sum(Bogoliubov(sb+T+1:sb+2*T,sb+i)**2) 
            ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            ! Check if this particular qp is a blocked, partner or not special
            if(allocated(blocked_qps)) then
              Bstr = '-'
              Pstr = '-'
              ov   = 0.0d0 
              do k=1,size(blocked_qps)
                if(N2.eq.0) then
                  if(blocked_qps(k) .eq. si+i-T) then
                    Bstr = 'B'
                    Pstr = 'P'
                    ov = 1.0
                  endif
                else
                  if(blocked_qps(k) .eq. si+i-T) then
                     Bstr = 'B'
                     Pstr = '-'
                     ov = partner_overlaps(k)
                  endif
                  if(partner_qps(k) .eq. si+i-T) then
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
            if(print_advanced) then
              print 2, i, QPenergies(sb+i), configmatrix(sb+i),                & 
              &           qpdispersions(sb+i),                                 &
              &           U(1), V(1), u2, v2, Bstr, Pstr, ov,                  &
              &           qp_JTR(1,sb+i), qp_JTI(2,sb+i), QP_J(3,sb+i)
            else
              print 3, i, QPenergies(sb+i), configmatrix(sb+i),                & 
              &           qpdispersions(sb+i),                                 &
              &           U(1), V(1), u2, v2, Bstr, Pstr, ov
            endif
          enddo

        case(1)
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
          ! BCS case
          ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
          allocate(indices(N))
          indices = order(BCSqps(si+1:si+N), N)
          do i=1, N
            ind  = indices(i)
            if(print_advanced) then
              print 2, i, BCSqps(si+ind), BCSf(si+ind),0.0d0, 0,0, 0.0d0,0.0d0,&
              &        '-', '-',0.0d0, 0.0d0, 0.0d0,0.0d0
            else
              print 3, i, BCSqps(si+ind), BCSf(si+ind),0.0d0, 0,0, 0.0d0,0.0d0,&
              &        '-', '-',0.0d0
            endif
          enddo
          deallocate(indices)
        end select
        si = si +   T
        sb = sb + 2*T
        print *
    enddo
    print 11
  end subroutine printqps

  subroutine print_qp_header(B)
    !---------------------------------------------------------------------------
    ! Print a header for the table of quasiparticle properties.
    !
    !---------------------------------------------------------------------------
    integer, intent(in) :: B  

    1  format ('Block ', i1, ':  P=',a1,'1',2x,  a8)
    2  format ( '  N      Eqp       f_n         disp     |   U   V', &
    &           '    sum u^2 sum v^2 |   B  P', '  ov_TR |'        , &
    &           4x,'JxT',4x,'JyT',4x,'Jz')
    3  format (110 ('_'))
  
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
end module
