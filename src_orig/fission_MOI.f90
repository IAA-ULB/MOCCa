module fission_MOI
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
 ! Module governing the calculation of the collective inertia tensor for 
 ! fission calculations. 
 !
 ! Currently implemented: 
 !  (i) perturbative cranking approximation
 !
 ! References used for the construction of this module:
 !  (a) S. A. Giuliani and L. M. Robledo,
 !      Non-perturbative collective inertias for fission: a comparative study
 !      Physics Letters B 787, 134-140 (2018).
 !
 !  (b) A. Baran, J. A. Sheikh, J. Dobaczewski, W. Nazarewicz and A. Staszczak
 !      Quadrupole collective inertia in nuclear fission: cranking approximation
 !      Phys. Rev. C 84, 054321 (2011)
 !
 !------------------------------------------------------------------------------
 ! Hephaestos keywords
 ! 
 ! PBROKEN : $PBROKEN
 !
 !------------------------------------------------------------------------------

  use geninfo
  use densities
  use parameterization
  use moments
  use timing

  implicit none

  !-----------------------------------------------------------------------------
  ! Multipole moments for which to construct the inertia tensor. 
  ! Hardcoded at the moment, maybe a runtime parameter in the future.
!  integer, parameter :: N_inertia            = 5
!  integer, parameter :: inertia_l(N_inertia) = (/1,2,2,3,4/) !,3/) 
!  integer, parameter :: inertia_m(N_inertia) = (/0,0,2,0,0/) !,0/)

  integer, parameter :: N_inertia            = 3
  integer, parameter :: inertia_l(N_inertia) = (/2,2,4/) !,3/) 
  integer, parameter :: inertia_m(N_inertia) = (/0,2,0/) !,0/)

  !-----------------------------------------------------------------------------
  ! Contains the full inertia tensor 
  real(KIND = dp), allocatable :: collective_inertia(:,:)
  
contains 

  subroutine print_collective_inertia()
    !---------------------------------------------------------------------------
    ! Print all entries in the collective inertia tensor.
    !---------------------------------------------------------------------------

    1 format (' Collective inertia tensor')
    2 format ('--------------------------')
    3 format ('     I_Q',2i1)
    4 format ('     I_Q',2i1, 1x,'|', 1x, 99f10.3)
    5 format (12('_'))

    character(len=80) :: header, sep
    character(len=12) :: tmp
    integer :: i
    
    print 2
    print 1
    print 2

    header = ''
    do i=1,N_inertia
			write(tmp, 3) inertia_l(i),inertia_m(i) 
			header = adjustl(trim(header)//tmp)
    enddo

    print *, '                 ', header   

    write(sep,5)
    tmp = ''
    do i=1,N_inertia
      write(tmp,5)    
      sep = adjustl(trim(sep)//tmp)  
    enddo
    print *, sep 
    do i=1, N_inertia
      print 4, inertia_l(i),inertia_m(i), collective_inertia(i,1:N_inertia)
    enddo
    print *,sep
    print *
        
  end subroutine print_collective_inertia

  subroutine calc_collective_inertia()
    !---------------------------------------------------------------------------
    ! Calculate the collective inertia tensor for the set of multipole moments
    ! Qlm determined in inertia_l and inertia_m. This tensor in the perturbative
    ! cranking approximation is given by
    !
    !     M_c = M_1^{-1} M_3 M_1^{-1}
    !
    ! where the matrices M_n are determined by
    !
    !                          Q^{20}_{i,ab} Q^{20}_{j,ab} 
    ! M_{n,ij} = Re sum_{ab}  ---------------------------
    !                               (E_a + E_b)^n
    !
    ! where the sum is over all quasiparticle states and E_a and E_b are 
    ! quasiparticle energies. 
    !
    ! Steps:
    !  (1) Calculate all the single-particle matrix elements of the Qlm
    !      in the HF-basis with routine Qlm_spme
    !  (2) Transform these matrix elements to the qp basis
    !      with routine calc_Q20
    !  (3) Sum the matrix elements, weighted with the appropriate power of 
    !      the quasiparticle energies, using Ksum_Mij
    !  (4) Invert M_1 with Lapack routines
    !  (5) Obtain M_c
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note: this routine is not yet ready to deal with blocked HFB vacua!
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable :: Mat(:,:,:,:), Qsp(:,:,:), Q20(:,:,:)
    real(KIND=dp), allocatable :: M1(:,:), M3(:,:), work(:)
    integer :: i, j, la, ma, lb, mb, l, m, info, lwork
    integer, allocatable :: ipiv(:)
        
    call start_timer(T_collective_MOI)
        
    if(.not.allocated(collective_inertia)) then
      allocate(collective_inertia(N_inertia, N_inertia))
    endif
    collective_inertia = 0
    
    allocate(Mat(N_inertia, N_inertia, 2,2)) ;  Mat   = 0.0d0
    allocate(Qsp(nwt,nwt,N_inertia))         ;  Qsp = 0.0d0
    
    if(pairingtype.eq.2) then
      allocate(Q20(nwt,nwt,N_inertia))         ;  Q20 = 0.0d0
    endif  
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Step 1 & 2: construct relevant sp matrices
    do i=1, N_inertia
      l = inertia_l(i)
      m = inertia_m(i)

      ! Calculate all relevant single-particle matrix elements              
      Qsp(:,:,i) = Qlm_spme(l,m,.false.) ! Hardcoded to consider only real parts
                                         ! at the moment

      if(pairingtype.eq.2) then
        ! Transform to the quasiparticle basis if needed
        Q20(:,:,i) = calc_Q20(Qsp(:,:,i), bogoliubov, l)
      endif
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Step 3: we calculate the sums for every combination of collective DOF
    do i=1,N_inertia
      la = inertia_l(i)
      do j=1, N_inertia
        lb = inertia_l(j)
        ! Perform the sums to obtain M_k for k=1,3
        select case(pairingtype)
        case(1)
          Mat(i,j,:,:) = Ksum_Mij_BCS(Qsp(:,:,i), Qsp(:,:,j), &
          &                                                      la, lb,(/1,3/))
        case(2)
          Mat(i,j,:,:) = Ksum_Mij(Q20(:,:,i), Q20(:,:,j), la, lb,  (/1,3/))
        case DEFAULT
          print *, 'NOT IMPLEMENTED.'
          stop
        end select
      enddo
    enddo

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Constructing explicitly the matrices M_1 and M_3 by summing proton and 
    ! neutron contributions    
    allocate(M1(N_inertia, N_inertia)) 
    allocate(M3(N_inertia, N_inertia)) 
    M1 = Mat(:,:,1,1) + Mat(:,:,1,2) 
    M3 = Mat(:,:,2,1) + Mat(:,:,2,2) 

    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Step 4: use LAPACK routines to invert M1
    ! 
    ! Ask for a workspace size
    allocate(work(1), ipiv(N_inertia))
    call dsytrf('U', N_inertia, M1, N_inertia, ipiv, work, -1, info)
    lwork = work(1)
    deallocate(work)
    allocate(work(lwork))
    ! Factorize M1
    call dsytrf('U', N_inertia, M1, N_inertia, ipiv, work, lwork, info)
    deallocate(work)
    ! Invert M1
    allocate(Work(N_inertia))
    call dsytri('U', N_inertia, M1, N_inertia,ipiv,work, info)

    if(info.ne.0) then
       print *, 'Problem for DSYTRI during the calculation of collective inertia.'
       print *, 'INFO = ', info
       stop
    endif
    ! Note that after DSYTRI, only the top half of M1 is guaranteed to be right
    ! Thus, we populate the other half here to avoid any surprises
    do i=1,N_inertia
      do j=i+1,N_inertia
        M1(j,i) = M1(i,j)
      enddo
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Step 5: calculate cranking tensor 
    !  M_c = M1^{-1} M3 M1^{-1}
    collective_inertia = matmul(matmul(M1, M3), M1)

    call stop_timer(T_collective_MOI)

  end subroutine calc_collective_inertia
  
  function Ksum_Mij(Qa, Qb, la, lb,  Ks) result(Ksum)
    !---------------------------------------------------------------------------
    ! Perform the relevant sums over the quasiparticle space for the calculation
    ! of the collective inertia, i.e.
    ! 
    !  M_k = sum_ij (<0| Qa | ij >< ij | Qb | 0 >)/[(Ei + Ej)^k]
    !
    ! This routine bunches the summations for the same multipole moments with
    ! all different powers k.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     Qa, Qb    : two-quasiparticle representations of both multipole 
    !                 operators
    !     Ks        : set of powers to use in the inverted calculation
    ! Output:
    !     Ksum      : result of the summations, array with the size of Ks
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Qa(:,:), Qb(:,:)
    real(KIND=dp), allocatable:: Ksum(:,:)
    integer, intent(in)       :: Ks(:), la, lb
    
    integer :: Na, N2a, Ta, Ba, Bb, Tb, Nb, N2b, ita, itb, k, i,j, Nk
    integer :: sai, sbi, sab, sbb
    real(KIND=dp) :: num, denom

    Nk = size(Ks)
    allocate(Ksum(Nk,2)) ; Ksum = 0.0d0

    ! If parity is conserved, there is a parity selection rule    
$PCONSERVED if(mod(la,2) .ne. mod(lb,2)) return     

    sai = 0 ; sab = 0
    do Ba=1,8,2
      Na = HFBlocks(Ba) ; if(Na.eq.0) cycle
      N2a= HFBlocks(Ba+1)
      Ta = Na + N2a
      ita = 1 ; if(Ba.gt.4) ita=2
      
      sbi = 0 ; sbb = 0
      do Bb=1,8,2
        Nb = HFBlocks(Bb) ; if(Nb.eq.0) cycle
        N2b= HFBlocks(Bb+1)
        Tb = Nb + N2b
        itb= 1 ; if(Bb.gt.4) itb=2
        
        !Gain some CPU time
        if(ita.ne.itb) then
          sbi = sbi +   Tb
          sbb = sbb + 2*Tb
          cycle
        endif
        
        do i=1,Ta
          do j=1,Tb
            num = Qa(sai+i,sbi+j) * Qb(sai+i,sbi+j)
            do k=1,Nk
              denom = (qpenergies(sab+Ta+i) + qpenergies(sbb+Tb+j))**Ks(k)
              Ksum(k,ita) = Ksum(k,ita) + num/denom
!              if(k.eq.1 .and. abs(num) .gt.1d-6 )  then
!                print *,Ba, Bb, i,j, num/denom, qpenergies(sab+Ta+i),qpenergies(sbb+Tb+j)
!              endif
            enddo
          enddo
        enddo
        sbi = sbi +   Tb
        sbb = sbb + 2*Tb
      enddo
      sai = sai +   Ta
      sab = sab + 2*Ta
    enddo
    
    ! Factor two for the time-reversal partners
    Ksum = 2*Ksum
    
  end function Ksum_Mij
    
  function Ksum_Mij_BCS(Qa, Qb, la, lb,  Ks) result(Ksum)
    !---------------------------------------------------------------------------
    ! Perform the relevant sums over the quasiparticle space for the calculation
    ! of the collective inertia in the case of a BCS calculation, i.e.
    ! 
    !  M_k = sum_ij (<i|Qa|j><j|Q^\dagger_b|i >)/[(Ei + Ej)^k] {eta^+_ij}^2
    !
    ! where 
    !  (*) the <|Q|> are single-particle matrix elements
    !  (*) the Ei and Ej are a BCS quasiparticle energies
    !  (*) eta^+_ij = u_i v_j + u_j v_i 
    !  (*) the sum runs over all single-particle states
    !
    ! This expression is Eq. 58 in 
    ! 
    !  A. Baran et al,  Phys. Rev. C 84, 054321 (2011).
    !   
    ! This routine bunches the summations for the same multipole moments with
    ! all different powers k.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     Qa, Qb    : two-quasiparticle representations of both multipole 
    !                 operators
    !     Ks        : set of powers to use in the inverted calculation
    ! Output:
    !     Ksum      : result of the summations, array with the size of Ks
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Qa(:,:), Qb(:,:)
    real(KIND=dp), allocatable:: Ksum(:,:)
    integer, intent(in)       :: Ks(:), la, lb
    real(KIND=dp)             :: num, denom, eta, ui, vi, uj, vj
    integer                   :: i,j, it, itb, k, Nk

    Nk = size(Ks)
    allocate(Ksum(Nk,2)) ; Ksum = 0.0d0

    ! If parity is conserved, there is a parity selection rule    
$PCONSERVED if(mod(la,2) .ne. mod(lb,2)) return     

    do i=1,nwt
      call uv_from_occupation(BCSoccupations(i), ui, vi)
     
      it = 1 ;  if(i.gt. nwn) it = 2
      do j=1,nwt
        call uv_from_occupation(BCSoccupations(j), uj, vj)
        
        itb = 1 ;  if(j.gt. nwn) itb = 2
        if(it .ne. itb) cycle

        eta   = ui * vj + vi * uj  
        num   = Qa(i,j) * Qb(i,j) * eta**2
        do k=1,Nk
         denom      = (BCSqps(i) + BCSqps(j))**Ks(k)
         Ksum(k,it) = Ksum(k,it) + num/denom
        enddo
      enddo
    enddo    
    ! Factor two for the time-reversal partners
    Ksum = 2*Ksum
    
  end function Ksum_Mij_BCS
  
  function calc_Q20(Qsp, bogo, l) result(Q20)
    !---------------------------------------------------------------------------
    ! Function that calculates the two-quasiparticle matrix for a multipole 
    ! moment operator in the case of a HFB calculation.
    ! 
    !  Q20 = U^\dagger Q V^* - V^\dagger Q^t U^*
    !
    ! where Q is the matrix of single-particle matrix elements of the 
    ! multipole operator.
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     Qsp       : single-particle matrix elements 
    !     Bogo      : Bogoliubov transformation from the sp basis to the 
    !                 quasiparticle basis
    !     l         : ell of the multipole moment
    ! Output:
    !     Q20       : two-quasiparticle component of Qlm
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note: 
    !
    !  * No attempt is made to take into account the effect of 
    !    (non-)conservation of time-reversal symmetry. This is not an issue 
    !    as long as we only consider real parts of multipole moments with 
    !    even m, on the condition we take the 'default' orientation of the 
    !    multipole moments. 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: Qsp(nwt,nwt), bogo(2*nwt,2*nwt)
    integer, intent(in)        :: l
    real(KIND=dp)              :: Q20(nwt,nwt), tmp(nwt,nwt)
    real(KIND=dp), allocatable :: U(:,:), V(:,:)    
    integer                    :: N, N2, T, si, sb, B

    Q20 = 0.0d0

    !---------------------------------------------------------------------------
    ! This situation is:
    !  (a) it concerns a multipole moment with even parity; or
    !  (b) parity is broken and there are no more parity blocks
    ! In both cases we can simply multiply matrices straightforwardly.
$PCONSERVED    if(mod(l,2).eq.0) then
      si = 0 ; sb = 0
      do B=1,8,2
        N = HFBlocks(B) ; if(N.eq.0) cycle
        N2= HFBlocks(B+1)
        T = N + N2

        U = bogo(sb+  1:sb+  T,sb+ T+1:sb+2*T)
        V = bogo(sb+T+1:sb+2*T,sb+ T+1:sb+2*T)
        
        Q20(si+1:si+T,si+1:si+T) = matmul(transpose(U),Qsp(si+1:si+T,si+1:si+T))
        ! U^dagger Q V^*
        Q20(si+1:si+T,si+1:si+T) = matmul(Q20(si+1:si+T,si+1:si+T), V)
        
        ! Q V^* 
        tmp(si+1:si+T,si+1:si+T) = matmul(Qsp(si+1:si+T,si+1:si+T), V)
        ! V^\dagger Q^t
        tmp(si+1:si+T,si+1:si+T) = transpose(tmp(si+1:si+T,si+1:si+T))
        ! V^dagger Q^t U^*
        tmp(si+1:si+T,si+1:si+T) = matmul(tmp(si+1:si+T,si+1:si+T), U)

        ! Put both parts together      
        !
        ! ATTENTION TO THE SIGN DUE TO TIMEREVERSAL
        !
        Q20(si+1:si+T,si+1:si+T) = - Q20(si+1:si+T,si+1:si+T) &
        &                          - tmp(si+1:si+T,si+1:si+T)

        si = si +   N +   N2
        sb = sb + 2*N + 2*N2
      enddo
$PCONSERVED    endif
    
    !---------------------------------------------------------------------------
    ! The following situation concerns only the case for a multipole moment 
    ! with odd l and parity is conserved.
    !
    ! Using straightforward notation, this case corresponds to
    !
    !  U = ( U+ 0 )  V = ( V+ 0 )  Q = (0   Q+-)
    !      ( 0  U-)      ( 0  V-)      (Q+- 0  )
    !
    ! And so
    !
    ! U^dagger Q V^* = (  0              U+^dagger Q V-  )
    !                  ( U-^dagger Q V+          0       )
    if(mod(l,2).eq.1) then
!      do B=1,8,4 ! Just an isospin loop ...;
!      
!      enddo
      ![TO BE IMPLEMENTED]
!      stop
    endif

  end function calc_Q20

  function Qlm_spme(l,m,Imaginary) result(me)
    !---------------------------------------------------------------------------
    ! Routine to calculate single-particle matrix elements of a given multipole
    ! moment in the (a) HARTREE-FOCK basis and (b) in units of b^(ell/2) with
    ! 1 b = 100 fm^2.
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !     l,m       : integers, 
    !     Imaginary : to calculate < i | Re Q_lm | j > or < i |Im Q_lm | j >
    !                 logical
    ! Output:
    !     me        : single-particle matrix elements of the multipole moment
    !                 in the canonical basis, real(nwt,nwt), in units of 
    !                 b^(ell/2).
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Notes:
    ! *) does not allow for the calculation of matrix elements with odd values 
    !    of m yet; those connect single-particle states with different 
    !    signature, for the traditional definition of spherical harmonics at 
    !    least.
    !---------------------------------------------------------------------------

    integer, intent(in)        :: l, m
    logical, intent(in)        :: Imaginary
    real(KIND=dp)              :: me(nwt,nwt)
    
    real(KIND=dp), allocatable         :: SpherHarmMesh(:,:,:,:,:,:)
    real(KIND=dp), allocatable, target :: Qlm(:)
    real(KIND=dp), pointer             :: harm_3D(:,:,:)
    
    integer :: i, j, Bi, Bj, Ni, Nj, si, sj, im, k

    if(mod(m,2) .eq. 1) then
      print *, 'Qlm_spme does not allow for the calculation of multipole matrix elements with odd m YET.'
      stop
    endif
    if(Imaginary) then 
      print *, 'Qlm_spme does not allow for the calculation of imaginary multipole matrix elements YET.'
      stop
    endif
    if((m .lt. 0) .or. (m .gt. l)) then
      print *, 'Invalid value of m in Qlm_spe.'      
      stop
    endif
    
    ! Initialize to zero
    me = 0.0d0

    ! Generate the necessary spherical harmonic
    allocate(SpherHarmMesh(nx,ny,nz,0:l,0:l,2), Qlm(nx*ny*nz))
    call GenSphericalHarmonics(l,nx,ny,nz, &
    &                          meshx_shifted,meshy_shifted,meshz_shifted,    & 
    &                          SpherHarmMesh,quantisationaxis,secondaryaxis) 

    if(Imaginary) then
      im  = 2
    else
      im  = 1
    endif
    
    harm_3D(1:nx,1:ny,1:nz) => Qlm
    harm_3D                 = SpherHarmMesh(:,:,:,l,m,im)
    deallocate(SpherHarmMesh)
      
    !--------------------------------------------------------------------------- 
    ! Loop over the neutron single-particle states
    !---------------------------------------------------------------------------
    si = 0
    do Bi = 1, 4
      Ni =  HFBlocks(Bi) ; if(Ni.eq.0) cycle
      do Bj = 1, 4
        Nj = HFBlocks(Bj) ; if(Nj.eq.0) cycle

        if(bj.eq.1) then
          sj = 0
        else
          sj = sum(HFblocks(1:Bj-1))
        endif

        ! Parity selection rule: if l = even, then only single-particle states
        !                        of identical parity (and signature) contribute    
$PCONSERVED        if(mod(l,2) .eq. 0 .and.  Bi .ne. Bj)        cycle
        ! Parity selection rule: if l = odd, then only single-particle states
        !                        of different parity contribute
$PCONSERVED        if(mod(l,2) .eq. 1 .and. abs(Bi - Bj).ne.2 ) cycle
        ! (THIS is only needed if parity conserved of course)
$PBROKEN if( Bi .ne. Bj ) cycle

        do i=1,NI    
         do j=1,NJ
          ! me = Int d^3r Sum_sigma psi^*_i(r,sigma) psi_j(r,sigma) Qlm(r)
          me(si+i,sj+j) = 0
          do k=1,4
            me(si+i,sj+j) = me(si+i,sj+j) & 
            &                    + sum(HFpsi(:,k,si+i)*HFpsi(:,k,sj+j)*Qlm(:))
          enddo
          ! All these matrix elements are real if 
          ! (i)  we consider only real multipole moments
          ! (ii) time simplex is conserved
          ! 
          ! which means the matrix we store them in is symmetric
          me(si+i,sj+j) = me(si+i,sj+j) * dv 
          me(sj+j,si+i) = me(si+i,sj+j)
         enddo
        enddo
      enddo
      si = si + NI
    enddo

    !--------------------------------------------------------------------------- 
    ! Loop over the proton single-particle states
    !---------------------------------------------------------------------------
    si = nwn
    do Bi = 5, 8
      Ni =  HFBlocks(Bi) ; if(Ni.eq.0) cycle
      do Bj = 5, 8
        Nj = HFBlocks(Bj) ; if(Nj.eq.0) cycle
        sj = sum(HFblocks(1:Bj-1))

        ! Parity selection rule: if l = even, then only single-particle states
        !                        of identical parity (and signature) contribute    
$PCONSERVED        if(mod(l,2) .eq. 0 .and.  Bi .ne. Bj)        cycle
        ! Parity selection rule: if l = odd, then only single-particle states
        !                        of different parity contribute
$PCONSERVED        if(mod(l,2) .eq. 1 .and. abs(Bi - Bj).ne.2 ) cycle
        ! (THIS is only needed if parity conserved of course)
$PBROKEN if( Bi .ne. Bj ) cycle

        do i=1,NI    
         do j=1,NJ
          ! me = Int d^3r Sum_sigma psi^*_i(r,sigma) psi_j(r,sigma) Qlm(r)
          me(si+i,sj+j) = 0
          do k=1,4
            me(si+i,sj+j) = me(si+i,sj+j) & 
            &                    + sum(HFpsi(:,k,si+i)*HFpsi(:,k,sj+j)*Qlm(:))
          enddo
          ! All these matrix elements are real if 
          ! (i)  we consider only real multipole moments
          ! (ii) time simplex is conserved
          ! 
          ! which means the matrix we store them in is symmetric
          me(si+i,sj+j) = me(si+i,sj+j) * dv 
          me(sj+j,si+i) = me(si+i,sj+j)
         enddo
        enddo
      enddo
      si = si + NI
    enddo
    
    ! Rescale with the units of b^(ell/2) with 1 b = 100 fm^2.
    me = me/(100**(l/2.0))
      
  end function Qlm_spme

end module fission_MOI

! Code zoo
    
!    allocate(coll_me(nwt,nwt))
!    
!    coll_me = Qlm_spme(2,0,.false.)
!    
!    s = 0 
!    do i=1, nwt 
!      do j=1,nwt
!        s = s + rho_pairing(i,j) * coll_me(j,i) * 2
!      enddo
!    enddo
!    print *, 'Collective Q20', s
!    
!    coll_me = Qlm_spme(3,0,.false.)
!    
!    s = 0 
!    do i=1, nwt 
!      do j=1,nwt
!        s = s + rho_pairing(i,j) * coll_me(j,i) * 2
!      enddo
!    enddo
!    print *, 'Collective Q30', s
!    
!    coll_me = Qlm_spme(4,0,.false.)
!    
!    s = 0 
!    do i=1, nwt 
!      do j=1,nwt
!        s = s + rho_pairing(i,j) * coll_me(j,i) * 2
!      enddo
!    enddo
!    print *, 'Collective Q40', s
!!    
!    
