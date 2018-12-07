module HFB
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
 !
 ! Module that solves the HFB problem.
 !
 !==============================================================================
 !
 ! TODO 
 !   * Think about time-reversal, it will be some index juggling to get correct.
 !
 !==============================================================================
  use diag
  use geninfo
  use wavefunctions
  use pairingcutoffs
  
  implicit none
  
  !-----------------------------------------------------------------------------
  ! Gaps and quasiparticle energies
  real(KIND=dp), allocatable :: HFBGaps(:,:)

  !-----------------------------------------------------------------------------
  ! Maximum amount of iterations for finding a Fermi energy
  integer :: maxHFBiter  = 200
  integer :: HFBsizes(4) = 0
  !------------------------------------------------------------------------------
  real(KIND=dp) :: HFBdispersion(2)
  
  procedure(delta_action_dummy), pointer :: delta_action_HFB

  real(KIND=dp)  :: Fermiprec = 1d-9

contains

 subroutine initHFB
  !-----------------------------------------------------------------------------
  ! Initialize things in this module.
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! As the structure of the matrices depends significantly on the symmetries 
  ! assumed, the size of the relevant matrices. 
  !
  !-----------------------------------------------------------------------------
  
  ! The factors two are due to time-reversal.
  HFBSizes(1) = HFblocks(1)
  HFBSizes(2) = HFblocks(3)
  HFBSizes(3) = HFblocks(5)
  HFBSizes(4) = HFblocks(7)
  
  ! Otherwise the sizes of the blocks are simply the same as the HFBlocks
  
 end subroutine initHFB

 subroutine solvepairing_HFB(fermi, rho_pairing, kappa_pairing, configmatrix,  &
 &                         qpenergies,HFBmix, HFBmixtype,BlockType,Blockindices)
  !-----------------------------------------------------------------------------
  ! Driver routine for the solving of the HFB equations.
  !
  ! |----
  ! | Do until the number of particles is precise enough
  ! |   1) Diagonalize the HFB Hamiltonian for given lambda
  ! |   2) Construct rho_pairing 
  ! |   3) Calculate number of particles
  ! |   4) Update Lambda
  ! |----
  !     5) Return rho_pairing and kappa_pairing
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------  
  ! History of the pairing matrices
  real(KIND=dp), allocatable ::  rho_history(:,:), kappa_history(:,:)
  real(KIND=dp), allocatable ::  configmatrix_history(:)

  ! Options for the mixing of the HFB configurations
  real(KIND=dp), intent(in)  :: HFBmix
  integer, intent(in)        :: HFBmixtype


  !-----------------------------------------------------------------------------
  real(KIND=dp) :: configmatrix(2*nwt), c

  ! Configuration for the blocking
  integer, intent(in)        :: Blockindices(:)
  integer, intent(in)        :: BlockType

  real(KIND=dp)              :: fermi(2)
  real(KIND=dp)              :: rho_pairing(nwt,nwt), qpenergies(nwt)
  real(KIND=dp)              :: kappa_pairing(nwt,nwt)
  real(KIND=dp)              :: df(2), dn(2,2)

  real(KIND=dp)              :: sphamil(nwt,nwt),vect(2*nwt,2*nwt)
  real(KIND=dp)              :: eigen(2*nwt),  work(nwt), tempp
  real(KIND=dp)              :: particles(2), temp(nwt), chi(nwt,nwt)
  real(KIND=dp), allocatable :: HFBHamil(:,:)
  
  logical                    :: converged(2)
  
  integer                    :: si, sb, N, B, iter, wave1, it, i,j,k
  
  !-----------------------------------------------------------------------------
  if(.not.allocated(rho_history)) then
        allocate(rho_history(nwt,nwt))     ; rho_history   = 0.0
  endif
  if(.not. allocated(kappa_history)) then
        allocate(kappa_history(nwt,nwt))   ; kappa_history = 0.0
  endif  
  if(.not.allocated(configmatrix_history)) then
        allocate(configmatrix_history(2*nwt)) ;  configmatrix_history = 0.0
  endif
  !-----------------------------------------------------------------------------
  ! Saving the history
  rho_history     = rho_pairing
  kappa_history   = kappa_pairing
  configmatrix_history = configmatrix

  dn = 0.0
  
  ! Guess a new Fermi energy if none is there
  if(all(Fermi.eq.0.0))   Fermi = -5
  !-----------------------------------------------------------------------------
  ! Construct the single-particle hamiltonian from the sp.energes
  sphamil = 0
  si      = 0
  do B=1,8
    N = HFBlocks(B)
    do wave1=1,N
      sphamil(si+wave1,si+wave1) = spenergies(si+wave1)
    enddo
    si = si +   N
  enddo
  
  df = 0.0
        
  !-----------------------------------------------------------------------------
  ! Start iterations over the Fermi energy
  do iter=1,maxHFBiter
     vect      = 0.0
     eigen     = 0.0
     particles = 0.0
     !--------------------------------------------------------------------------
     ! Diagonalize every block. 
     si = 0 ; sb = 0
     do B=1,4
        
        it = 1
        if(B>2) it = 2
     
        N = HFBsizes(B)
        allocate(HFBHamil(2*N,2*N)) 
        HFBHamil = ConstructHFBHamil(sphamil(si+1:si+N,si+1:si+N),       &
        &                         HFBgaps(si+1:si+N,si+1:si+N), Fermi(it))
  
        call diagon (HFBHamil,2*N,2*N,vect(sb+1:sb+2*N, sb+1:sb+2*N),    &
        &                                         eigen(sb+1:sb+2*N),work)
        deallocate(HFBHamil)

!
!       Vect now contains the eigenvectors of the HFB hamiltonian in the 
!       given parity-isospin block, ordered by increasing E_qp.
!

!        !-----------------------------------------------------------------------
!        ! Switching half of the eigenvectors
!
!        This is not needed in a time-reversal invariant code.
!         
!        eigen(sb+1:sb+N) = -eigen(sb+1:sb+N)
!        do i = 1,N
!          temp                      = vect(sb+1  :sb+N,   sb+i)
!          vect(sb  +1:sb+N,   sb+i) = vect(sb+N+1:sb+2*N, sb+i)
!          vect(sb+N+1:sb+2*N, sb+i) = temp 
!        enddo
!        
        !-----------------------------------------------------------------------
        ! Saving the quasiparticle excitation energies
        qpenergies(si+1:si+N) = eigen(sb+N+1:sb+2*N)

        ! Startindex (si) for the next block.
        si = si +  N
        sb = sb +2*N
     enddo
     !--------------------------------------------------------------------------
     ! Construct the configuration matrix C
     configmatrix=ConstructConfiguration(Vect,QPenergies,BlockType,BlockIndices)
     !--------------------------------------------------------------------------

     si = 0 ; sb = 0
     do B=1,4                
        it = 1
        if(B>2) it = 2
        N = HFBsizes(B)
        !-----------------------------------------------------------------------
        ! Calculate the number of particles in here  
        ! Sum_i rho_ii =  Sum_ii   U   f U^{\dagger} + V^{*}(1 - f)V^{T}
        !            sum_(ij>N) f_(j) V^*_ij V^T_ji = sum_ij f_(j) V^*_ij V_ij
        !          + sum_(ij<N) f_(j) U^*_ij U^T_ji = sum_ij f_(j) U^*_ij U_ij  
        
!        print *, 'B', B, N, configmatrix(sb+1:sb+2*N)

        do i=1,N
            particles(it) = particles(it)                                      &
            &     + 2*configmatrix(sb+N+i) * sum(vect(sb+N+1:sb+2*N, sb+N+i)**2)            
            ! Time reversal is responsible for the factor 2
        enddo
        do i=1,N
            particles(it) = particles(it)                                      &
            &     + 2*configmatrix(sb  +i) * sum(vect(sb  +1:sb+  N, sb+N+i)**2)            
            ! Time reversal is responsible for the factor 2
        enddo
        
        !-----------------------------------------------------------------------
        ! Startindex (si) for the next block.
        si = si +  N
        sb = sb +2*N
     enddo
     !--------------------------------------------------------------------------
     ! Adjust Fermi energy, using a secant method for the moment
     dn(:,2) = dn(:,1)
     
     dn(1,1) = particles(1) - neutrons
     dn(2,1) = particles(2) - protons

     converged = .false.
     if(abs(dn(1,1)) .lt. FermiPrec) then
      converged(1) = .true.
     endif
     if(abs(dn(2,1)) .lt. FermiPrec) then
      converged(2) = .true.
     endif
     if(all(converged)) exit
     
     if(iter.eq.1) then
       ! Try a new value for the next iteration
       Fermi = Fermi + 0.1
       df = 0.1
     else
       ! Secant update
       df = -dn(:,1) * df(:)/(dn(:,1) - dn(:,2))
       
       ! Safeguard against large steps
       do it=1,2
         if(abs(df(it)) .gt. 1.0) then
            df(it) = df(it)/abs(df(it))
         endif
       enddo
       ! Update
       do it=1,2
          if(converged(it)) cycle            ! Do not iterate when close enough                         
          Fermi(it) = Fermi(it) + df(it) 
       enddo
     endif
     !--------------------------------------------------------------------------
  enddo

  !-----------------------------------------------------------------------------
  ! Linear mixing of the generalized density matrix, if asked for.
  if(.not.all(configmatrix_history.eq.0.0)) then
        if(HFBmixtype .eq. 1) then
            configmatrix=HFBmix*configmatrix+(1-HFBmix)*configmatrix_history
        endif   
  endif  
  !-----------------------------------------------------------------------------
  ! Construct the density and anomalous density matrix.
  ! 
  ! rho   =  U   f U^{\dagger} + V^{*}(1 - f)V^{T}
  ! kappa =  U   f V^{\dagger} + V^{*}(1 - f)U^{T} 
  rho_pairing = 0.0 ; kappa_pairing = 0.0

  si = 0 ; sb = 0
  do B=1,4
    N = HFBsizes(B)
    do i=1,N
      do j=1,N
        do k=1,N
          !---------------------------------------------------------------------
          !                                   U      f  U^{\dagger}
          rho_pairing(si+i,si+j)  = rho_pairing(si+i,si+j) +                   &
          &       configmatrix(sb+  k)*vect(sb+  i,sb+N+k) * vect(sb+  j,sb+N+k)
          !                                   V^* (1-f) V^{T}
          rho_pairing(si+i,si+j)  = rho_pairing(si+i,si+j) +                   &
          &       configmatrix(sb+N+k)*vect(sb+N+i,sb+N+k) * vect(sb+N+j,sb+N+k)

          !---------------------------------------------------------------------
          !                                   U   f        V^{\dagger}     
          ! Note the minus sign due to the hidden time-reversal!
          kappa_pairing(si+i,si+j)  = kappa_pairing(si+i,si+j) -               &
          &       configmatrix(sb+  k)*vect(sb+  i,sb+N+k) * vect(sb+N+j,sb+N+k)
          !                                   V^{*}(1 - f) U^{T} 
          kappa_pairing(si+i,si+j)  = kappa_pairing(si+i,si+j) +               &
          &       configmatrix(sb+N+k)*vect(sb+N+i,sb+N+k) * vect(sb  +j,sb+N+k)
        enddo
      enddo
    enddo

    si = si +   N
    sb = sb + 2*N
  enddo

  !-----------------------------------------------------------------------------
  ! Linear mixing of rho and kappa, if asked for.  
  if(.not.all(rho_history.eq.0.0)) then
        if(HFBmixtype .eq. 0) then
            rho_pairing   =  HFBmix * rho_pairing   + (1-HFBmix) * rho_history
            kappa_pairing =  HFBmix * kappa_pairing + (1-HFBmix) * kappa_history
        endif
  endif  

  !-----------------------------------------------------------------------------
  ! Side-effect, calculate the dispersion
  si            = 0
  HFBdispersion = 0.0
  do B=1,4
    N = HFBsizes(B) 
    it  = 1
    if(B > 2) it = 2
    
    ! rho squared
    chi(si+1:si+N, si+1:si+N) = matmul(rho_pairing(si+1:si+N, si+1:si+N),      &
    &                                        rho_pairing(si+1:si+N, si+1:si+N) )    

    !                                       Tr rho - Tr rho^2
    do i=1,N
        HFBdispersion(it) = HFBdispersion(it) + rho_pairing(si+i, si+i)        &
        &                                     - chi(si+i, si+i) 
    enddo
    si = si + N
  enddo
  ! Time-reversal
  HFBdispersion = 4 * HFBdispersion 
  !-----------------------------------------------------------------------------
 end subroutine solvepairing_HFB

 function ConstructConfiguration(Vect, QPenergies, BlockType,BlockIndices) &
                                      & result(R)
    !---------------------------------------------------------------------------
    ! Construct the configuration matrix C, to determine what kind of HFB  
    ! state we are aiming to construct.
    !---------------------------------------------------------------------------
    integer, intent(in)       :: BlockType
    integer, intent(in)       :: BlockIndices(:)
    real(KIND=dp), intent(in) :: qpenergies(:), vect(:,:)
    real(KIND=dp), allocatable:: R(:)
    
    integer :: N, B, sb, i, NB, j, block,  qblock, ind, si
    real(KIND=dp) :: compare, occ

    N = size(qpenergies) 
    allocate(R(2*N)) ;  R = 0
    
    !---------------------------------------------------------------------------
    ! Construct the DEFAULT configuration, corresponding to all positive energy
    ! quasiparticles.
    !---------------------------------------------------------------------------
    sb = 0
    si = 0    
    do B=1,4
        N = HFBsizes(B)   
        do i=1,N
            if(inversetemp .gt. 0.0_dp) then
                !---------------------------------------------------------------
                ! At finite temperature, things can get partially occupied and 
                ! we are dealing with a statistical mixture.
                occ = exp(inversetemp * Qpenergies(si+i))
                occ = 1.0/occ
                R(sb+N+i) = 1.0_dp - occ
                R(sb  +i) =          occ
            else
                !---------------------------------------------------------------
                !Completely empty or full, we want pure HFB states.
                R(sb+N+i) = 1.0_dp
                R(sb  +i) = 0.0_dp       
            endif
        enddo
        si = si +   N 
        sb = sb + 2*N        
    enddo

    !---------------------------------------------------------------------------
    ! Modify this default configuration when needed.
    select case(Blocktype)
    case(0)
        !-----------------------------------------------------------------------
        ! No blocking asked for. 
    case(1)
        !-----------------------------------------------------------------------
        ! Ordinary blocking.
        print *, 'Time-reversal breaking needed and not implemented.'
        stop
    case(2)
        !-----------------------------------------------------------------------
        ! EFA blocking
        NB = size(BlockIndices)

        do j=1,NB
            compare = 0.0
            ind     = 0

            ! Check which block we are dealing with
            call Identify(Blockindices(j), block, qblock)            
            !-------------------------------------------------------------------
            ! Look for the column in the second half of the eigenvectors with
            ! the largest overlap with asked for state.
            sb = 0
            do B=1,4
                N = HFBsizes(B)
                if(B.eq.qblock) then
                    do i=N+1,2*N
                        if(vect(sb+block,sb+i)**2 .gt. compare) then
                            compare = vect(sb+block+N,sb+i)**2 
                            ind     = i
                        endif
                    enddo
                endif
                sb = sb +2*N
            enddo 

!            print *, 'Blocking column', ind, ind - HFBsizes(qblock)  
            
            sb = 0 ; si =0 
            do B=1,4
                N = HFBsizes(B)
                if(qblock.eq.B) then
                    R(sb+ind-N)       = 0.5
                    R(sb+ind)         = 0.5
                endif
                sb = sb + 2*N
                si = si + N
            enddo    
        enddo
    end select

!    sb = 0
!    do B=1,4
!        N = HFBsizes(B)
!        do i=1,2*N
!            print *,' Config', i, R(sb+i)
!        enddo
!        print *
!        sb = sb + 2*N
!    enddo

 end function constructconfiguration
 
 subroutine Identify(i, bi, qblock)
    !---------------------------------------------------------------------------
    ! Identifies both the symmetry block(qblock) and index in said symmetry 
    ! block (bi), based on the index i in the HF basis. 
    !---------------------------------------------------------------------------
    integer, intent(in)  :: i
    integer, intent(out) :: bi, qblock    
        
    integer :: sb, N, B

    sb = 0
    do B=1,4
        N = HFBSizes(B)
        if( i .gt. sb .and. i.le.sb+N) then
            qblock = B
            bi     = i - sb
            return
        endif
        sb = sb + N
    enddo
 end subroutine Identify

 function ConstructHFBHamil(sphamil, gaps, Fermi) result(H)
    !---------------------------------------------------------------------------
    ! Construct the HFB hamiltonian
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in) :: sphamil(:,:), gaps(:,:), Fermi
    real(KIND=dp), allocatable  :: H(:,:)
    
    integer :: N, i
    
    N = size(sphamil,1)
    allocate(H(2*N,2*N)) ; H = 0
    
    H(1:N, 1:N)         =  sphamil
    H(N+1:2*N, N+1:2*N) = -sphamil
    
    do i=1,N
      H(i,i)      = H(i,i)      - Fermi
      H(i+N, i+N) = H(i+N, i+N) + Fermi
    enddo
    
    H(1:N, N+1:2*N)     =  gaps
    H(N+1:2*N, 1:N)     = -gaps
    
 end function ConstructHFBHamil

 subroutine Canonical(rho_pairing, kappa_pairing, rho_can, kappa_can,          &
 &                                                     rhotransfo, kappatransfo)
  !-----------------------------------------------------------------------------
  ! a) Diagonalize  Rho
  ! b) Canonicalize Kappa
  ! c) Return the diagonal elements of rho, and the offdiagonal elements of 
  !    kappa, as well as the transformation.
  !-----------------------------------------------------------------------------
  
  real(KIND=dp), intent(in)  :: rho_pairing(nwt,nwt)
  real(KIND=dp), intent(in)  :: kappa_pairing(nwt,nwt)
  real(KIND=dp), intent(out) :: rho_can(nwt), kappa_can(nwt)
  real(KIND=dp), intent(out) :: rhotransfo(nwt,nwt), kappatransfo(nwt,nwt)
  real(KIND=dp)              :: temptransfo(nwt,nwt)
  
  real(KIND=dp) :: work(2*nwt), E, X(nwt), Y(nwt)
  real(KIND=dp), allocatable :: tmp(:,:), cpy(:)
  
  integer :: si,N, B, i, sb, j
  
  !-----------------------------------------------------------------------------
  ! a) Diagonalize rho
  !
  ! Transforms as rho'  = D^dagger rho D
  !-----------------------------------------------------------------------------
  si         = 0
  rhotransfo = 0 ; kappatransfo = 0
  do B=1,4
    N =HFBsizes(B)
    
    allocate(tmp(N,N))
    
    tmp = rho_pairing(si+1:si+N, si+1:si+N)
    
    ! Diagonalize rho in this block
    call diagon(tmp,N,N,rhotransfo(si+1:si+N,si+1:si+N),rho_can(si+1:si+N),work)

!    ! DEBUG
!    do i=1,N
!       print ('(20f10.3)'), rho_pairing(si+i,si+1:si+N)        
!    enddo
!    print *
!    print ('(a6, i3, 20f10.3)'), 'Occ B=',B, rho_can(si+1:si+N)
!    print *

    si = si + N
    deallocate(tmp)
  enddo
  
  ! Time-reversal
  rho_can = 2*rho_can
  
  do i=1,nwt
    if(rho_can(i).gt.2.0) rho_can(i) = 2.0
    if(rho_can(i).lt.0.0) rho_can(i) = 0.0
  enddo
  
  !-----------------------------------------------------------------------------
  ! b) Bring kappa (with cutoffs) into canonical form
  !
  ! Transforms as kappa'  = D^dagger kappa D^*
  !-----------------------------------------------------------------------------
  si = 0
  do B=1,4
  
    N = HFBsizes(B)
    
    allocate(tmp(N,N))
    tmp = kappa_pairing(si+1:si+N, si+1:si+N)   

    tmp = matmul(transpose(rhotransfo(si+1:si+N, si+1:si+N)), tmp)
    tmp = matmul(tmp,rhotransfo(si+1:si+N, si+1:si+N))
  
    ! With the assumption of time-reversal, the diagonal matrix elements in this
    ! transformed kappa matrix are the matrix elements (i, ibar).
    !
    ! They are used in the BCS case only:
    do i=1,N
        kappa_can(si+i) = tmp(i,i)
    enddo
    
    deallocate(tmp)
    si = si + N
  enddo    
  

   !-----------------------------------------------------------------------------
!  ! b) Bring kappa (with cutoffs) into canonical form
!  !
!  ! Transforms as kappa'  = D^dagger kappa D^*
!  !-----------------------------------------------------------------------------
!  si = 0
!  do B=1,4
!  
!    N = HFBsizes(B)
!    
!    allocate(tmp(N,N))
!    
!!    print *, 'Kappa'
!!    do i=1,N
!!      print ('(99f8.2)'), kappa_pairing(si+i, si+1:si+N)
!!    enddo
!!    print *

!    !---------------------------------------------------------------------------
!    ! Construct kappa with cutoffs in this block
!    tmp = kappa_pairing(si+1:si+N, si+1:si+N) 
!    do i=1,N
!      do j=1,N
!        tmp(i,j) = tmp(i,j) * Pcutoffs(si+i) * Pcutoffs(si+j)
!      enddo
!    enddo
!    
!    ! Lets diagonalize -(kappa)^2
!    tmp = - matmul(tmp, tmp)
!    call diagon(tmp,N,N,temptransfo(si+1:si+N,si+1:si+N),kappa_can(si+1:si+N), & 
!    &                                                                      work)
!    
!    print *, 'kappa 2'
!    print *, kappa_can(si+1:si+N)
!    
!!    do i=1,N
!!        print ('(90f8.2)') , temptransfo(si+1:si+N,i)
!!    enddo
!!    stop
!    
!    ! Reconstruct kappa
!    tmp = kappa_pairing(si+1:si+N, si+1:si+N) 
!    do i=1,N
!      do j=1,N
!        tmp(i,j) = tmp(i,j) * Pcutoffs(si+i) * Pcutoffs(si+j)
!      enddo
!    enddo
!    
!    do i=1,N,2
!      ! Every pair of eigenvalues K^2 
!      E = sqrt(kappa_can(si+i))
!      kappa_can(si+i)   =  E
!      kappa_can(si+i+1) = -E
!    
!      X(1:N) = temptransfo(si+1:si+N, si+i) 
!    
!      kappatransfo(si+1:si+N, si+i)   =  X(1:N)
!      Y(1:N) =                             matmul(tmp,X(1:N))
!      kappatransfo(si+1:si+N, si+i+1) = Y(1:N)/sum(Y(1:N)**2)       
!    enddo
!    
!    do i=1,N
!        print ('(90f8.2)') , kappatransfo(si+1:si+N,i)
!    enddo
!    stop
!    
!    stop
!    !---------------------------------------------------------------------------

!    si = si + N
!    deallocate(tmp)
!  enddo
!  print *, 'Canonical'
!  stop
 end subroutine Canonical
 
 subroutine ConstructCanonicalBasis(Transfo)
  !-----------------------------------------------------------------------------
  ! Transform the spwf wavefunctions in the HFBasis into Canbasis, with the 
  ! passed in Transfo. 
  !-----------------------------------------------------------------------------

  integer                   :: wave1, wave2, B, N, si
  real(KIND=dp), intent(in) :: Transfo(nwt,nwt)
  
  if(.not.allocated(CanPsi)) then
    allocate(CanPsi(mv,4,nwt)) ; CanPsi      = 0.0
    allocate(Canenergies(nwt)) ; canenergies = 0.0
  endif

  si     = 0   
  CanPsi = 0.0 ; canenergies = 0.0
  do B=1,8
    N = HFBlocks(B)
    !---------------------------------------------------------------------------
    ! Apply the transformation in this symmetry block
    do wave1=1, N 
      do wave2=1,N
        CanPsi(:,:,si+wave1)  = CanPsi(:,:,si+wave1) +                         &
        &                       Transfo(si+wave2,si+wave1) * HFPsi(:,:,si+wave2) 
        
        canenergies(si+wave1) = canenergies(si+wave1) +                        &
        &               abs(Transfo(si+wave2,si+wave1)**2) *spenergies(si+wave2) 
      enddo 
    enddo
    
    si = si +  N
    !---------------------------------------------------------------------------
  enddo
  
 end subroutine ConstructCanonicalBasis

 subroutine calcHFBgaps(Fermi)
  !-----------------------------------------------------------------------------
  ! Calculates the HFB gaps for use in the HFB solver.
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(in) :: Fermi(2)
  integer                   :: wave1, wave2, iso, si, B, N
  real(KIND=dp)             :: deltapsi(mv,4)
  
  if(.not.associated(Delta_action_HFB)) stop
  !-----------------------------------------------------------------------------
  ! Use the delta_action to calculate the elements in the gaps
  si = 0
  HFBgaps = 0
  do B=1,8
    N = HFBlocks(B)
    
    iso = -1
    if(B>4) iso = 1
  
    do wave1=1,N
      deltapsi = delta_action_HFB(  hfpsi(:,:,si+wave1)  ,                   &
      &                            hfdpsi(:,:,:,si+wave1),                   &
      &                           hfddpsi(:,:,:,si+wave1),                   &
      &                          hfdddpsi(:,:,:,si+wave1),                   &
      &              sx(:,si+wave1), sy(:,si+wave1), sz(:,si+wave1),iso,.false.)
      
      ! Mystery factor 0.5 in here
      do wave2=wave1,N
        HFBgaps(si+wave2,si+wave1) = 0.5*sum(hfpsi(:,:,si+wave2)*deltapsi)*dv *&
        &                                  Pcutoffs(si+wave1)*Pcutoffs(si+wave2)
        HFBgaps(si+wave1,si+wave2) = - HFBgaps(si+wave2,si+wave1)
      enddo
    enddo
    si = si + N
  enddo

 end subroutine calcHFBgaps
 
!===============================================================================
!  Never to be used function to define an interface for delta_action
!===============================================================================   
 function delta_action_dummy(psi, dpsi, ddpsi, dddpsi, sx,sy,sz,iso, onthefly) &
                                                                result(deltapsi)
    !---------------------------------------------------------------------------
    ! Dummy function to allow this module to acces the functional.f90 module 
    ! to acces the information on the acces of deltas.
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)    :: psi(mv,4)  
    real(KIND=dp), intent(inout) :: dpsi(mv,3,4),ddpsi(mv,6,4), dddpsi(mv,10,4)
    integer, intent(in)       :: sx(4),sy(4),sz(4),   iso
    real(KIND=dp)             :: deltapsi(mv,4)
    real(KIND=dp)             :: temp(mv,4)
    real(KIND=dp)             ::   dtemp(mv,3,4)
    real(KIND=dp)             ::  ddtemp(mv,3,3,4)
    real(KIND=dp)             :: dddtemp(mv,3,3,3,4)
    real(KIND=dp)             :: laptemp(mv,4)
    logical, intent(in)       :: onthefly
 end function

end module HFB


!===============================================================================
! CODE ZOO
!===============================================================================
!  si = 0
!  do B=1,4
!  
!    N = HFBsizes(B)
!    
!    allocate(tmp(N,N))
!    
!    print *, 'Kappa'
!    do i=1,N
!      print ('(99f8.2)'), kappa_pairing(si+i, si+1:si+N)
!    enddo
!    print *

!    !---------------------------------------------------------------------------
!    ! We transform kappa**2, this is slightly easier to manipulate
!    tmp = kappa_pairing(si+1:si+N, si+1:si+N) !matmul(kappa_pairing(si+1:si+N, si+1:si+N), &
!    !&                                        kappa_pairing(si+1:si+N,si+1:si+N))
!    
!    ! Apply the transformation
!    tmp = matmul(tmp,transfo(si+1:si+N,si+1:si+N))
!    tmp = matmul(transpose(transfo(si+1:si+N,si+1:si+N)), tmp)

!    print *
!    print *, 'Transformed kappa'
!    do i=1,N
!      print ('(99f8.2)'), tmp(i, 1:N)
!    enddo
!    print *

!!    do i=1,N
!!      kappa_can(si + i) = sqrt(abs(tmp(i,i)))
!!    enddo
!!    print *, 'Canonical elements of kappa'
!!    print *, kappa_can(si+1:si+N) 
!!    print *

!    ! Construct kappa with cutoffs
!    tmp = kappa_pairing(si+1:si+N, si+1:si+N)
!    do i=1,N
!      do j=1,N
!        tmp(i,j) = tmp(i,j) * Pcutoffs(si+i) * Pcutoffs(si+j)
!      enddo
!    enddo
!    
!    print *
!    print *, 'Kappa with cutoffs'
!    do i=1,N
!      print ('(99f8.2)'), tmp(i, 1:N)
!    enddo
!    print *
!    !---------------------------------------------------------------------------
!    ! Apply the transformation
!    tmp = matmul(tmp,transfo(si+1:si+N,si+1:si+N))
!    tmp = matmul(transpose(transfo(si+1:si+N,si+1:si+N)), tmp)

!    print *
!    print *, 'Transformed kappa with cutoffs'
!    do i=1,N
!      print ('(99f8.2)'), tmp(i, 1:N)
!    enddo
!    print *
!    print *, '-----------------------------------'


!    si = si + N
!    deallocate(tmp)
!  enddo
!  print *, 'Canonical'
!  stop

!
