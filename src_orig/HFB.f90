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
 ! Module that is useful for solving the HFB problem. 
 !
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
  integer :: maxHFBiter  = 100
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

 subroutine solvepairing_HFB(fermi, rho_pairing, kappa_pairing, qpenergies)
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
  
  real(KIND=dp)              :: fermi(2)
  real(KIND=dp)              :: rho_pairing(nwt,nwt), qpenergies(nwt)
  real(KIND=dp)              :: kappa_pairing(nwt,nwt)
  real(KIND=dp)              :: df(2), dn(2,2)

  real(KIND=dp)              :: sphamil(nwt,nwt),vect(2*nwt,2*nwt)
  real(KIND=dp)              :: eigen(2*nwt),  work(nwt)
  real(KIND=dp)              :: particles(2), temp(nwt), chi(nwt,nwt)
  real(KIND=dp), allocatable :: HFBHamil(:,:)
  
  logical                    :: converged(2)
  
  integer                    :: si, sb, N, B, iter, wave1, it, i,j,k
  
  dn = 0.0
  
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

        !-----------------------------------------------------------------------
        ! Switching half of the eigenvectors
        eigen(sb+1:sb+N) = -eigen(sb+1:sb+N)
        do i = 1,N
          temp                      = vect(sb+1  :sb+N,   sb+i)
          vect(sb  +1:sb+N,   sb+i) = vect(sb+N+1:sb+2*N, sb+i)
          vect(sb+N+1:sb+2*N, sb+i) = temp 
        enddo

        !-----------------------------------------------------------------------
        ! Saving the quasiparticle excitation energies
        qpenergies(si+1:si+N) = eigen(sb+N+1:sb+2*N)

        !-----------------------------------------------------------------------
        ! Calculate the number of particles in here  
        ! Sum_i rho_ii = sum_ij V^*_ij V^T_ji = sum_ij V^*_ij V_ij
        particles(it) =  particles(it)  & 
        &                   + 2*sum(vect(sb+N+1:sb+2*N, sb+N+1:sb+2*N)**2)
        ! Time reversal is responsible for the factor 2

        deallocate(HFBHamil)
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
  ! Construct the density and anomalous density matrix.
  rho_pairing = 0.0 ; kappa_pairing = 0.0

  si = 0 ; sb = 0
  particles = 0
  do B=1,4
    N = HFBsizes(B)
    do i=1,N
      do j=1,N
        do k=1,N
          rho_pairing(si+i,si+j)  = rho_pairing(si+i,si+j) +                   &
          &                            vect(sb+N+i,sb+N+k) * vect(sb+N+j,sb+N+k)
          kappa_pairing(si+i,si+j)  = kappa_pairing(si+i,si+j) +               &
          &                            vect(sb+N+i,sb+N+k) * vect(sb  +j,sb+N+k) 
        enddo
      enddo
    enddo
 
    si = si +   N
    sb = sb + 2*N
  enddo
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
  enddo
  !-----------------------------------------------------------------------------
 end subroutine solvepairing_HFB
 
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

    si = si + N
    deallocate(tmp)
  enddo
  
  ! Time-reversal
  rho_can = 2*rho_can
  
  do i=1,nwt
    if(rho_can(i).gt.2.0) rho_can(i) = 2.0
    if(rho_can(i).lt.0.0) rho_can(i) = 0.0
  enddo
  
!  !-----------------------------------------------------------------------------
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
