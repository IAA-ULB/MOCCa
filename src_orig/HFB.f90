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
  
  implicit none
  
  !-----------------------------------------------------------------------------
  ! Gaps and quasiparticle energies
  real(KIND=dp), allocatable :: HFBGaps(:,:)
  real(KIND=dp), allocatable :: HFBqps(:)
  
  ! Maximum amount of iterations for finding a Fermi energy
  integer :: maxHFBiter = 200
  integer :: HFBsizes(4) = 0
  
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
  HFBSizes(1) = 2*HFblocks(1)
  HFBSizes(2) = 2*HFblocks(3)
  HFBSizes(3) = 2*HFblocks(5)
  HFBSizes(4) = 2*HFblocks(7)
  
  ! Otherwise the sizes of the blocks are simply the same as the HFBlocks
  
 end subroutine initHFB

 subroutine solvepairing_HFB(fermi, rho_pairing, kappa_pairing)
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
  real(KIND=dp)              :: rho_pairing(2*nwt,2*nwt)
  real(KIND=dp)              :: kappa_pairing(2*nwt,2*nwt)
  real(KIND=dp)              :: df(2), dn(2,2)

  real(KIND=dp)              :: sphamil(2*nwt,2*nwt),  vect(4*nwt,4*nwt)
  real(KIND=dp)              :: eigen(4*nwt),  work(2*nwt)
  real(KIND=dp)              :: particles(2)
  real(KIND=dp), allocatable :: HFBHamil(:,:)
  
  logical                    :: converged(2)
  
  integer                    :: si, sb, N, B, iter, wave1, it, i,j,k
  
  dn = 0.0
  
  Fermi = -10
  
  !-----------------------------------------------------------------------------
  ! Construct the single-particle hamiltonian from the sp.energes
  sphamil = 0
  si      = 0
  sb      = 0
  do B=1,8,2
    N = HFBlocks(B)
    do wave1=1,N
      sphamil(sb+wave1  ,sb+wave1  ) = spenergies(si+wave1)
      sphamil(sb+wave1+N,sb+wave1+N) = spenergies(si+wave1)  ! Time-reversed 
    enddo
    si = si +   N
    sb = sb + 2*N
  enddo
  
  df = 0.0
  converged = .false.
  
  !-----------------------------------------------------------------------------
  ! Start iterations over the Fermi energy
  do iter=1,maxHFBiter
      
     vect      = 0.0
     eigen     = 0.0
     particles = 0.0
     !--------------------------------------------------------------------------
     ! Diagonalize every block. 
     si = 0
     sb = 0
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
        ! Calculate the number of particles in here  
        ! Sum_i rho_ii = sum_ij V^*_ij V^T_ji = sum_ij V^*_ij V_ij
        particles(it) =  particles(it)  & 
        &                     + sum(vect(sb+N+1:sb+2*N, sb+N+1:sb+2*N)**2)

        deallocate(HFBHamil)
        ! Startindex (si) for the next block.
        si = si + N
        sb = sb + 2*N

     enddo
     !--------------------------------------------------------------------------
     ! Adjust Fermi energy, using a secant method for the moment
     dn(:,2) = dn(:,1)
     
     dn(1,1) = particles(1) - neutrons
     dn(2,1) = particles(2) - protons
     
     print *, iter, particles

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
         if(abs(df(it)) .gt. 1) then
            df(it) = df(it)/abs(df(it))
         endif
       enddo
       ! Update
       do it=1,2
        if(converged(it)) cycle            ! Do not iterate when close enough                         
        Fermi(it) = Fermi(it) + df(it) 
       enddo
     endif
     
  enddo
  !-----------------------------------------------------------------------------
  ! Construct the density and anomalous density matrix.
  rho_pairing = 0.0 ; kappa_pairing = 0.0

  si = 0
  sb = 0
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
    
    do i=1,N
        print ('(99f8.2)'), kappa_pairing(si+i, si+1:si+N)
    enddo
    print *

    si = si +  N
    sb = sb +2*N
  enddo
  
  do i=1,2*nwt
        print ('(99f8.2)'), kappa_pairing(i, 1:2*nwt)
  enddo
  print *

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

 subroutine Canonical(rho_pairing, kappa_pairing, rho_can, kappa_can, transfo)
  !-----------------------------------------------------------------------------
  ! a) Diagonalize  Rho
  ! b) Canonicalize Kappa
  ! c) Return the diagonal elements of rho, and the offdiagonal elements of 
  !    kappa, as well as the transformation.
  !-----------------------------------------------------------------------------
  
  real(KIND=dp), intent(in)  :: rho_pairing(2*nwt,2*nwt)
  real(KIND=dp), intent(in)  :: kappa_pairing(2*nwt,2*nwt)
  real(KIND=dp), intent(out) :: rho_can(2*nwt), kappa_can(2*nwt)
  real(KIND=dp), intent(out) :: transfo(2*nwt,2*nwt)
  
  real(KIND=dp) :: work(2*nwt)
  real(KIND=dp), allocatable :: tmp(:,:)
  
  integer :: si,N, B, i
  
  !-----------------------------------------------------------------------------
  ! a) Diagonalize rho
  !
  ! Transforms as rho'  = D^dagger rho D
  !-----------------------------------------------------------------------------
  si      = 0
  transfo = 0
  do B=1,4
    N =HFBsizes(B)
    
    allocate(tmp(N,N))
    
    tmp = rho_pairing(si+1:si+N, si+1:si+N)


    ! Diagonalize first part
    call diagon(tmp(1:N/2, 1:N/2),N/2,N/2,transfo(si+1:si+N/2,si+1:si+N/2),    &
    &                                                 rho_can(si+1:si+N/2),work)
    
    ! Diagonalize second part
    ! (We could have put this because of time-reversal symmetry)
    call diagon(tmp(N/2+1:N, N/2+1:N),N/2,N/2,                             &
    &        transfo(si+N/2+1:si+N,si+N/2+1:si+N),rho_can(si+N/2+1:si+N),  work)
    
    si = si + N
    deallocate(tmp)
  enddo
  !-----------------------------------------------------------------------------
  ! b) Bring kappa into canonical form
  !
  ! This is currently applying the cantransfo deduced from the diagonalisation
  ! of rho for a mean-field calculation. 
  !
  ! When no antilinear, antihermitian symmetry is conserved, we need to take
  ! out an additional phase here!
  !
  ! Transforms as kappa'  = D^dagger rho D^*
  !-----------------------------------------------------------------------------
  si = 0
  do B=1,4
  
    N = HFBsizes(B)
    
    allocate(tmp(N,N))
    
    !---------------------------------------------------------------------------
    ! We transform kappa**2, this is slightly easier to manipulate
    tmp = matmul(kappa_pairing(si+1:si+N, si+1:si+N), &
    &                                        kappa_pairing(si+1:si+N,si+1:si+N))
    ! Apply the transformation
    tmp = matmul(tmp,transfo(si+1:si+N,si+1:si+N))
    tmp = matmul(transpose(transfo(si+1:si+N,si+1:si+N)), tmp)

    do i=1,N
      kappa_can(si + i) = sqrt(abs(tmp(i,i)))
    enddo
    si = si + N
    deallocate(tmp)
  enddo
  
 end subroutine Canonical
 
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
