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

  use diag
  use geninfo
  use wavefunctions
  use pairingcutoffs

  implicit none

  !-----------------------------------------------------------------------------
  ! Matrix of all qp energies
  real(KIND=dp), allocatable :: HFBGaps(:,:)
  ! Maximum amount of iterations for finding a Fermi energy
  integer :: maxHFBiter  = 200
  ! Dimensions of all the blocks in the HFBhamiltonian.
  integer, allocatable :: HFBsizes(:)
  ! Dispersion of the particle number
  real(KIND=dp) :: HFBdispersion(2)
  ! Pointer to relink procedures
  procedure(delta_action_dummy), pointer :: delta_action_HFB
  !-----------------------------------------------------------------------------

  contains

  subroutine initHFB
    !---------------------------------------------------------------------------
    ! Determine the matrix sizes of the HFB problem. 
    !---------------------------------------------------------------------------

    integer :: i

    allocate(HFBsizes(8))

    do i=1,8
      HFBSizes(i) = HFblocks(i)
    enddo

  end subroutine initHFB

  subroutine solvepairing_HFB(fermi, Bogoliubov, rho_pairing, kappa_pairing,   &
  &                           configmatrix, qpenergies,HFBmix, HFBmixtype,     &
  &                                         BlockType,Blockindices, blocklowest)

    !---------------------------------------------------------------------------
    ! Driver routine for the solving of the HFB equations, represented in the 
    ! Hartree-Fock basis.
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(inout) :: Fermi(2)
    real(KIND=dp), intent(inout) :: Bogoliubov(:,:)
    real(KIND=dp), intent(inout) :: kappa_pairing(:,:), rho_pairing(:,:)
    real(KIND=dp), intent(inout) :: configmatrix(:), qpenergies(:)
    
    !---------------------------------------------------------------------------
    ! History of the pairing matrices, for mixing purposes.
    real(KIND=dp), allocatable ::  rho_history(:,:), kappa_history(:,:)
    real(KIND=dp), allocatable ::  configmatrix_history(:) 
    real(KIND=dp), allocatable ::  Bogoliubov_history(:,:)

    ! Options for the mixing of the HFB configurations
    real(KIND=dp), intent(in)    :: HFBmix
    integer, intent(in)          :: HFBmixtype

    ! Configuration for the blocking
    integer, intent(in)          :: Blockindices(:)
    integer, intent(in)          :: BlockType
    character(len=2), intent(in) :: BlockLowest(:)
    integer, allocatable         :: neutron_block(:), proton_block(:)

    ! Quantities for the HFB hamiltonian
    real(KIND=dp)              :: sphamil(nwt,nwt), HFBHamil(2*nwt, 2*nwt)
    
    integer                    :: si, sb, N, B, iter, wave1, it, i,j,k, np, nn
    integer                     :: n_ind, p_ind
    real(KIND=dp), allocatable :: chi(:,:)
    !-----------------END OF DECLARATIONS --------------------------------------

    if(.not.allocated(rho_history)) then
          allocate(rho_history(nwt,nwt))            ; rho_history   = 0.0
          allocate(kappa_history(nwt,nwt))          ; kappa_history = 0.0
          allocate(configmatrix_history(2*nwt))     ; configmatrix_history = 0.0
          allocate(Bogoliubov_history(2*nwt, 2*nwt)); Bogoliubov_history = 0.0
    endif

    ! Saving the history
    rho_history          = rho_pairing
    kappa_history        = kappa_pairing
    configmatrix_history = configmatrix
    Bogoliubov_history   = Bogoliubov

    ! Guess a new Fermi energy if none is there
    if(all(Fermi.eq.0.0))   Fermi = -5

    ! Preparation for blocking, we need to separate configurations by isospin.
    select case (Blocktype)
    case(0)
      allocate(proton_block(1))
      allocate(neutron_block(1))
    case(1,3)
      ! We pass to the configuration routine the indices of all quasiparticles
      ! to excite.
      nn = 0 ; np = 0
      do i=1,size(blockindices)
        if(blockindices(i).le.nwn) then
          nn = nn + 1 
        else
          np = np + 1
        endif
      enddo
        
      allocate(proton_block(np))
      allocate(neutron_block(nn))

      n_ind = 1 ; p_ind = 1
      do i=1,size(blockindices)
        if(blockindices(i).le.nwn) then
          neutron_block(n_ind) = blockindices(i)
          n_ind = n_ind + 1
        else
          proton_block(p_ind)  = blockindices(i) - nwn
          p_ind = p_ind + 1
        endif
      enddo

    case(2,4) 
      ! We pass to the configuration routine the number of quasiparticles to
      ! excite in every symmetry block. (The fifth entry tells the code to 
      ! construct the overall lowest.)
      allocate(proton_block(5))  ; proton_block  = 0 
      allocate(neutron_block(5)) ; neutron_block = 0
         
      do i=1, size(Blocklowest)
        select case (Blocklowest(i))
        case('n+')
            neutron_block(1) = neutron_block(1) + 1 
        case('n-')
            neutron_block(3) = neutron_block(3) + 1 
        case('p+')
            proton_block(1)  = proton_block(1)  + 1
        case('p-')
            proton_block(3)  = proton_block(3)  + 1
        case('n0')
            neutron_block(5) = neutron_block(5) + 1
        case('p0')
            proton_block(5)  = proton_block(5)  + 1
        end select
      enddo
    end select

    !---------------------------------------------------------------------------
    ! a) We construct the HFB-hamiltonian in every block. 
    sphamil = 0
    si      = 0 ; sb = 0
    do B=1,8
      N = HFBlocks(B) ;  if(N .eq. 0) cycle 

      do wave1=1,N
        sphamil(si+wave1,si+wave1) = spenergies(si+wave1)
      enddo
      HFBHamil(sb+1:sb+2*N, sb+1:sb+2*N) = ConstructHFBHamil(                  & 
      &               sphamil(si+1:si+N,si+1:si+N),HFBgaps(si+1:si+N,si+1:si+N))  
      
      si = si +   N
      sb = sb + 2*N
    enddo
    !---------------------------------------------------------------------------
    ! b) We repeatedly diagonalize the matrix to find a suitable Fermi energy, 
    !    for every isospin. 
    call FindFermi(HFBHamil(      1:2*nwn,      1:2*nwn), HFBsizes(1:4),       &
    &              neutrons, configmatrix(  1:2*nwn),                          &
    &              Bogoliubov(1:2*nwn, 1:2*nwn),                               &
    &              qpenergies(1:nwn),    Fermi(1), maxhfbiter,                 &
    &              blocktype, neutron_block)   

    call FindFermi(HFBHamil(2*nwn+1:2*nwt,2*nwn+1:2*nwt), HFBsizes(5:8),       &
    &              protons, configmatrix(2*nwn+1:2*nwt),                       &
    &              Bogoliubov(2*nwn+1:2*nwt, 2*nwn+1:2*nwt),                   &
    &              qpenergies(nwn+1:nwt),Fermi(2), maxhfbiter,                 &
    &              blocktype, proton_block)     

    !---------------------------------------------------------------------------
    ! c) Optionally mix the configuration matrices.  
    if(.not.all(configmatrix_history.eq.0.0)) then
      if(HFBmixtype .eq. 1) then
        configmatrix=HFBmix*configmatrix+(1-HFBmix)*configmatrix_history
      endif   
    endif  

    !---------------------------------------------------------------------------
    ! d) Construct the density and anomalous density matrices, based on the 
    !    configmatrix and the Bogoliubov transformation
    call PairingMatrices(configmatrix, bogoliubov, rho_pairing, kappa_pairing)
    
    !---------------------------------------------------------------------------
    ! e) Optionally mix these densities
    if(.not.all(rho_history.eq.0.0)) then
      if(HFBmixtype .eq. 0) then
        rho_pairing   =  HFBmix * rho_pairing   + (1-HFBmix) * rho_history
        kappa_pairing =  HFBmix * kappa_pairing + (1-HFBmix) * kappa_history
      endif
    endif  
    
    !---------------------------------------------------------------------------
    ! f) Calculate the dispersion, while we are at it. 
    si            = 0
    HFBdispersion = 0.0
    allocate(chi(nwt, nwt))
    do B=1,8
      N = HFBsizes(B)  ;  if(N .eq. 0) cycle 

      it  = 1
      if(B .gt. 4) it = 2
      
      ! rho squared
      chi(si+1:si+N, si+1:si+N) = matmul(rho_pairing(si+1:si+N, si+1:si+N),    &
      &                                      rho_pairing(si+1:si+N, si+1:si+N) )    

      !                                       Tr rho - Tr rho^2
      do i=1,N
          HFBdispersion(it) = HFBdispersion(it) + rho_pairing(si+i, si+i)      &
          &                                     - chi(si+i, si+i) 
      enddo
      si = si + N
    enddo
    deallocate(chi)
    ! Factor 2 from the formula
    HFBdispersion = 2 * HFBdispersion 
    ! Time-reversal introduces a factor of two
$TR    HFBdispersion = 2 * HFBdispersion 

  end subroutine solvepairing_HFB

  function ConstructConfiguration(Bogo, Eqp, blocks, blocktype, blockconf)     &
  &                              result(R)
    !---------------------------------------------------------------------------
    ! Construct the configuration matrix, based on the various user options.
    !
    !---------------------------------------------------------------------------

    integer, intent(in)          :: BlockType
    integer, intent(in)          :: blocks(4) 
    integer, intent(in)          :: blockconf(:)

    real(KIND=dp), intent(in)    :: Eqp(:), Bogo(:,:)
    real(KIND=dp), allocatable   :: R(:)
      
    integer                      :: N, B, sb, i, NB, j, qblock, ind, si, bi
    real(KIND=dp)                :: compare, occ, qpmin, qpb
    integer                      :: toblock(4)

    N = size(Eqp) 
    allocate(R(2*N)) ;  R = 0

    !---------------------------------------------------------------------------
    ! Construct the DEFAULT configuration, corresponding to all positive energy
    ! quasiparticles.
    sb = 0 ; si = 0      

    do B=1,4
        N = blocks(B) ; if (N.eq. 0) cycle
        do i=1,N
            if(inversetemp .gt. 0.0_dp) then
                !---------------------------------------------------------------
                ! At finite temperature, things can get partially occupied and 
                ! we are dealing with a statistical mixture.
                occ = exp(inversetemp * Eqp(si+i))
                occ = 1.0/(1 + occ)
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
    occ = 0
    select case(Blocktype)
    case(1,2)
        ! Full blocking
        occ = 1.0_dp
        print *, 'Time-reversal breaking is needed and not implemented.'
        stop
    case(3,4)
        ! EFA blocking
        occ = 0.5_dp
    end select
    !---------------------------------------------------------------------------
    ! Modify this default configuration when needed.
    select case(Blocktype)
    case(0)
        !-----------------------------------------------------------------------
        ! No blocking asked for. 
    case(1,3)
        !-----------------------------------------------------------------------
        ! The user asked for a specific configuration that needs to be 
        ! identified. The array blockconf now contains the indices in the 
        ! HF-basis.

        NB = size(blockconf)

        do j=1,NB
            compare = 0.0
            ind     = 0

            ! Check which block the requested index is in.
            call Identify(blockconf(j),blocks, bi, qblock)
            !-------------------------------------------------------------------
            ! Look for the column in the second half of the eigenvectors with
            ! the largest overlap with asked for state.
            sb = 0
            do B=1,4
                N = blocks(B) ; if (N.eq.0) cycle
                if(B.eq.qblock) then
                    do i=N+1,2*N
                        if(Bogo(sb+bi,sb+i)**2 .gt. compare) then
                            compare = Bogo(sb+bi+N,sb+i)**2 
                            ind     = i
                        endif
                    enddo
                endif
                sb = sb +2*N
            enddo 
            !-------------------------------------------------------------------
            !  Change the occupation of this particular qp
            sb = 0 ; si =0 
            do B=1,4
                N = HFBsizes(B)
                if(qblock.eq.B) then
                    R(sb+ind-N)       = occ 
                    R(sb+ind)         = 1 - occ
                endif
                sb = sb + 2*N
                si = si + N
            enddo    

        enddo    
    case(2,4)
        !-----------------------------------------------------------------------
        ! The user asked for a the lowest configuration of a specific type.
        ! In this case, blockconf contains the number of qp excitations to  
        ! construct in every block.

        toblock = blockconf(1:4)

        if(blockconf(5).ne.0) then
          do i = 1, blockconf(5)
            qpmin = 10000000
            si = 0
            do B=1,4
              N = blocks(B) ; if (N.eq.0) cycle

              if(Eqp(si+toblock(B)+1) .lt. qpmin) then
                qpmin = Eqp(si+toblock(B)+1)
                qpb   = B
              endif
              si = si +   N
            enddo
            toblock(qpb) = toblock(qpb) + 1
          enddo
        endif

        !  For every block, we flip the required number of qps.
        sb = 0 ; si = 0
        do B=1,4
          N = blocks(B) ; if (N.eq.0) cycle
          do j=1,toblock(B)
            ! The qps are ordered in energy from the diagonalization
            !  So we simply flip the first ones
            R(sb + N + j ) = 1 - occ
            R(sb     + j ) =     occ
          enddo
          si = si +   N
          sb = sb + 2*N
        enddo
    end select
   
  end function ConstructConfiguration

  subroutine PairingMatrices(config, bogo, rho, kappa)
    !---------------------------------------------------------------------------
    ! Calculate the matrices rho and kappa for a given configuration matrix 
    ! and Bogoliubov transformation.
    !
    !    rho   = U   f U^\dagger + V^* (1 - f) V^T
    !    kappa = U   f V^\dagger + V^* (1 - f) U^T  
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: config(:), Bogo(:,:)
    real(KIND=dp), intent(out):: rho(:,:), kappa(:,:)
    integer                   :: si, sb, B, i,j,k, N

    rho = 0.0 ; kappa = 0.0

    si = 0 ; sb = 0
    do B=1,8
      N = HFBsizes(B) ;  if(N .eq. 0) cycle 
      do i=1,N
        do j=1,N
          do k=1,N
            !-------------------------------------------------------------------
            !                                   U      f  U^{\dagger}
            rho(si+i,si+j)  = rho(si+i,si+j) +                 &
            &           config(sb+  k)*bogo(sb+  i,sb+N+k) * bogo(sb+  j,sb+N+k)
            !                                   V^* (1-f) V^{T}
            rho(si+i,si+j)  = rho(si+i,si+j) +                 &
            &           config(sb+N+k)*bogo(sb+N+i,sb+N+k) * bogo(sb+N+j,sb+N+k)

            !-------------------------------------------------------------------
            !                                   U   f        V^{\dagger}     
            ! Note the minus sign due to the hidden time-reversal!
            kappa(si+i,si+j)  = kappa(si+i,si+j) -             &
            &           config(sb+  k)*bogo(sb+  i,sb+N+k) * bogo(sb+N+j,sb+N+k)
            !                                   V^{*}(1 - f) U^{T} 
            kappa(si+i,si+j)  = kappa(si+i,si+j) +             &
            &           config(sb+N+k)*bogo(sb+N+i,sb+N+k) * bogo(sb  +j,sb+N+k)
          enddo
        enddo
      enddo
   
      si = si +   N
      sb = sb + 2*N
    enddo

  end subroutine PairingMatrices

  subroutine FindFermi(H, blocks, targetparticles, config, Bogo, Eqp,  lambda, &
  &                    maxhfbiter, blocktype, blockconf)
      !-------------------------------------------------------------------------
      ! Subroutine that diagonalizes the HFB hamiltonian (repeatedly) to find  
      ! the correct Fermi energy that fixes the average number of particles.
      ! The routine only solves this for one particular isospin.
      !
      ! Input
      !   H        : HFB hamiltonian, without Fermi energy
      !   blocks   : Sizes of the symmetry blocks that can be used to simplify 
      !              the problem.
      !   particles: Average number of particles to target. 
      !   lambda   : Initial guess for the Fermi energy
      !   maxhfbiter: Maximum number of iterations to perform
      !
      ! Output
      !   config   : Configuration matrix of the final solution
      !   Eqp      : Quasiparticle energies of the final solution
      !   Bogo     : Bogoliubov transformation that diagonalizes H
      !   Lambda   : Final fermi energy
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in)    :: H(:,:), targetparticles
      real(KIND=dp), intent(out)   :: config(:), Bogo(:,:), Eqp(:)
      real(KIND=dp), intent(inout) :: lambda
      integer, intent(in)          :: blocks(4), maxhfbiter, blocktype
      integer, intent(in)          :: blockconf(:)

      real(KIND=dp), allocatable   :: eigen(:), work(:), A(:,:)
      real(KIND=dp)                :: df, dn(2), particles

      integer :: iter, sb, si, N, B, i

      ! Initialization
      df = 0 ; dn = 0.0
      allocate(work(2*sum(blocks)))
      allocate(eigen(2*sum(blocks)))
      do iter=1, maxHFBiter

        !-----------------------------------------------------------------------
        ! a) Diagonalization of the HFB Hamiltonian by block. 
        si = 0 ; sb = 0
        do B=1,4
          N = blocks(B) ; if(N .eq. 0) cycle
  
          ! Construct the blocks of H including the Fermi energy
          allocate(A(2*N,2*N)) ; A = 0
          
          A = H(sb+1:sb+2*N, sb+1:sb+2*N)
          do i=1,N
            A(i  ,i  ) = A(i  , i  ) - lambda
            A(i+N,i+N) = A(i+N, i+N) + lambda
          enddo
                          
          ! Diagonalize every block
          call diagon (A,2*N,2*N,Bogo(sb+1:sb+2*N, sb+1:sb+2*N),               &
          &                       eigen(sb+1:sb+2*N),work)

          Eqp(si+1:si+N) = eigen(sb+N+1:sb+2*N)
          deallocate(A)
          ! Indices for the next block
          si = si +   N
          sb = sb + 2*N
        enddo
        !-----------------------------------------------------------------------
        ! b) We construct the configuration matrix that was asked for
        config = ConstructConfiguration(Bogo, Eqp, blocks, blocktype, blockconf) 
        !-----------------------------------------------------------------------
        ! c) Count the total number of particles that we have.
        si = 0 ; sb = 0
        particles   = 0
        do B=1,4                
            N = HFBsizes(B) ;  if(N .eq. 0) cycle 
            !-------------------------------------------------------------------
            ! Calculate the number of particles in here  
            ! Sum_i rho_ii =  Sum_ii   U   f U^{\dagger} + V^{*}(1 - f)V^{T}
            !          sum_(ij>N) f_(j) V^*_ij V^T_ji = sum_ij f_(j) V^*_ij V_ij
            !        + sum_(ij<N) f_(j) U^*_ij U^T_ji = sum_ij f_(j) U^*_ij U_ij         
            do i=1,N
                particles = particles                                          &
                &         +   config(sb+N+i)*sum(bogo(sb+N+1:sb+2*N, sb+N+i)**2)            
            enddo
            do i=1,N
                particles = particles                                          &
                &         +   config(sb  +i)*sum(bogo(sb  +1:sb+  N, sb+N+i)**2)            
            enddo

            ! indices
            si = si +  N
            sb = sb +2*N
        enddo
        ! When Time-reversal is conserved, we need an extra factor of two
$TR     particles = 2 * particles            

        ! Return if we do not want to readjust the Fermi energy
        if(MaxHFBiter.eq.1) return
        
        !-----------------------------------------------------------------------
        ! d) Readjust the Fermi energy based on the number of particles.
        !    We use the secant method.
        dn(2) = dn(1)
        dn(1) = particles - targetparticles
      
!        print ('(i3,f8.3,f8.3,3x, f4.1, f8.3, 3x, e12.3)'),                   &
!        &  iter, lambda, particles, targetparticles, df, dn(1)

        if(abs(dn(1)).lt.pairing_prec) return
          
        if(iter.eq.1) then
          ! We try lambda + 0.1 for the first iteration
          lambda = lambda + 0.1
          df     =          0.1
        else
          df     = - dn(1) * df/(dn(1) - dn(2))

          if(abs(df).gt.1.0) df = 0.1 * df/abs(df)
          lambda = lambda + df
        endif
      enddo

      deallocate(work)  ; deallocate(eigen)
  end subroutine FindFermi

  function ConstructHFBHamil(sphamil, gaps) result(H)
    !---------------------------------------------------------------------------
    ! Construct the HFB hamiltonian
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in) :: sphamil(:,:), gaps(:,:)
    real(KIND=dp), allocatable  :: H(:,:)
    
    integer :: N, i
    
    N = size(sphamil,1)
    allocate(H(2*N,2*N)) ; H = 0
    
    H(1:N, 1:N)         =  sphamil
    H(N+1:2*N, N+1:2*N) = -sphamil
    
    H(1:N, N+1:2*N)     =  gaps
    H(N+1:2*N, 1:N)     =  gaps 
    
  end function ConstructHFBHamil 

  subroutine calcHFBgaps(Fermi)
    !---------------------------------------------------------------------------
    ! Calculates the HFB gaps for use in the HFB solver.
    !
    ! Fermi is a dummy argument in this routine, since we need that argument 
    ! for the BCS solver.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! When time-reversal is conserved, the matrix calculated, is not the 
    ! full pairing gap matrix. Instead it is only part of it:
    !
    !                (  0         Delta )
    ! Delta_full  =  (                  )
    !                ( -Delta^T     0   ) 
    !
    ! where the separation in blocks is effected by a conserved antihermitian,
    ! linear symmetry. In EV8, this would have been z-signature.
    !
    ! For this smaller matrix Delta, not all of the single-particle 
    ! wavefunctions are actually in storage. We hence calculate only matrix 
    ! elements where
    !
    !               Delta_{i\bar{j}}
    !
    ! where \bar{j} is the time-reversed partner of the wavefunction j that is
    ! in storage. Note that Delta_{i\bar{j}} = Delta{ \bar{i} j}, hence the 
    ! smaller matrix is symmetric.
    !
    ! The implicit time-reversal in this action is included in the code 
    ! generated for the subroutine delta_action in functional.f90. 
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    ! When time-reversal is not conserved, it is indeed the full matrix that 
    ! is stored. This full matrix is antisymmetric, not symmetric!
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Fermi(2)
    integer                   :: wave1, wave2, iso, si,  B, N, inda, indb
    real(KIND=dp)             :: deltapsi(mv,4)
    
    if(.not.associated(Delta_action_HFB)) stop
    !---------------------------------------------------------------------------
    ! Use the delta_action to calculate the elements in the gaps
    si      = 0
    HFBgaps = 0

    ! Loop over the first of all blocks linked by the antihermitian symmetry
    ! 1,3,5,7
    do B=1,8,2   
      N   = HFBlocks(B)

      iso = -1
      if(B .gt. 4) iso = 1
    
      do wave1=1,N
        
        ! If time-reversal is not conserved:
        !     The first index comes from the second symmetry block.
        inda = si + wave1 + N
        ! If time-reversal is conserved:
        !     The first index comes from the first block, and an implicit
        !     time-reversal operation is performed in delta_action_HFB.
$TR        inda = si + wave1

        deltapsi = delta_action_HFB(  hfpsi(:,:,  inda),                       &
        &                            hfdpsi(:,:,:,inda),                       &
        &                           hfddpsi(:,:,:,inda),                       &
        &                          hfdddpsi(:,:,:,inda),                       &
        &                        sx(:,inda), sy(:,inda), sz(:,inda),iso,.false.)
        
        do wave2=wave1,N
          ! The second index is always in the first block. 
          indb = si + wave2 

          ! Mystery factor 0.5 in here
          HFBgaps(indb,inda) = 0.5*sum(hfpsi(:,:,indb)*deltapsi)*dv *          &
          &                                Pcutoffs(inda)*Pcutoffs(indb)

          ! The full matrix Delta is antisymmetric...
$NTR      HFBgaps(inda,indb) =  - HFBgaps(indb,inda)
          ! ... but the stored matrix is symmetric when time-reversal is 
          ! conserved.
$TR       HFBgaps(inda,indb) = HFBgaps(indb,inda)

        enddo
      enddo
      si = si + N
  enddo
  end subroutine calcHFBgaps

  subroutine PrintHFBconvergence(rho_pairing, kappa_pairing)
    !----------------------------------------------------------------------------
    ! Prints out some convergence info on the HFB subproblem.
    !   a) size of rho*rho - rho + kapppa * kappa^T
    !   b) size of rho*kappa - kappa * rho
    ! These should be small at convergence.
    !----------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: rho_pairing(:,:), kappa_pairing(:,:)
    real(KIND=dp)              :: test1(4), test2(4)
    real(KIND=dp), allocatable :: A(:,:) , r(:,:), k(:,:)
    integer                    :: N,i,j, P, B, si

    1 format (' HFB convergence:   (N,+)    (N,-)    (P,+)    (P,-)')
    2 format ('  r^2-r+k*k^T    = ',  4es9.2)
    3 format ('  r*k-k*r        = ',  4es9.2)

    print *
    print 1
    
    si = 0
    do B=1,8
        N = HFBlocks(B) ; if (N.eq. 0) cycle
 
        allocate(A(N,N), r(N,N), k(N,N))
        r = rho_pairing(si+1:si+N,si+1:si+N)
        k = kappa_pairing(si+1:si+N,si+1:si+N)
        
        ! A = rho^2 - rho + kappa * kappa^T
        A = matmul(r,r) - r + matmul(k, transpose(k))
        test1(B) = sqrt(sum(A**2))

        ! A = rho * kappa - kappa * rho
        A = matmul(r, k) - matmul(k,r)

        test2(B) = sqrt(sum(A**2))
        si = si + N
        deallocate(A,r,k)
    enddo

    print 2, test1
    print 3, test2

  end subroutine PrintHFBconvergence

  subroutine Canonical(rho_pairing, kappa_pairing, rho_can, kappa_can,         &
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
    do B=1,8
      N =HFBsizes(B) ;  if(N .eq. 0) cycle 
      
      allocate(tmp(N,N)) 
      
      tmp = rho_pairing(si+1:si+N, si+1:si+N)
      
      ! Diagonalize rho in this block
      call diagon(tmp,N,N,rhotransfo(si+1:si+N,si+1:si+N),rho_can(si+1:si+N),work)

!      ! DEBUG
!      do i=1,N
!         print ('(20f10.3)'), rho_pairing(si+i,si+1:si+N)        
!      enddo
!      print *
!      print ('(a6, i3, 20f10.3)'), 'Occ B=',B, rho_can(si+1:si+N)
!      print *

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
    do B=1,8
    
      N = HFBsizes(B) ;  if(N .eq. 0) cycle
      
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
      N = HFBlocks(B)  ;  if(N .eq. 0) cycle 
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

   subroutine Identify(i, blocks, bi, qblock)
    !---------------------------------------------------------------------------
    ! Identifies both the symmetry block(qblock) and index in said symmetry 
    ! block (bi), based on the index i in the HF basis. 
    !---------------------------------------------------------------------------
    integer, intent(in)  :: i, blocks(4)
    integer, intent(out) :: bi, qblock    
        
    integer :: sb, N, B

    sb = 0; bi = 0; qblock = 0
    do B=1,4
        N = Blocks(B) ; if (N.eq.0) cycle

        if( i .gt. sb .and. i.le.sb+N) then
            qblock = B
            bi     = i - sb
            return
        endif
        sb = sb + N
    enddo
   end subroutine Identify

!===============================================================================
!  Never to be used function to define an interface for delta_action
!===============================================================================   
   function delta_action_dummy(psi,dpsi,ddpsi, dddpsi, sx,sy,sz,iso, onthefly) &
                                                                result(deltapsi)
      !-------------------------------------------------------------------------
      ! Dummy function to allow this module to acces the functional.f90 module 
      ! to acces the information on the acces of deltas.
      !-------------------------------------------------------------------------
      
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
 
end module
