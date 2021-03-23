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

  use geninfo
  use wavefunctions
  use pairingcutoffs

  use HFB_direct
  use HFB_gradient

  implicit none

  !-----------------------------------------------------------------------------
  ! Matrix of all qp energies
  real(KIND=dp), allocatable :: HFBGaps(:,:)
  ! Maximum amount of iterations for finding a Fermi energy
  integer :: maxHFBiter  = 200
  ! Dispersion of the particle number
  real(KIND=dp) :: HFBdispersion(2)
  ! Integer indexing the conjugate partners in a HFB calculation
  integer, allocatable :: conjugp(:)
  ! Pointer to relink procedures
  procedure(delta_action_dummy), pointer :: delta_action_HFB
  ! HFB gauge parameter
  real(KIND=dp) :: HFBGauge(2) = 0.0
  ! Norm of the gradient for the gradient solver
  real(KIND=dp) :: HFBgradnorm(4)    = 0.0
  real(KIND=dp) :: oldHFBgradnorm(4) = 0.0
  real(KIND=dp) :: expectedDE(2)     = 0.0
  real(KIND=dp), allocatable :: prev_update(:,:)
  !---------------------------------------------------------------------------
  ! History of the pairing matrices, for mixing purposes.
  real(KIND=dp), allocatable ::  rho_history(:,:), kappa_history(:,:)
  real(KIND=dp), allocatable ::  configmatrix_history(:) 
  real(KIND=dp), allocatable ::  Bogoliubov_history(:,:)
  real(KIND=dp), allocatable ::  overrho, overkap

  integer, save :: effblocks(8) = 0


  interface
   function delta_action_dummy(psi,dpsi,ddpsi, dddpsi, sx,sy,sz,iso, onthefly) &
                                                                result(deltapsi)
      !-------------------------------------------------------------------------
      ! Dummy function to allow this module to acces the functional.f90 module 
      ! to acces the information on the acces of deltas.
      !-------------------------------------------------------------------------
      real*8, intent(in)    :: psi(:,:)  
      real*8, intent(inout) :: dpsi(:,:,:),ddpsi(:,:,:), dddpsi(:,:,:)
      integer, intent(in)   :: sx(:),sy(:),sz(:),iso
      real*8, allocatable   :: deltapsi(:,:)
      logical, intent(in)   :: onthefly
   end function
  end interface
  !-----------------------------------------------------------------------------
  ! Which routine to use to find the Fermi energy
  procedure(FindFermi_Brent), pointer  :: FindFermi

contains
    
  subroutine solvepairing_HFB_direct(sphamil, gaps, fermi, Bogoliubov,         & 
  &                           rho_pairing, kappa_pairing, configmatrix,        & 
  &                           qpenergies, BlockType,Blockindices,              &
  &                           blocklowest, blocked_qps, ifail)

    !---------------------------------------------------------------------------
    ! Driver routine for the solving of the HFB equations in a direct fashion.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(inout) :: Fermi(2)
    real(KIND=dp), intent(inout) :: Bogoliubov(:,:)
    real(KIND=dp), intent(inout) :: kappa_pairing(:,:), rho_pairing(:,:)
    real(KIND=dp), intent(inout) :: configmatrix(:), qpenergies(:)
    ! Quantities for the HFB hamiltonian
    real(KIND=dp), intent(in)    :: sphamil(:,:),gaps(:,:)
    real(KIND=dp)                :: HFBHamil(2*nwt, 2*nwt)
    
    ! Configuration for the blocking
    integer, intent(in)          :: Blockindices(:)
    integer, intent(in)          :: BlockType
    integer, intent(out)         :: ifail
    character(len=2), intent(in) :: BlockLowest(:)
    integer, allocatable         :: neutron_block(:), proton_block(:)
    integer, allocatable         :: blocked_qps(:), p_blocked(:), n_blocked(:)

    
    integer                     :: si, sb, N, N2, B,  wave1, it, i, np, nn
    integer                     :: n_ind, p_ind, NB
    real(KIND=dp), allocatable :: chi(:,:), temp(:,:), tempc(:), tempqe(:)
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
$NTR        if(blocktype.eq.4) neutron_block(2) = neutron_block(2) + 1 
        case('n-')
            neutron_block(3) = neutron_block(3) + 1 
$NTR        if(blocktype.eq.4) neutron_block(4) = neutron_block(4) + 1 
        case('p+')
            proton_block(1)  = proton_block(1)  + 1
$NTR        if(blocktype.eq.4) proton_block(2)  = proton_block(2)  + 1
        case('p-')
            proton_block(3)  = proton_block(3)  + 1
$NTR        if(blocktype.eq.4) proton_block(4)  = proton_block(4)  + 1
        case('n0')
            neutron_block(5) = neutron_block(5) + 1
        case('p0')
            proton_block(5)  = proton_block(5)  + 1
        end select
      enddo
    end select

    !---------------------------------------------------------------------------
    ! a) We construct the HFB-hamiltonian for every block. 
    !    We pass in everything to the routine by PAIRS of blocks
    si      = 0 ; sb = 0
    do B=1,8,2
      N  = HFBlocks(B)    ! Size of the first partner block
      N2 = HFBlocks(B+1)  ! Size of the second partner block
      
      it = 1 ; if (B .gt. 4) it = 2

      HFBHamil(sb+1:sb+2*N+2*N2, sb+1:sb+2*N+2*N2) = ConstructHFBHamil(        &
      &                           sphamil(si+1:si+N+N2,si+1:si+N+N2),          &
      &                           gaps(si+1:si+N+N2,si+1:si+N+N2), N, N2,      &
      &                           rho_history(si+1:si+N+N2,si+1:si+N+N2),      &
      &                           kappa_history(si+1:si+N+N2,si+1:si+N+N2),    &
      &                           HFBgauge(it))  

      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
    enddo

    !---------------------------------------------------------------------------
    !    We repeatedly diagonalize the matrix to find a suitable Fermi energy, 
    !    for every isospin. 
    call FindFermi(HFBHamil(      1:2*nwn,      1:2*nwn), HFBlocks(1:4),       &
    &              neutrons, configmatrix(  1:2*nwn),                          &
    &              Bogoliubov(1:2*nwn, 1:2*nwn),                               &
    &              qpenergies(1:2*nwn),  Fermi(1), maxhfbiter,                 &
    &              blocktype, neutron_block, n_blocked, ifail)   

    call FindFermi(HFBHamil(2*nwn+1:2*nwt,2*nwn+1:2*nwt), HFBlocks(5:8),       &
    &              protons, configmatrix(2*nwn+1:2*nwt),                       &
    &              Bogoliubov(2*nwn+1:2*nwt, 2*nwn+1:2*nwt),                   &
    &              qpenergies(2*nwn+1:2*nwt),Fermi(2), maxhfbiter,             &
    &              blocktype, proton_block, p_blocked, ifail)     

    ! 
    if(allocated(blocked_qps)) deallocate(blocked_qps)
    NB = 0 ;  NN = 0 ; NP = 0
    if(allocated(n_blocked)) NN = size(n_blocked)
    if(allocated(p_blocked)) NP = size(p_blocked)
    NB = NP + NN

    allocate(blocked_qps(NB))
    if(allocated(n_blocked))  blocked_qps(   1:NN) = n_blocked
    if(allocated(p_blocked)) then
       ! We need to offset stuff by the number of neutron qps
       do i=1, NP
         blocked_qps(NN+i) = p_blocked(i) + sum(HFBlocks(1:4))
       enddo
    endif
    !---------------------------------------------------------------------------
    !    Reorganise the matrices into the block-form used by the rest of the 
    !    program.
    call reorganise_matrices(Bogoliubov,qpenergies,configmatrix)
  
    !---------------------------------------------------------------------------
    !    Construct the density and anomalous density matrices, based on the 
    !    configmatrix and the Bogoliubov transformation
    call PairingMatrices(configmatrix, bogoliubov, rho_pairing, kappa_pairing)

    HFBdispersion = calc_dispersion_HFB(rho_pairing, kappa_pairing)
  end subroutine solvepairing_HFB_direct

  subroutine solvepairing_HFB_gradient(stepsize, sphamil, gaps, fermi, Bogo,   & 
  &                          rho_pairing, kappa_pairing,configmatrix,          & 
  &                          qpenergies)
    !---------------------------------------------------------------------------
    ! Driver routine for solving the HFB equations by gradient stepping in the 
    ! Bogoliubov manifold.
    !---------------------------------------------------------------------------
    use HFB_gradient

    real(KIND=dp), intent(inout) :: Fermi(2)
    real(KIND=dp), intent(inout) :: Bogo(:,:)
    real(KIND=dp), intent(inout) :: kappa_pairing(:,:), rho_pairing(:,:)
    real(KIND=dp), intent(inout) :: configmatrix(:), qpenergies(:) 
    real(KIND=dp), intent(in)    :: sphamil(:,:),gaps(:,:),stepsize

    real(KIND=dp)                :: oldnorm(4)
    real(KIND=dp), allocatable   :: tempEqp(:), tempdisp(:), tempBogo(:,:)
    integer, allocatable         :: indices(:) 

    integer :: si, sb, B, N, N2, T,i, Np, Nm, ind, ind2
  
    if(.not.allocated(rho_history)) then
      allocate(rho_history(nwt,nwt))            ; rho_history   = 0.0
      allocate(kappa_history(nwt,nwt))          ; kappa_history = 0.0
      allocate(configmatrix_history(2*nwt))     ; configmatrix_history = 0.0
      allocate(Bogoliubov_history(2*nwt, 2*nwt)); Bogoliubov_history = 0.0
    endif

    ! Saving the history
    !rho_history          = rho_pairing
    !kappa_history        = kappa_pairing
    !configmatrix_history = configmatrix
    Bogoliubov_history   = Bogo

    ! Guess a new Fermi energy if none is there
    if(all(Fermi.eq.0.0))   Fermi = -5

    allocate(tempEqp(nwt))          ; tempEqp  = 0.0d0
    allocate(tempBogo(2*nwt, 2*nwt)); tempBogo = 0.0d0

    if(all(effblocks.eq.0)) then    
      effblocks = HFblocks
      sb = 0
      do B=1,8,2
        N = HFBlocks(B)   ; if(N.eq.0) cycle
        N2= HFBlocks(B+1)
        T = N+N2

        ! Counting the number of positive signature states
        Np = 0
        do i=1,N
          if(configmatrix(sb+T+i).eq.1.0d0) Np = Np+1
        enddo
        do i=N+1,T
          if(configmatrix(sb+T+i).eq.0.0d0) Np = Np+1
        enddo
        ! Getting the indices right
        allocate(indices(T)) ; ind = 0 ; ind2 = 0
        do i=1,N
          if(configmatrix(sb+T+i).eq.1.0d0) then
            ind          = ind + 1
            indices(ind) = i
          else
            ind2             = ind2 + 1
            indices(Np+ind2) = i
          endif
        enddo
        do i=N+1,T
          if(configmatrix(sb+T+i).eq.1.0d0) then
            ind2             = ind2 + 1
            indices(Np+ind2) = i
          else
            ind              = ind + 1
            indices(ind )    = i
          endif
        enddo

        do i=1,T
          if(configmatrix(sb+T+indices(i)).eq.1.0d0) then
            tempBogo(sb+1:sb+2*T ,sb+T+i)  = Bogo(sb+1:sb+2*T,sb+T+indices(i))  
          else
            tempBogo(sb+1  :sb+  T,sb+T+i) = Bogo(sb+T+1:sb+2*T,sb+T+indices(i))    
            tempBogo(sb+T+1:sb+2*T,sb+T+i) = Bogo(sb+  1:sb+  T,sb+T+indices(i))    
          endif
        enddo    
        deallocate(indices)

        effBlocks(B)   = Np
        effBlocks(B+1) = T - Np

        sb = sb + 2*T
      enddo
    else
        tempBogo = Bogo
    endif

    !---------------------------------------------------------------------------
    ! Stepping for the neutrons
    call gradient_step(sphamil(1:nwn,1:nwn),gaps(1:nwn,1:nwn),                 & 
    &                  effblocks(1:4), neutrons,                               &
    &                  tempBogo(1:2*nwn,1:2*nwn),                              &
    &                  tempEqp(1:nwn), stepsize,Fermi(1),                      &
    &                  HFBgradnorm(1:2), expectedDE(1),                        &
    &                  rho_pairing(1:nwn,1:nwn),kappa_pairing(1:nwn,1:nwn),    & 
    &                  0.0d0,1)
    ! and for the protons
    call gradient_step(sphamil(nwn+1:nwt,nwn+1:nwt),gaps(nwn+1:nwt,nwn+1:nwt), & 
    &                 effblocks(5:8),protons,                                  &
    &                 tempBogo(2*nwn+1:2*nwt,2*nwn+1:2*nwt),                   &
    &                 tempEqp(nwn+1:nwt),stepsize,Fermi(2),                    &
    &                 HFBgradnorm(3:4), expectedDE(2),                         & 
    &                 rho_pairing(nwp+1:nwt,nwp+1:nwt),                        & 
    &                 kappa_pairing(nwn+1:nwt,nwn+1:nwt),                      & 
    &                 0.0d0,1)
  
    Bogo         = tempBogo
    sb = 0
    configmatrix = 0.0d0
    do B=1,8,2
      N = HFBlocks(B)   ; if(N.eq.0) cycle
      N2= HFBlocks(B+1)
      T = N+N2

      configmatrix(sb+1:sb+T)     = 0.0d0
      configmatrix(sb+T+1:sb+2*T) = 1.0d0
      sb = sb + 2*T
    enddo

    call PairingMatrices(configmatrix, bogo, rho_pairing, kappa_pairing)
    
    ! Final organisation
    sb = 0 ; si = 0
    do B=1,8,2
      N = HFBlocks(B) ; N2 = HFblocks(B+1)       
      T = N + N2        
      do i=1,T
        qpenergies(sb  +i) = -tempEqp(si+T-i+1)
        qpenergies(sb+T+i) =  tempEqp(si+i)
      enddo
      si = si +   T
      sb = sb + 2*T
    enddo
    
    deallocate(tempEqp, tempBogo)
  end subroutine solvepairing_HFB_gradient

  function calc_dispersion_HFB(rho,kappa) result(dispersion)
    !---------------------------------------------------------------------------
    ! We calculate the dispersion of the particle number
    ! 
    !  <N^2> = sum_{ab} <a^{\dagger}_{a} a_{a} a^{\dagger}_{b} a_{b} > 
    !        = sum_{ab} rho_{aa} rho_{bb} 
    !                +  rho_{ab} ( 1 - rho^*_{ab})
    !                +  kappa_{ab}^* kappa_{ab}
    !
    ! So <N^2> - <N>^2 = Tr(rho ( 1 -rho)) +  Tr(kappa * kappa^{\dagger})
    ! 
    !-------------------------------------------------------------------------
    ! Note, that at T = 0, we have that (kappa * kappa^{\dagger}) = rho(1-rho).
    ! So in that case, we have 
    !  < Delta N^2 > = < N^2 > - <N>^2 = 2 * Tr(rho(1-rho))
    ! which is the old formula from EV8, CR8, etc...
    !
    ! We implement however the formula above, since this is the one that 
    ! correctly generalizes to T != 0.
    !-------------------------------------------------------------------------
    real(KIND=dp), intent(inout) :: rho(:,:), kappa(:,:)
    real(KIND=dp), allocatable   :: chi(:,:)
    real(KIND=dp)                :: dispersion
    integer                      :: N, B, si, i, it

    si            = 0
    dispersion = 0.0
    allocate(chi(nwt, nwt))
    do B=1,8
      N = HFBlocks(B)  ;  if(N .eq. 0) cycle 

      it  = 1
      if(B .gt. 4) it = 2
      
      ! rho squared
      chi = matmul(rho(si+1:si+N, si+1:si+N), rho(si+1:si+N, si+1:si+N) )    
      !                                       Tr rho - Tr rho^2
      do i=1,N
          HFBdispersion(it) = HFBdispersion(it) + rho(si+i, si+i)     &
          &                                     - chi(i,i)            &
          &                                     + kappa(si+i, si+i)**2
      enddo
      si = si + N
    enddo
    deallocate(chi)

    ! Time-reversal introduces a factor of two
$TR    dispersion = 2 * dispersion 

  end function calc_dispersion_HFB

  subroutine mix_pairing(mix, rho, kappa) 
    !---------------------------------------------------------------------------
    ! Linearly mix rho and kappa
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: mix
    real(KIND=dp), intent(inout) :: rho(:,:), kappa(:,:)

    rho   =  mix * rho   + (1.0d0-mix) * rho_history
    kappa =  mix * kappa + (1.0d0-mix) * kappa_history

  end subroutine mix_pairing 

  subroutine PairingMatrices(config, bogo, rho, kappa)
    !---------------------------------------------------------------------------
    ! Calculate the matrices rho and kappa for a given configuration matrix 
    ! and Bogoliubov transformation.
    !
    !    rho   = U   f U^\dagger + V^* (1 - f) V^T
    !    kappa = U   f V^\dagger + V^* (1 - f) U^T  
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !  Notice how this routine only uses the last half of the columns of the  
    !  Bogoliubov transformation.
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: config(:), Bogo(:,:)
    real(KIND=dp), intent(out):: rho(:,:), kappa(:,:)
    integer                   :: si, sb, B, i,j,k, N, N2, column

    si = 0 ; sb = 0 ; kappa = 0.0d0 ; rho = 0.0d0
    do B=1,8,2
      N = HFBlocks(B)   ;  if(N .eq. 0) cycle 
      N2= HFBlocks(B+1)
      do i=1,N+N2
        do j=1,N+N2
          do k=1,N+N2
            column = sb+N+N2+k
            !-------------------------------------------------------------------
            !                                   U      f  U^{\dagger}
            rho(si+i,si+j)  = rho(si+i,si+j) +                 &
            &     config(sb+  k)*bogo(sb+  i,column) * bogo(sb+  j,column)
            !                                   V^* (1-f) V^{T}
            rho(si+i,si+j)  = rho(si+i,si+j) +                 &
            &     config(column)*bogo(sb+N+N2+i,column) * bogo(sb+N+N2+j,column)

            !-------------------------------------------------------------------
            !                                   U   f        V^{\dagger}     
            ! Note the minus sign due to the hidden time-reversal!
$TR            kappa(si+i,si+j)  = kappa(si+i,si+j) -             &
$TR            &     config(sb+  k)*bogo(sb+  i,column) * bogo(sb+N+N2+j,column)


$NTR            kappa(si+i,si+j)  = kappa(si+i,si+j) +             &
$NTR            &     config(sb+  k)*bogo(sb+  i,column) * bogo(sb+N+N2+j,column)
            !                                   V^{*}(1 - f) U^{T} 
            kappa(si+i,si+j)  = kappa(si+i,si+j) +             &
            &     config(column)*bogo(sb+N+N2+i,column) * bogo(sb  +j,column)
          enddo
        enddo
      enddo
   
      si = si +   N +  N2
      sb = sb + 2*N +2*N2
    enddo

  end subroutine PairingMatrices

  function ConstructHFBHamil(sphamil, gaps, N, N2, r, k, gauge) result(H)
    !---------------------------------------------------------------------------
    !  Construction of the HFB Hamiltonian in the form
    !   
    !  ( h+     d+-  0     0   )
    !  ( d+-   -h-   0     0   )   = H_constructed 
    !  ( 0      0    h-    d-+ )
    !  ( 0      0    d-+  -h+  )
    !
    ! which is a reshuffling of 
    !
    !  ( h+     0    0    d+- )
    !  ( 0      h-   d-+  0   )    = H_formal
    !  ( 0      d-+ -h+   0   )
    !  ( d+-    0    0   -h-  )
    !
    ! H_constructed can be diagonalised in blocks. 
    !
    ! The + and - blocks are characterized by a conserved, linear, antihermitian
    ! symmetry. The usual suspect is R_z, z-signature symmetry. Their dimensions
    ! are N and N2. Note that this routine should trivially work correctly in
    ! the case that N2=0, i.e. if there is no separation in blocks because 
    ! either
    !           (i)  the symmetry is not conserved
    !           (ii) an antilinear, antihermitian, symmetry is conserved that
    !                guarantees the equality of both blocks. (The usual 
    !                suspect is time-reversal). 
    !   
    ! Note that all conserved symmetries that are easier to handle, in 
    ! particular linear hermitian ones (such as parity) are handled outside
    ! of this routine. 
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! This routine also gives you the option to add the generalised density
    ! matrix to the HFB Hamiltonian with a gauge parameter. 
    !
    !   H => H  + alpha * R
    !
    !   R = (  rho         kappa)
    !       (- kappa^*    1 - rho^*)
    !
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !
    ! Input:   
    !   sphamil: single-particle hamiltonian; in block-structure, i.e.
    !   
    !               h = ( h+ 0 )
    !                   ( 0  h-)
    !
    !  gaps    : pairing gaps; in block-structure, i.e.
    !               d = ( 0   d+-)
    !                   (d-+  0  )
    !
    !  N, N2   : sizes of the respective blocks
    !
    !  r, k    : rho and kappa pairing matrices for the addition of the 
    !            generalized density matrix. Should also be in correct
    !            block structure.
    !   
    !  gauge   : real parameter alpha
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)   :: sphamil(:,:), gaps(:,:), r(:,:), k(:,:)
    real(KIND=dp)               :: gauge
    real(KIND=dp), allocatable  :: H(:,:), genR(:,:)
    integer, intent(in)         :: N, N2
    integer                     :: T, i

    T = N + N2
    allocate(H(2*T,2*T), genR(2*T, 2*T)) 
    H = 0 ; genR = 0

    genR(  1:N   ,  1:N)    = r(1:N, 1:N)
    genR(  1:N   , N+1:2*N) =-k(1:N, 1:N)
    genR(N+1:2*N , 1:N)     =-k(1:N, 1:N)

    genR(N+1:2*N   ,N+ 1:2*N)=-r(1:N, 1:N)
    do i=1, N
        genR(N+i, N+i) = genR(N+i, N+i) + 1
    enddo
    !---------------------------------------------------------------------------
    ! SP hamil
    ! Block 1
    H(  1:N   ,  1:N)    = &
    &                sphamil(  1:N   ,   1:N)    + gauge * r(  1:N   ,   1:N)
    H(N+1:N+N2,N+1:N+N2) = &
    &               -sphamil(N+1:N+N2, N+1:N+N2) - gauge * r(N+1:N+N2, N+1:N+N2)
    do i=1, N2
        H(N+i, N+i) =  H(N+i, N+i) + gauge
    enddo

    ! Block 2
    H(N+  N2+1:  N+2*N2,N+  N2+1:  N+2*N2) = & 
    &                 sphamil(N+1:N+N2,N+1:N+N2)+  gauge * r(N+1:N+N2, N+1:N+N2)
    H(N+2*N2+1:2*N+2*N2,N+2*N2+1:2*N+2*N2) = & 
    &                -sphamil(  1:N   ,   1:N)  -  gauge * r(  1:N  ,    1:N)
    do i=1, N
        H(N+2*N2+i,N+2*N2+i) =  H(N+2*N2+i,N+2*N2+i) + gauge
    enddo

    !---------------------------------------------------------------------------
    ! Gaps
    ! Block 1
    H(  1:N   ,N+1:2*N)  = gaps(  1:N   , N2+1:N2+N) + gauge * k(1:N, N2+1:N2+N)
    H(N+1:2*N ,  1:N  )  = gaps(  1:N   , N2+1:N2+N) + gauge * k(1:N, N2+1:N2+N)

    ! Block 2
    H(2*N   +1:2*N+N2  ,2*N+N2+1:2*N+2*N2) = &
    &                            gaps(N+1:N+N2, 1:N2) + gauge * k(N+1:N+N2,1:N2)
    H(2*N+N2+1:2*N+2*N2,2*N   +1:2*N+  N2) = & 
    &                            gaps(N+1:N+N2, 1:N2) + gauge * k(N+1:N+N2,1:N2)

$TR if(gauge.ne.0.0) then
$TR     print *, 'something is fishy with the HFBgauge when T is conserved.'
$TR     print *, 'Investigate sign of kappa'
$TR     stop
$TR endif 

  end function ConstructHFBHamil 

  subroutine calcHFBgaps(Fermi, stabfactor)
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
    real(KIND=dp), intent(in) :: Fermi(2), stabfactor(2)
    integer                   :: wave1, wave2, iso, si,  B, N, inda, indb, N2,T
    real(KIND=dp)             :: deltapsi(mv,4), val(2), stabfac
    
    val = Fermi ! To avoid the unused dummy argument warning from the compiler

    if(.not.associated(Delta_action_HFB)) stop

    !---------------------------------------------------------------------------
    ! Use the delta_action to calculate the elements in the gaps
    si      = 0
    HFBgaps = 0

    ! Loop over the first of all blocks linked by the antihermitian symmetry
    ! 1,3,5,7
    do B=1,8,2   
      N   = HFBlocks(B)  ; if(N .eq. 0) cycle
      N2  = HFBlocks(B+1)
      T   = N + N2
      iso = -1
      if(B .gt. 4) iso = 1
    
      do wave1=1,N
        
        ! If time-reversal is not conserved:
        !     The first index comes from the second symmetry block.
$NTR        inda = si + wave1 + N
        ! If time-reversal is conserved:
        !     The first index comes from the first block, and an implicit
        !     time-like symmetry operation is performed in delta_action_HFB.
$TR        inda = si + wave1

        deltapsi = delta_action_HFB(  hfpsi(:,:,  inda),                       &
        &                            hfdpsi(:,:,:,inda),                       &
        &                           hfddpsi(:,:,:,inda),                       &
        &                          hfdddpsi(:,:,:,inda),                       &
        &                        sx(:,inda), sy(:,inda), sz(:,inda),iso,.false.)
        
        do wave2=1,N
          ! The second index is always in the first block. 
          indb = si + wave2 

          ! Add the stabilisation factor
          ! (1 if the pairing functional is not stabilized)
          stabfac = 1 + stabfactor((iso+3)/2)

          HFBgaps(inda,indb) =     sum(hfpsi(:,:,indb)*deltapsi)*dv *stabfac

          ! The full matrix Delta is antisymmetric...
$NTR      HFBgaps(indb,inda) =  - HFBgaps(inda,indb)
          ! ... but the stored matrix is symmetric when time-reversal is 
          ! conserved; but this is not exploited at the moment!
!$TR       HFBgaps(indb,inda) = HFBgaps(inda,indb)

        enddo
      enddo
      !-------------------------------------------------------------------------      
      ! We have now calculated the gaps (without cutoffs) in the basis that 
      ! is currently in storage. This can either be the HF basis or not!
      if(.not.diagsphamil .and. allocated(HFtransfo)) then
        ! Transform to the HF basis
        HFBgaps(si+1:si+T,si+1:si+T) = &
        &                     matmul(transpose(HFtransfo(si+1:si+T,si+1:si+T)),&
        &                                          HFBgaps(si+1:si+T,si+1:si+T))
        HFBgaps(si+1:si+T,si+1:si+T) = &
        &    matmul(HFBgaps(si+1:si+T,si+1:si+T),HFtransfo(si+1:si+T,si+1:si+T))
      endif
      !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -     
      ! Multiply with the cutoff factors
      do wave1=1,N
$NTR    inda = si + wave1 + N
$TR     inda = si + wave1
        do wave2=1,N
          indb = si + wave2 
          HFBgaps(inda,indb) = HFBgaps(inda,indb)*Pcutoffs(inda)*Pcutoffs(indb) 
$NTR      HFBgaps(indb,inda) = HFBgaps(indb,inda)*Pcutoffs(inda)*Pcutoffs(indb) 
        enddo
      enddo
      if(.not.diagsphamil .and. allocated(HFtransfo)) then
        ! Transform back
        HFBgaps(si+1:si+T,si+1:si+T) = &
        & matmul((HFtransfo(si+1:si+T,si+1:si+T)),HFBgaps(si+1:si+T,si+1:si+T))
        HFBgaps(si+1:si+T,si+1:si+T) = &
        & matmul(HFBgaps(si+1:si+T,si+1:si+T),                                 &
        &                             transpose(HFtransfo(si+1:si+T,si+1:si+T)))
      endif

      si = si + N + N2
    enddo

  end subroutine calcHFBgaps

  subroutine PrintHFBconvergence(rho_pairing, kappa_pairing)
    !---------------------------------------------------------------------------
    ! Prints out some convergence info on the HFB subproblem.
    !   a) sqrt(sum( (rho*rho - rho + kapppa * kappa^T)**2)
    !   b) sqrt(sum( (rho*kappa - kappa * rho)**2))
    ! 
    ! These should be (numerically) vanishingly small at convergence, when we
    ! have solved the HFB problem consistently.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note that these numbers do not vanish when we do either
    !   (i)  finite-temperature calculations (they increase as T increases)
    !   (ii) Equal Filling-style blocking
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: rho_pairing(:,:), kappa_pairing(:,:)
    real(KIND=dp)              :: test1(8), test2(8)
    real(KIND=dp), allocatable :: A(:,:) , r(:,:), k(:,:)
    integer                    :: N, B, si, N2

    1 format (' HFB convergence:   (N,+)    (N,-)    (P,+)    (P,-)')
    2 format ('  r^2-r+k*k^T    = ',  4es9.2)
    3 format ('  r*k-k*r        = ',  4es9.2)
    4 format (' Grad. norm      = ',  4es9.2)

    print *
    print 1
    
    si = 0
    do B=1,8,2
        N = HFBlocks(B) ; if (N.eq. 0) cycle
        N2= HFBlocks(B+1)

        allocate(A(N+N2,N+N2), r(N+N2,N+N2), k(N+N2,N+N2))
        r = rho_pairing(si+1:si+N+N2  ,si+1:si+N+N2)
        k = kappa_pairing(si+1:si+N+N2,si+1:si+N+N2)
        
        ! A = rho^2 - rho + kappa * kappa^T
        A = matmul(r,r) - r + matmul(k, transpose(k))
        test1(B) = sqrt(sum(A**2))

        ! A = rho * kappa - kappa * rho
        A = matmul(r, k) - matmul(k,r)

        test2(B) = sqrt(sum(A**2))
        si = si + N + N2
        deallocate(A,r,k)
    enddo

    print 2, test1(1), test1(3), test1(5), test1(7)
    print 3, test2(1), test2(3), test2(5), test2(7)
    print 4, HFBgradnorm

  end subroutine PrintHFBconvergence

  subroutine Canonical(rho_pairing, kappa_pairing, rho_can, kappa_can,         &
  &                                             rhotransfo, kappatransfo, ifail)
    !---------------------------------------------------------------------------
    ! a) Diagonalize  Rho
    ! b) Canonicalize Kappa
    ! c) Return the diagonal elements of rho, and the offdiagonal elements of 
    !    kappa, as well as the transformation.
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)  :: rho_pairing(nwt,nwt)
    real(KIND=dp), intent(in)  :: kappa_pairing(nwt,nwt)
    real(KIND=dp), intent(out) :: rho_can(nwt), kappa_can(nwt)
    real(KIND=dp), intent(out) :: rhotransfo(nwt,nwt), kappatransfo(nwt,nwt)
 
    integer, intent(out) :: ifail
   
    real(KIND=dp), allocatable :: tmp(:,:), work(:)
    
    integer :: si, N, N2, B, i, lwork
$NTR integer :: j
    
    !---------------------------------------------------------------------------
    ! a) Diagonalize rho
    !
    ! Transforms as rho'  = D^dagger rho D
    !---------------------------------------------------------------------------
    si         = 0
    rhotransfo = 0 ; kappatransfo = 0
    do B=1,8
      N =HFBlocks(B) ;  if(N .eq. 0) cycle 
      
      allocate(tmp(N,N)) 
      
      tmp = rho_pairing(si+1:si+N, si+1:si+N)
      
      lwork = -1 ; allocate(work(1))
      call DSYEV( 'V', 'U', N, tmp, N, rho_can(si+1:si+N), work, lwork, ifail)
      lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
      call DSYEV( 'V', 'U', N, tmp, N, rho_can(si+1:si+N), work, lwork, ifail)
      deallocate(work)
  
      rhotransfo(si+1:si+N,si+1:si+N) = tmp

      if(ifail.ne.0) then
        print *, 'WARNING: diagon failed in subroutine Canonical.'
        print *, '         Problematic block B = ', B
        deallocate(tmp)
        return
      endif

      si = si + N
      deallocate(tmp)
    enddo
    
    ! Time-reversal
$TR    rho_can = 2*rho_can
    
    do i=1,nwt
$TR      if(rho_can(i).gt.2.0) rho_can(i) = 2.0
$NTR     if(rho_can(i).gt.1.0) rho_can(i) = 1.0
      if(rho_can(i).lt.0.0) rho_can(i) = 0.0
    enddo
    
    !---------------------------------------------------------------------------
    ! b) Bring kappa into canonical form
    !
    ! Transforms as kappa'  = D^T kappa D^*
    !---------------------------------------------------------------------------
    si = 0

$NTR if(.not. allocated(conjugp))  allocate(conjugp(nwt)) 
$NTR conjugp = 0

    do B=1,8,2
    
      N  = HFBlocks(B)   ;  if(N .eq. 0) cycle
      N2 = HFBlocks(B+1)

      allocate(tmp(N+N2,N+N2))
      tmp = kappa_pairing(si+1:si+N+N2, si+1:si+N+N2)   

      tmp = matmul(transpose(rhotransfo(si+1:si+N+N2, si+1:si+N+N2)), tmp)
      tmp = matmul(tmp,rhotransfo(si+1:si+N+N2, si+1:si+N+N2))

      ! With the assumption of time-reversal, the diagonal matrix elements in 
      ! this transformed kappa matrix are the matrix elements (i, ibar).
$TR      do i=1,N
$TR          kappa_can(si+i) = tmp(i  ,i)
$TR      enddo

      ! Without the assumption of time-reversal we cannot be guaranteed to know
      ! the canonical partners beforehand.
$NTR  conjugp(si+1:si+N+N2) = 0

$NTR  do i=1, N+N2
$NTR    do j=i+1, N+N2
$NTR      if(abs(rho_can(si+i) - rho_can(si+j)).lt.1d-12) then
$NTR        if( abs(tmp(i,j)).gt. 1d-14) then
$NTR            conjugp(si+i) =  si+j
$NTR            conjugp(si+j) =  si+i
$NTR            kappa_can(si+i) = tmp(i,j)
$NTR            kappa_can(si+j) = tmp(j,i)   
$NTR        endif
$NTR      endif    
$NTR    enddo
$NTR  enddo

      deallocate(tmp)
      si = si + N + N2
    enddo    
    
   end subroutine Canonical
   
   subroutine ConstructCanonicalBasis(Transfo)
    !---------------------------------------------------------------------------
    ! Transform the spwf wavefunctions in the HFBasis into Canbasis, with the 
    ! passed in Transfo. 
    !---------------------------------------------------------------------------

    integer                   :: wave1, wave2, B, N, si, wave3
    real(KIND=dp), intent(in) :: Transfo(nwt,nwt)
    
    if(.not.allocated(CanPsi)) then
      allocate(CanPsi(mv,4,nwt)) ; CanPsi      = 0.0
      allocate(Canenergies(nwt)) ; canenergies = 0.0
    endif

    si     = 0   
    CanPsi = 0.0 ; canenergies = 0.0
    do B=1,8
      N = HFBlocks(B)  ;  if(N .eq. 0) cycle 
      !-------------------------------------------------------------------------
      ! Apply the transformation in this symmetry block
      do wave1=1, N 
        do wave2=1,N
          CanPsi(:,:,si+wave1)  = CanPsi(:,:,si+wave1) +                       &
          &                     Transfo(si+wave2,si+wave1) * HFPsi(:,:,si+wave2) 
          
          if(.not. allocated(current_sph)) then
            canenergies(si+wave1) = canenergies(si+wave1) +                    &
            &          abs(Transfo(si+wave2,si+wave1)**2) *spenergies(si+wave2) 
          else
            do wave3=1,N
              canenergies(si+wave1) = canenergies(si+wave1) +                  &
            &         Transfo(si+wave2,si+wave1) *  Transfo(si+wave3,si+wave1) & 
            &                       * current_sph(si+wave2,si+wave3) 
            enddo
          endif
        enddo 
      enddo
      
      si = si +  N
      !-------------------------------------------------------------------------
    enddo
    
   end subroutine ConstructCanonicalBasis

   function construct_generalized_density(rho, kappa) result(R)
    !---------------------------------------------------------------------------
    ! Construct the generalized density matrix R from the matrices rho and kappa
    !
    !   R = (  rho     kappa    )
    !       (-kappa^*  1 - rho^*)
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in)  :: rho(:,:), kappa(:,:)
    real(KIND=dp), allocatable :: R(:,:)
    integer                    :: N,i 

    N = size(rho, 1)

    allocate(R(2*N, 2*N))

    R(  1:  N,   1:  N) = rho
    R(N+1:2*N, N+1:2*N) =-rho

    do i=1, N
        R(N+i, N+i) = R(N+i, N+i) + 1.0d0
    enddo
    
    R(  1:  N,N+1:2*N) = kappa
    R(N+1:2*N,  1:  N) =-kappa

   end function construct_generalized_density

   subroutine clean_HFB
    if(allocated(HFBgaps))  deallocate(HFBGaps)
   end subroutine clean_HFB

end module
