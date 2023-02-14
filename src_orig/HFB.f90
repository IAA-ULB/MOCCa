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
 ! This module is in charge of the solution of the HFB pairing subproblem 
 ! in the space spanned by the single-particle wavefunctions currently in 
 ! storage. This effort is chiefly handled by the routines
 ! 
 !     solvepairing_HFB_direct
 !     solvepairing_HFB_gradient
 !
 ! which make several calls each to the (more low-level) routines in the 
 ! HFB_direct and HFB_gradient modules, respectively. 
 !
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 !
 ! Hephaestos keywords:
 !
 !         TR : $TR
 !        NTR : $NTR
 ! PCONSERVED : $PCONSERVED
 !    PBROKEN : $PBROKEN
 !
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 !
 ! Implemented routines: 
 !
 !  - subroutine solvepairing_HFB_direct
 !  - subroutine solvepairing_HFB_gradient
 !  - subroutine reorganise_bogo_gradient
 !  - subroutine mix_pairing
 !  - subroutine pairingmatrices
 !  - subroutine calcHFBgaps
 !  - subroutine PrintHFBconvergence
 !  - subroutine Canonical
 !  - subroutine clean_HFB
 !  - subroutine correct_ordering_eqp
 !  - subroutine figure_out_blocking_structure
 !  - subroutine figure_out_blocking_structure_EFA
 !  - function   obtain_eqp
 !  - function   calc_dispersion_HFB
 !  - function   ConstructHFBHamil
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
  ! Angular momentum "expectation values" of the Bogoliubov quasiparticles
  real(KIND=dp), allocatable :: qp_J(:,:), qp_JTR(:,:), qp_JTI(:,:)
  ! Quasiparticle dispersions: <H^2> - E_qp^2
  real(KIND=dp), allocatable :: qpdispersions(:)
  !---------------------------------------------------------------------------
  ! History of the pairing matrices, for mixing purposes.
  real(KIND=dp), allocatable ::  rho_history(:,:), kappa_history(:,:)
  real(KIND=dp), allocatable ::  configmatrix_history(:) 
  real(KIND=dp), allocatable ::  Bogoliubov_history(:,:)
  !-----------------------------------------------------------------------------
  ! Cutoff parameter to judge whether or not levels are participating in the 
  ! pairing.
  real(KIND=dp), parameter            :: rho_cutoff = 1e-12

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
  &                          rho_pairing, kappa_pairing, configmatrix,         & 
  &                          qpenergies, BlockType,Blockindices,               &
  &                          blocklowest, blocked_qps, partner_qps,qp_overlaps,&
  &                          ifail)

    !---------------------------------------------------------------------------
    ! Driver routine for the solving of the HFB equations in a direct fashion,
    ! meaning by diagonalising the HFB Hamiltonian and explicitly constructing
    ! the desired Bogoliubov vacuum state out of its eigenvectors.
    !
    ! Input: 
    !    sphamil      : matrix of the single-particle hamiltonian
    !    gaps         : matrix of the pairing gaps
    !    fermi        : initial guess for the Fermi energies of both nucleon 
    !                   species
    !    Blocktype    : type of blocking procedure desired, to be passed into
    !                   the construct_configuration routine in HFB_direct
    !    Blockindices : when blocktype=1,3,5 contains the indices for the 
    !                   overlap-calculation for blocking. To be passed into 
    !                   construct_configuration.
    !    Blocklowest  : when blocktype=2,4,6 contains the type of excitations 
    !                   we want to build. To be passed into 
    !                   construct_configuration.
    !
    ! Ouput:
    !    Fermi        : final value obtained by the solver for the Fermi 
    !                   energies of both nucleon species. 
    !    Bogoliubov   : Bogoliubov transformation, i.e. eigenvectors of the 
    !                   HFB hamiltonian. Ordered according to the eigenvalues
    !                   = the qp energies. Note, that this is not necessarily 
    !                   ordered in terms of "selected" U&V columns!
    !    rho_pairing, : Normal and anomalous density matrices constructed
    !    kappa_pairing| from the Bogoliubov transformation and configmatrix.
    !    configmatrix : Matrix containing the quasiparticle occupations, what
    !                   I call the qp configuration collectively.
    !    qpenergies   : quasiparticle energies corresponding to the columns
    !                   of the Bogoliubov transformation.
    !    blocked_qps  : Indices of the quasiparticles that have been selected
    !                   to be blocked.
    !    partner_qps  : Indices of the partner qps to the blocked qps, i.e. the 
    !                   ones that are closest to being their time-reversal 
    !                   partners.
    !    qp_overlaps  : overlaps between blocked qps and their partners
    !    ifail        : Signals the appearance of problems. If 0, no problem
    !                   has been encountered. (This problem-signalling is not
    !                   entirely operational yet.)
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
    character(len=2), intent(in), allocatable :: BlockLowest(:)
    real(KIND=dp), allocatable, intent(out)   :: qp_overlaps(:)

    integer, allocatable         :: neutron_block(:), proton_block(:)
    integer, allocatable         :: p_blocked(:), n_blocked(:)
    integer, allocatable         :: blocked_qps(:), partner_qps(:)
    integer, allocatable         :: n_partners(:), p_partners(:)
    real(KIND=dp), allocatable   :: p_overlaps(:), n_overlaps(:)
    
    integer                     :: si, sb, N, N2, B, it, i, np, nn
    integer                     :: n_ind, p_ind, NB
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

    !---------------------------------------------------------------------------
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

    case(5,6)
      ! We simply pass which isospin needs to be compared to the modelspwf 
      allocate(proton_block(5))  ; proton_block  = 0 
      allocate(neutron_block(5)) ; neutron_block = 0
      
      if(modelblock .gt. 4) then
          proton_block(modelblock-4)  = 1
      else
          neutron_block(modelblock)   = 1
      endif
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
      &                           gaps(si+1:si+N+N2,si+1:si+N+N2), N, N2)

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
    &              blocktype, neutron_block, n_blocked, n_partners, n_overlaps,&
    &              ifail)   

    call FindFermi(HFBHamil(2*nwn+1:2*nwt,2*nwn+1:2*nwt), HFBlocks(5:8),       &
    &              protons, configmatrix(2*nwn+1:2*nwt),                       &
    &              Bogoliubov(2*nwn+1:2*nwt, 2*nwn+1:2*nwt),                   &
    &              qpenergies(2*nwn+1:2*nwt),Fermi(2), maxhfbiter,             &
    &              blocktype, proton_block, p_blocked, p_partners, p_overlaps, &
    &              ifail)     

    ! 
    
    if(allocated(blocked_qps)) deallocate(blocked_qps)
    if(allocated(partner_qps)) deallocate(partner_qps)
    
    NB = 0 ;  NN = 0 ; NP = 0
    if(allocated(n_blocked)) NN = size(n_blocked)
    if(allocated(p_blocked)) NP = size(p_blocked)
    NB = NP + NN

    if(NB.ne.0) then
      allocate(blocked_qps(NB), partner_qps(NB), qp_overlaps(NB))
      if(allocated(n_blocked)) then
         blocked_qps(   1:NN) = n_blocked
         partner_qps(   1:NN) = n_partners
         qp_overlaps(   1:NN) = n_overlaps
      endif
      if(allocated(p_blocked)) then
         ! We need to offset stuff by the number of neutron qps
         do i=1, NP
           blocked_qps(NN+i) = p_blocked(i)  + sum(HFBlocks(1:4))
           partner_qps(NN+i) = p_partners(i) + sum(HFBlocks(1:4))
           qp_overlaps(NN+i) = p_overlaps(i)
         enddo
      endif
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

  subroutine solvepairing_HFB_gradient( sphamil, gaps, fermi,  Bogo,           & 
  &                          rho_pairing, kappa_pairing, configmatrix,         & 
  &                          qpenergies, BlockType, Blockindices,              &
  &                          blocklowest, blocked_qps, partner_qps,            & 
  &                          p_overlaps, move, maxhfbiter,  ifail) 
    !---------------------------------------------------------------------------
    ! Driver routine for solving the HFB equations by heavy-ball evolution in 
    ! the manifold of Bogoliubov states connected by a Thouless transformation.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! 
    ! Input:
    !   sphamil      : single-particle hamiltonian
    !   gaps         : pairing gaps Delta
    !   fermi        : current guess for the Fermi energy of both isospins
    !   lambda       : Lagrange multiplier for the constraint on the particle
    !                  number dispersion
    !   Bogo         : current Bogoliubov transformation
    !   configmatrix : configuration matrix of the current Bogoliubov vacuum
    !   blocktype    : -|-> Input on blocking that is not used to actually solve
    !   blocklowest  : -|   the HFB problem in this routine. It is only used 
    !   blockindices : -|   to attempt to reconstruct the blocked_qp indices
    !                       after the heavy-ball step.
    !   move         : logical. 
    !                   .true. : update the Bogoliubov transformation
    !                   .false.: don't move (useful for initialisation)                    
    !
    ! Output         :
    !   Bogo         : evolved Bogoliubov transformation
    !   configmatrix : new configuration (see remark below)
    !   fermi        : new Fermi energy
    !   rho_pairing  : density matrix in the current basis
    !   kappa_pairing: anomalous density matrix in the current basis
    !   qpenergies   : quasi-particle energies (or at least an estimate)
    !   blocked_qps  : indices of blocked quasiparticles
    !                  (which are necessary to calculate rotational correction)
    !   partner_qps  : and their partners under time-reversal 
    !                  (well, the closest ones we can find)
    !   p_overlaps   : overlaps between time-reversed blocked qp and the partner
    !   ifail        : 0 = no problem, 1 = something went wrong
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    ! TODO: update these comments, they no longer reflect reality
    !
    ! Note that the gradient solver is conceptually different from the direct 
    ! HFB solver: it requires an initial Bogoliubov transformation to get 
    ! started and requires the knowledge of the "selected" quasiparticles. 
    ! For this reason 
    !
    ! (i) this routine starts by counting "effective" block sizes and 
    !     reorganising the Bogoliubov transformation, such that the "selected"
    !     quasiparticles are on the rhs
    !
    !          B = ( V^*  U )       (1)
    !              ( U^*  V )
    !
    !     with the configmatrix
    !         
    !          C = ( 0   0)         (2)
    !              ( 0   1)
    !
    !     This requires a reorganisation of the quasiparticle columns when the
    !     input corresponds to a blocked configuration for the direct HFB 
    !     solver. In particular, this reorganisation seems to break the 
    !     signature symmetry of the U and V matrices, as because of the blocking
    !     we exchange one (U_k,V_k)^T column for one of the opposite signature.
    !     If we start out with N (=N2) columns of positive (negative) signature,
    !     then on the r.hs. of (1) above, we seem to have N+1, N-1 (or N-1,N+1)
    !     columns after reorganisation. This of course is but an appearance, 
    !     at the end of the day we still have 2*N = 2*N2 columns of every 
    !     type in the complete matrix.
    ! 
    ! (ii) the routine ends by reformatting the Bogoliubov transformation and
    !      quasiparticle energies for the rest of the program. The output will
    !      always in the form of (1) and (2), and hence does not necessarily 
    !      respect the symmetry structure associated with signature anymore.
    !
    !---------------------------------------------------------------------------
    use HFB_gradient

    real(KIND=dp), intent(inout) :: Fermi(2)          , Bogo(:,:)
    real(KIND=dp), intent(inout) :: kappa_pairing(:,:), rho_pairing(:,:)
    real(KIND=dp), intent(inout) :: configmatrix(:)   , qpenergies(:) 
    real(KIND=dp), intent(in)    :: sphamil(:,:),gaps(:,:)
    integer, intent(inout)       :: ifail
    integer, intent(in)          :: maxhfbiter
    logical, intent(in)          :: move

    integer, intent(in)                       :: BlockType
    integer, intent(in), allocatable          :: Blockindices(:)
    character(len=2), intent(in), allocatable :: BlockLowest(:)

    real(KIND=dp)                :: minqp, maxqp, condi, trash
    real(KIND=dp), allocatable   :: tempEqp(:), full_eqp(:), occ(:)
    integer, allocatable         :: blocked_qps(:), partner_qps(:)
    real(KIND=dp), allocatable   :: p_overlaps(:) 


    integer :: si,sb, B, N, N2, T,i, j, stind,endind !,NB ,X(1),Y(1)
  
    ! Statement to stop the compiler complaining about this dummy variable
    if(allocated(blockindices)) trash = 0.0d0
  
    if(.not.allocated(rho_history)) then
      allocate(rho_history(nwt,nwt))            ; rho_history   = 0.0
      allocate(kappa_history(nwt,nwt))          ; kappa_history = 0.0
      allocate(configmatrix_history(2*nwt))     ; configmatrix_history = 0.0
      allocate(Bogoliubov_history(2*nwt, 2*nwt)); Bogoliubov_history = 0.0
    endif

    if(.not.allocated(Z_updates)) then
      allocate(Z_updates(nwt,nwt,2)) ; Z_updates = 0.0d0
    endif
    allocate(tempEqp(nwt))  ; tempEqp = 0 

    ifail = 0
    !---------------------------------------------------------------------------
    ! Saving the history
    Bogoliubov_history   = Bogo
    rho_history          = rho_pairing
    kappa_history        = kappa_pairing 
    
    ! Guess a new Fermi energy if none is there
    if(all(Fermi.eq.0.0))   Fermi = -5
    ! Reorganise the Bogoliubov transformation if needed
    call reorganise_Bogo_gradient(Bogo, configmatrix, grad_blocks)

    if(estimategradparams) then     
      ! Obtain an estimate for the quasi-particle energies to estimate the 
      ! evolution parameters
      full_eqp      = obtain_eqp(sphamil, gaps, Fermi, HFblocks)
      
      minqp = +100000
      maxqp = -100000
      
      do i=1,2*nwt
       if(i.gt.2*nwn) then
        stind = 2*nwn+1
        endind= 2*nwt
       else
        stind = 1
        endind= 2*nwn
       endif
       do j=stind,endind
        ! Consider only positive quasiparticle energies
        if(full_eqp(i) .gt. 0.0d0 .and. full_eqp(j) .gt. 0.0d0) then
          if(full_eqp(i) +full_eqp(j) .lt. minqp) then
            minqp = full_eqp(i) + full_eqp(j) 
          endif        
          if(full_Eqp(i) + full_Eqp(j) .gt. maxqp) then
             maxqp = full_Eqp(i) + full_Eqp(j)
          endif
        endif
       enddo
     
$TR    if(2*full_Eqp(i) .gt. maxqp .and. full_eqp(i) .gt. 0.0d0) then
$TR      maxqp = 2*full_Eqp(i)
$TR    endif
     
$TR    if(2*full_Eqp(i) .lt. minqp .and. full_eqp(i) .gt. 0.0d0) then
$TR      minqp = 2*full_Eqp(i)
$TR    endif
      enddo

      minqp =  max(minqp, gradient_safety)
      condi =  maxqp/minqp
      gradient_mu       = ((sqrt(condi)-1)/(sqrt(condi)+1))**2
      gradient_stepsize =  2.0/maxqp * (  1 + gradient_mu) * 0.9 
      deallocate(full_eqp)
    endif    

    if(move) then
      !-------------------------------------------------------------------------
      ! For clarity, build the occupation factors for the gradient routines
      allocate(occ(nwt)) ; occ = 0.0
    
      si = 0 ; sb = 0
      do B=1,8,2
        N = HFblocks(B)   ; if(N.eq.0) cycle
        N2= HFblocks(B+1)
        T = N + N2
        
        do i=1,T
          occ(si+i) = 1.0d0 - configmatrix(sb+T+i)
        enddo
        si = si +  T
        sb = sb +2*T
      enddo    
    
      !-------------------------------------------------------------------------
      ! Heavy-ball stepping for the neutrons
      call gradient_step(sphamil(1:nwn,1:nwn),gaps(1:nwn,1:nwn),               & 
      &                  neutrons, Bogo(1:2*nwn,1:2*nwn),                      &
      &                  occ(1:nwn),                                           &
      &                  tempEqp(1:nwn),Fermi(1),                              &
      &                  gradient_stepsize, gradient_mu,                       &
      &                  Z_updates(1:nwn,1:nwn,:),                             &
      &                  gradient_precon, HFBgradnorm(1), grad_blocks(1:4),    &
      &                  maxhfbiter, ifail)
      ! and for the protons
      call gradient_step(sphamil(nwn+1:nwt,nwn+1:nwt),gaps(nwn+1:nwt,nwn+1:nwt),& 
      &                  protons,Bogo(2*nwn+1:2*nwt,2*nwn+1:2*nwt),            &
      &                  occ(nwn+1:nwt),                                       &
      &                  tempEqp(nwn+1:nwt),Fermi(2),                          &
      &                  gradient_stepsize, gradient_mu,                       &
      &                  Z_updates(nwn+1:nwt,nwn+1:nwt,:),                     &
      &                  gradient_precon, HFBgradnorm(2), grad_blocks(5:8),    &
      &                  maxhfbiter, ifail)

      deallocate(occ)
    endif
    !---------------------------------------------------------------------------
    ! Copying the Bogoliubov matrix into its 'left side'
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! W.R. 29/07/'22 : this is now handled inside the gradient_step routine.
    !    sb = 0
    !    do B=1,8,2
    !      N = HFblocks(B)   ; if(N.eq.0) cycle
    !      N2= HFblocks(B+1)
    !      T = N + N2
    !  
    !      do i=1,T
    !        ! Populate the columns of the Bogoliubov transformation that have not 
    !        ! been evolved. Note the extra minus sign when time-reversal is conserved.
    !$TR     Bogo(sb  +1:sb  +T, sb+T+1-i) =-Bogo(sb+T+1:sb+2*T, sb+T+i)   
    !$NTR    Bogo(sb  +1:sb  +T, sb+T+1-i) = Bogo(sb+T+1:sb+2*T, sb+T+i)
    !        Bogo(sb+T+1:sb+2*T, sb+T+1-i) = Bogo(sb  +1:sb+  T, sb+T+i)   
    !      enddo 
    !      sb = sb + 2*T
    !    enddo

    !---------------------------------------------------------------------------
    ! Calculate the density and anomalous density matrix
    call PairingMatrices(configmatrix, bogo, rho_pairing, kappa_pairing)
    !---------------------------------------------------------------------------
    ! Final organisation of the Bogoliubov transformation B and QP energies. 
    call correct_ordering_eqp(sphamil,gaps,Fermi,bogo,HFblocks, &
    &                                                  qpenergies,qpdispersions)

    if(blocktype .ne. 4) then
      call figure_out_blocking_structure(sphamil, gaps, Fermi, bogo, &
      &              blocked_qps,partner_qps, p_overlaps,blocktype, blocklowest)
    else
      call figure_out_blocking_structure_EFA( &
      & configmatrix,blocked_qps,partner_qps, p_overlaps,blocktype, blocklowest)
    endif
    !---------------------------------------------------------------------------
    ! Calculate the number dispersion
    HFBdispersion = calc_dispersion_HFB(rho_pairing, kappa_pairing)

    deallocate(tempEqp)
  end subroutine solvepairing_HFB_gradient
  
  subroutine figure_out_blocking_structure(sphamil , gaps, lambda,             &
  &                                             Bogo_ref, bl_qps, part_qps,    &
  &                                             p_overlaps,                    &
  &                                             BlockType, blocklowest) 
    !---------------------------------------------------------------------------
    !  This subroutine attemps to figure out which quasiparticles are blocked
    !  in a given "ordered" Bogoliubov transformation matrix. In addition, it
    !  figures out what the "partner qps" are, the qps that are closest to 
    !  being the time-reversal partner of the blocked qps. 
    !
    !
    ! With ordered Bogoliubov matrix, I mean that it is built as 
    !
    !   Bogo_ref = ( V^T U )
    !              ( U^T V )
    !
    !  with a generalized density matrix in qp representation that is
    ! 
    !              ( 0 0 )
    !              ( 0 1 )
    !
    ! In practice, this means that the "selected" columns of the HFB matrix
    ! are all in the right half. 
    !
    ! In order to determine the blocked qps, we explicitly diagonalize the 
    ! the HFB Hamiltonian by itself, telling the code to build NO excitations.
    ! So
    !
    !   Bogo = ( V'^T U' )
    !          ( U'^T V' )
    !
    ! We then calculate the product of both Bogoliubov transformations, 
    !
    !                  M = Bogo^T Bogo_ref 
    !
    ! If Bogo_ref is close enough to diagonalising the HFB Hamiltonian, then
    ! this matrix is some permutation of the identity matrix, i.e. in every
    ! column there is (to good approximation) one matrix element that is close
    ! to 1. 
    !
    !       M = ( signature= +    cross    )
    !           (    cross       signature=-)
    !
    ! For all non-blocked qps, the 1's will be in the diagonal blocks, i.e. they
    ! will not mix signatures. All blocked qps however, will find their '1's in
    ! the cross parts of this matrix.
    !
    !---------------------------------------------------------------------------
    ! Input:
    ! ------
    !     sphamil     : matrix of the single-particle hamiltonian
    !     gaps        : matrix of the pairing gaps
    !     lambda      : fermi energies of both nucleon species
    !     Bogo_ref    : Bogoliubov transformation to be investigated. 
    !     Blocktype   : Type of blocking to be performed. 
    !                   Acceptable in this routine = 0,2,4
    !     Blocklowest : Set of types of excitation to build.
    !
    ! Output: 
    ! ------
    !     bl_qp       : indices of the blocked quasiparticles, determined by
    !                   the procedure discussed above.
    !     part_qps    : indices of the partner quasiparticles
    !     p_overlaps  : overlap between the time-reversed blocked qp and the
    !                   (detected) partner
    !---------------------------------------------------------------------------
  
    real(KIND=dp), intent(in) :: sphamil(:,:), gaps(:,:), lambda(2)
    real(KIND=dp), intent(in) :: bogo_ref(:,:)
    integer, intent(in)                            :: BlockType
    character(len=2), intent(in), allocatable      :: BlockLowest(:)
    integer, allocatable, intent(out)       :: bl_qps(:), part_qps(:)
    real(KIND=dp), allocatable, intent(out) :: p_overlaps(:)

    real(KIND=dp)             :: part, lambda_copy(2), maxov, tr_over
    real(KIND=dp), allocatable:: HFBHamil(:,:), config(:), bogo(:,:), eqp(:)
    real(KIND=dp), allocatable:: overlap(:,:), tr_qp(:), qpover(:)
    integer                   :: si, sb, N, N2, ifail, B, it, i, j, k, NB, ind
    integer                   :: column
    integer, allocatable      :: blocked_qp(:),blockblock(:), pqp(:)
    logical                   :: check
    
    if(blocktype.ne.2 .and. blocktype.ne.0) then
      print *, 'The blocking identification for the gradient solver is not '
      print *, 'yet capable of dealing with blocktype != 0,2.'
      stop
    endif
    if(.not.allocated(blocklowest)) return
    
    NB = size(blocklowest)
    
    allocate(bl_qps(NB))    ; bl_qps     = 0
    allocate(part_qps(NB))  ; part_qps   = 0
    allocate(p_overlaps(NB)); p_overlaps = 0.0d0
    
    allocate(HFBHamil(2*nwt, 2*nwt)) ; HFBHamil = 0.0d0
    allocate(bogo(2*nwt, 2*nwt))     ; bogo     = 0.0d0
    allocate(config(2*nwt))          ; config   = 0.0d0
    allocate(eqp(2*nwt))             ; eqp      = 0.0d0
    
    !---------------------------------------------------------------------------
    ! Part 1: determine the indices of the blocked quasiparticles
    !---------------------------------------------------------------------------
    
    ! Build the full HFB-hamiltonian
    si      = 0 ; sb = 0
    do B=1,8,2
      N  = HFBlocks(B)    ! Size of the first partner block
      N2 = HFBlocks(B+1)  ! Size of the second partner block
      
      it = 1 ; if (B .gt. 4) it = 2

      HFBHamil(sb+1:sb+2*N+2*N2, sb+1:sb+2*N+2*N2) = ConstructHFBHamil(        &
      &                           sphamil(si+1:si+N+N2,si+1:si+N+N2),          &
      &                           gaps(si+1:si+N+N2,si+1:si+N+N2), N, N2)

      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
    enddo
    !---------------------------------------------------------------------------
    ! now diagonalize the HFB Hamiltonian "as is", without 
    !  (1) any blocking 
    !  (2) any changing of the Fermi energy, i.e. the particle number will not 
    !      be correct.
    
    lambda_copy = lambda
    part = Diagbyblock(HFBHamil(1:2*nwn,1:2*nwn), HFblocks(1:4),               &
    &                  config(1:2*nwn),                                        &
    &                  Bogo(1:2*nwn,1:2*nwn),Eqp(1:2*nwn),      &
    &                  lambda_copy(1), 0 , (/0/), blocked_qp, pqp, qpover,ifail)
    
    part = Diagbyblock(HFBHamil(2*nwn+1:2*nwt,2*nwn+1:2*nwt), HFblocks(5:8),   &
    &                  config(2*nwn+1:2*nwt),                                  &
    &                  Bogo(2*nwn+1:2*nwt,2*nwn+1:2*nwt), Eqp(2*nwn+1:2*nwt),  &
    &                  lambda_copy(2), 0 , (/0/), blocked_qp, pqp, qpover,ifail)

    ! Don't forget to correct the structure of the matrices
    call reorganise_matrices(Bogo,Eqp, config)
    !---------------------------------------------------------------------------
    ! We have now in memory the "reference" Bogoliubov state, which has wrong
    ! particle number, but which should have the "unblocked" U and V columns.
    ! (Note that it is not necesarily the lowest energy HFB vacuum, but rather
    ! the lowest even-even vacuum constructed here...)

    ! We multiply the reference transformation with the transformation in 
    ! memory, to determine which qps have been blocked
    overlap = matmul(transpose(Bogo), Bogo_ref)

    allocate(blockblock(NB)); blockblock = 0
    do i=1,NB
      select case(blocklowest(i))
      case('n+')
        blockblock(i) = 1
      case('n-')
        blockblock(i) = 3
      case('p+')
        blockblock(i) = 5
      case('p-')
        blockblock(i) = 7
      case('n0')
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! If parity is conserved, the code does not know how to deal with 
        ! the possibility of the blocking option being in either of both 
        ! parity blocks. 
$PCONSERVED        print *, 'The blocking identification for the gradient solver is not '
$PCONSERVED        print *, 'yet capable of dealing with "n0", "p0" blocking options.'
$PCONSERVED        stop

        ! If parity is broken, then there is only one possible block for 
        ! the neutron qp excitation to be in 
$PBROKEN blockblock(i) = 1
      case('p0')
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! If parity is conserved, the code does not know how to deal with 
        ! the possibility of the blocking option being in either of both 
        ! parity blocks. 
$PCONSERVED        print *, 'The blocking identification for the gradient solver is not '
$PCONSERVED        print *, 'yet capable of dealing with "n0", "p0" blocking options.'
$PCONSERVED        stop
      
        ! If parity is broken, then there is only one possible block for 
        ! the proton qp excitation to be in 
$PBROKEN blockblock(i) = 5
      end select
    enddo
  
    si      = 0 ; sb = 0 ; ind = 1
    do B=1,8,2
      N  = HFBlocks(B)    ! Size of the first partner block
      N2 = HFBlocks(B+1)  ! Size of the second partner block
      
      it = 1 ; if (B .gt. 4) it = 2

      check = .false. 
      do k=1,NB
        if(blockblock(k) .eq. B) then 
          check = .true.
        endif
      enddo
      
      if(check) then
        ! We check the off-diagonal components of the overlap matrix to find 
        ! the one qp that is most like the non-selected part
        maxov = -100000.0d0
        do i=1,N+N2
          do j=1,N+N2
            if(abs(overlap(sb+i, sb+N+N2+j)) .gt. maxov) then
              maxov = abs(overlap(sb+i, sb+N+N2+j))
              bl_qps(ind) = si + j
            endif
          enddo
        enddo
        ind = ind + 1
      endif      
      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
    enddo
    
    !---------------------------------------------------------------------------
    ! Part 2: determine the closest thing to time-reversal partners that 
    !         we can find. 
    !---------------------------------------------------------------------------
    do k=1,NB
      si = 0 ; sb = 0
      do B=1,8,2
        N = HFBlocks(B)  ; if(N.eq.0) cycle
        N2= HFBlocks(B+1)
        if( bl_qps(k) .lt. si+N+N2 ) exit
      
        si = si +  N+  N2
        sb = sb +2*N+2*N2
      enddo     

      si =   sum(HFBlocks(1:B-1))
      sb = 2*sum(HFBlocks(1:B-1))

      N  = HFBlocks(B)
      N2 = HFBlocks(B+1)
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! At this point:  
      !   bl_qps(k) :  index of the blocked qp
      !
      ! which has a "blocked" structure, i.e. in a signature conserving  
      ! calculation, looks like
      !       
      !  ( 0    )
      !  ( U-_b )
      !  ( V-_b )
      !  ( 0    )
      !
      ! which can be seen as the result of switching U <-> V of a (non-blocked)
      ! qp of the form
      !
      !  ( U+  )    ( V-_b  )
      !  ( 0   ) =  (  0    )
      !  ( 0   )    (  0    )
      !  ( V+  )    ( U-_b  )
      !
      ! Note that we look for the time-reversal partner of the unblocked qp, 
      ! i.e. we look for the quasiparticle that is "close to"
      !
      ! (  0   )   (  0    )
      ! (  U+  ) = (  V-_b )
      ! ( -V+  )   ( -U-_b )
      ! (  0   )   (  0    )
      !
      ! This notation assumes a signature-conserved calculation and a 
      ! blocked quasiparticle of signature +i. The signature = -i blocked
      ! qp is similar, but NOT IMPLEMENTED YET!
      !
      allocate(tr_qp(2*N+2*N2))
      
      ! The right column in the Bogoliubov matrix is the following one.
      ! We are doing some gymnastics, as the qps are indexed with the 
      ! single-particle dimension, i.e. we count only the selected ones.
      ! Meanwhile, the whole Bogoliubov matrix has twice that many columns.
      column =  sb+bl_qps(k)-si+N+N2
      tr_qp(       1:  N     ) = + 0
      tr_qp(  N   +1:  N+  N2) = + Bogo_ref(sb+2*N   +1:sb+2*N+N2,column) 
      tr_qp(  N+N2+1:2*N+  N2) = - Bogo_ref(sb+  N   +1:sb+  N+N2,column) 
      tr_qp(2*N+N2+1:2*N+2*N2) = - 0 

      p_overlaps(k) = -1000000      
      do i=1,N+N2
        tr_over = 0
        do j=1,2*N+2*N2
          tr_over = tr_over + tr_qp(j) * Bogo_ref(sb+j, sb+N+N2+i)
        enddo      
        
        if(abs(tr_over).gt.p_overlaps(k)) then
          p_overlaps(k) = abs(tr_over)
          part_qps(k) = si + i
        endif
      enddo
      
      deallocate(tr_qp)
    enddo
    
  end subroutine figure_out_blocking_structure
  
  subroutine figure_out_blocking_structure_EFA(&
  &          config,bl_qps, part_qps, p_overlaps, BlockType, blocklowest) 
    !---------------------------------------------------------------------------
    ! Identify the blocked quasiparticle(s) and their (almost) time-reversal
    ! partner(s) in the case of Equal Filling blocking.
    !
    ! Input:
    ! ------
    ! config     : configuration (1-f) factors
    ! blocktype  : type of blocking to be done
    ! blocklowest: types of qps to be blocked
    ! 
    ! Output:
    ! -------
    ! bl_qps     : identified blocked quasiparticles
    ! part_qps   : identified partner quasiparticles
    ! p_overlaps : overlaps with time-reversal, i.e. <B|T|P>
    ! 
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)                      :: config(:)
    integer, allocatable, intent(out)              :: bl_qps(:), part_qps(:)
    real(KIND=dp), allocatable, intent(out)        :: p_overlaps(:)
    integer, intent(in)                            :: BlockType
    character(len=2), intent(in), allocatable      :: BlockLowest(:)

    integer :: i, NB, si, sb, N, N2, T, ind, B

    NB = size(blocklowest)

    if(blocktype.ne.4) then
      print *, 'figure_out_blocking_structure_EFA should only be called for'
      print *, 'blocktype.eq.40'
      stop
    endif

    allocate(bl_qps(NB))    ; bl_qps     = 0
    allocate(part_qps(NB))  ; part_qps   = 0
    allocate(p_overlaps(NB)); p_overlaps = 0.0d0
    
    si      = 0 ; sb = 0 ; ind = 0
    do B=1,8,2
      N  = HFBlocks(B)    ! Size of the first partner block
      N2 = HFBlocks(B+1)  ! Size of the second partner block
      T  = N + N2
      
      do i=1,T
        if(config(sb+i) .eq. 0.5d0) then
          ind             = ind + 1
          bl_qps(ind)     = si + i
          part_qps(ind)   = si + i
          p_overlaps(ind) = 1.0d0
        endif
      enddo
      
      si = si +   T
      sb = sb + 2*T
    enddo

  end subroutine figure_out_blocking_structure_EFA
  
  function obtain_eqp(sphamil, gaps, lambda, blocks) result(eigen)
    !---------------------------------------------------------------------------
    ! Obtain the quasiparticle energies by constructing and diagonalizing the 
    ! HFB hamiltonian. 
    !
    ! Input :
    !  sphamil : matrix of the single-particle hamiltonian
    !  gaps    : matrix of the pairing gaps
    !  lambda  : Fermi energies of both nucleon species
    !  blocks  : size of the symmetry blocks of the HFB Hamiltonian
    !
    ! Output:
    !  eigen   : the eigenvalues = qp energies of the HFB Hamiltonian 
    ! 
    !---------------------------------------------------------------------------
    ! Note 
    ! (1) this routine offers nothing beyond that, and does not offer any 
    !     other side-effects. For a more complete routine, see diagbyblock 
    !     in the HFB_direct module.
    ! (2) this routine produces the ACTUAL qp energies of a given HFB 
    !     Hamiltonian, as opposed to the subroutine correct_ordering_eqp
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in)   :: sphamil(:,:), gaps(:,:), lambda(2)
    real(KIND=dp), allocatable  :: HFBhamil(:,:), work(:), A(:,:)
    integer, intent(in)         :: blocks(8)
    real(KIND=dp),allocatable   :: eigen(:)

    integer       :: si, sb, B, N, N2, lwork, i, ifail, it
  
    si      = 0 ; sb = 0
    allocate(eigen(2*sum(blocks))) ;   eigen   = 0
    
    allocate(HFBHamil(2*sum(blocks), 2*sum(blocks)))
    do B=1,8,2
      N  = Blocks(B)    ! Size of the first partner block
      N2 = Blocks(B+1)  ! Size of the second partner block
        
      HFBHamil(sb+1:sb+2*N+2*N2, sb+1:sb+2*N+2*N2) = ConstructHFBHamil(        &
      &                           sphamil(si+1:si+N+N2,si+1:si+N+N2),          &
      &                           gaps(si+1:si+N+N2,si+1:si+N+N2), N, N2)

      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
   enddo
   
   si=0 ; sb =0
   do B=1,8
   
      N  = Blocks(B) ; if(N.eq.0) cycle
      it = 1 ; if(B.gt.4) it = 2

      A = HFBHamil(sb+1:sb+2*N, sb+1:sb+2*N)
      do i=1,N
        A(i  ,i  ) = A(i  , i  ) - lambda(it)
        A(i+N,i+N) = A(i+N, i+N) + lambda(it)
      enddo
                        
      ! Diagonalize
      lwork = -1; allocate(work(1))
      call DSYEV( 'V', 'U', 2*N, A, 2*N, eigen(sb+1:sb+2*N),work,lwork,ifail)
      lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
      call DSYEV( 'V', 'U', 2*N, A, 2*N, eigen(sb+1:sb+2*N),work,lwork,ifail)
      deallocate(work)

      si = si +   N 
      sb = sb + 2*N 
   enddo
  
  end function obtain_eqp
  
  subroutine correct_ordering_eqp(sphamil,gaps,lambda,bogo,blocks, eigen,disp) 
    !---------------------------------------------------------------------------
    ! This routine obtains estimates of the qp-energies of the HFB Hamiltonian H
    ! by taking the diagonal matrix elements of
    !
    !       W^T H W
    !
    ! where W is the Bogoliubov transformation. This ensures that the estimate 
    ! of an qp energy is correctly associated with the corresponding 
    ! (approximate) eigenvector of H. We also calculate the associated 
    ! dispersions.
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    !
    ! Input : 
    ! -------
    !   sphamil : matrix of the single-particle hamiltonian
    !   gaps    : matrix of the pairing gaps
    !   lambda  : Fermi energies for both nucleon species
    !   blocks  : sizes of the symmetry blocks of the HFB Hamiltonian
    !
    ! Output: 
    ! -------
    !   eigen   : the diagonal matrix elements of W^T H W
    !   disp    : the dispersion of the vectors in W
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)   :: sphamil(:,:), gaps(:,:), lambda(2), bogo(:,:)
    real(KIND=dp), allocatable  :: HFBhamil(:,:)
    integer, intent(in)         :: blocks(8)
    real(KIND=dp),  intent(out) :: eigen(:), disp(:)
    
    real(KIND=dp), allocatable  :: prod(:,:)
    integer       :: si, sb, B, N, N2, T, i,  it
  
    si    = 0 ; sb   = 0
    eigen = 0 ; disp = 0
    
    allocate(HFBHamil(2*sum(blocks), 2*sum(blocks)))
    !  ^
    !  |
    ! Actually inefficient memory useage: we can allocate this thing inside
    ! the loop with the size of individual blocks....

    do B=1,8,2
      N  = Blocks(B)    ! Size of the first partner block
      N2 = Blocks(B+1)  ! Size of the second partner block
      T  = N + N2
      
      it = 1 ; if(B.gt.4) it = 2
      
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! I construct the HFB Hamiltonian by hand as it needs to have the right 
      ! structure from the start
      HFBHamil(sb  +1:sb+  T,sb  +1:sb+  T) = +sphamil(si+1:si+T,si+1:si+T)
      HFBHamil(sb+T+1:sb+2*T,sb+T+1:sb+2*T) = -sphamil(si+1:si+T,si+1:si+T)
  
$NTR  HFBHamil(sb+T+1:sb+2*T,sb  +1:sb  +T) = -gaps(si+1:si+T,si+1:si+T)
$TR   HFBHamil(sb+T+1:sb+2*T,sb  +1:sb  +T) = +gaps(si+1:si+T,si+1:si+T)
       
      HFBHamil(sb+1  :sb+  T,sb+T+1:sb+2*T) = +gaps(si+1:si+T,si+1:si+T)
      
      do i=1,T
        HFBHamil(sb  +i,sb  +i) =  HFBHamil(sb  +i,sb  +i) -lambda(it)
        HFBHamil(sb+T+i,sb+T+i) =  HFBHamil(sb+T+i,sb+T+i) +lambda(it)
      enddo

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Multiply with the Bogliubov transformation from the right      
      prod = matmul(HFBHamil(sb+1:sb+2*T,sb+1:sb+2*T), &
      &             Bogo(sb+1:sb+2*T,sb+1:sb+2*T))
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! ..... and with its transpose from the left
      prod = matmul(transpose(Bogo(sb+1:sb+2*T,sb+1:sb+2*T)), prod)
      ! We don't actually need ALL these matrix elements, but it is a nice 
      ! debugging tool to just be able to 'print HFBHamiltonian'....
      
      do i=1, 2*T
        eigen(sb+i) = prod(i,i)
      enddo

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Squaring the HFB Hamiltonian
      HFBHamil(sb+1:sb+2*T,sb+1:sb+2*T) = matmul( &
      &    HFBHamil(sb+1:sb+2*T,sb+1:sb+2*T),HFBHamil(sb+1:sb+2*T,sb+1:sb+2*T))     

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Multiply with the Bogliubov transformation from the right      
      prod = matmul(HFBHamil(sb+1:sb+2*T,sb+1:sb+2*T), &
      &             Bogo(sb+1:sb+2*T,sb+1:sb+2*T))
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! ..... and with its transpose from the left
      prod = matmul(transpose(Bogo(sb+1:sb+2*T,sb+1:sb+2*T)), prod)
      
      do i=1, 2*T
        ! Dispersion = <H^2> - E^2
        disp(sb+i) = prod(i,i) - eigen(sb+i)**2
      enddo

      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
   enddo
  
  end subroutine correct_ordering_eqp

  pure subroutine reorganise_Bogo_gradient(Bogo, config, effblocks) 
    !---------------------------------------------------------------------------
    ! The gradient solver fundamentally only works with "ordered" Bogoliubov
    ! transformation as input, i.e. with the "selected" U and V matrices on the
    ! right-most half of the Bogoliubov transformation. This is equivalent  
    ! to having a Bogoliubov transformation with a trivial configuration matrix.
    !
    ! This routine does exactly that: it takes as input a transformation with
    ! a corresponding configuration matrix and orders them. In the meantime, 
    ! it counts the effective block sizes of the matrices to be passed in to the
    ! gradient solver.
    !
    ! Input: 
    !    Bogo      :  non-ordered Bogoliubov transformation
    !    config    :  configuration matrix of the vacuum
    !    effblocks :  If non-zero, the routine will not reorder things.
    !
    ! Output:
    !    Bogo      : ordered Bogoliubov transformation
    !    config    : trivial configuration matrix, i.e. (0,0,...,0,1,1,....,1)
    !                in every symmetry block
    !    effblocks : modified block sizes to pass into the gradient solver, 
    !                which might be larger or smaller than the HFBlocks.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(inout) :: Bogo(:,:)
    real(KIND=dp), intent(inout) :: config(:)
    integer, intent(inout)       :: effblocks(8)
    real(KIND=dp), allocatable   :: tempBogo(:,:), tempconfig(:)
   
    integer                      :: B, N, N2, T, sb,  NP, ind, ind2, i
    integer, allocatable         :: indices(:) 
    real(KIND=dp)                :: compare
    
    allocate(tempBogo(2*nwt, 2*nwt), tempconfig(2*nwt))
    tempBogo = 0.0d0
    !---------------------------------------------------------------------------
    ! if these variables have been set, then we know that the 
    ! organisation of the Bogoliubov transformation is okay
    if(any(effblocks.ne.0)) return

    !---------------------------------------------------------------------------
    ! For time-reversal invariant calculations, no further work is necessary
$TR effblocks = HFblocks    
$TR return
  
    !---------------------------------------------------------------------------
    ! If some qp excitations were made, then we need to move things around
    ! in the Bogoliubov transformation
    effblocks = HFblocks
    sb = 0
    do B=1,8,2
      N = HFBlocks(B)   ; if(N.eq.0) cycle
      N2= HFBlocks(B+1)
      T = N+N2

      ! Counting the number of positive signature states
      Np = 0
      do i=1,N
        if(config(sb+T+i).eq.1.0d0) Np = Np+1
      enddo
      do i=N+1,T
        if(config(sb+T+i).eq.0.0d0) Np = Np+1
      enddo
      ! Getting the indices right
      allocate(indices(T)) ; ind = 0 ; ind2 = 0
      do i=1,N
        if(config(sb+T+i).eq.1.0d0) then
          ind          = ind + 1
          indices(ind) = i
        else
          ind2             = ind2 + 1
          indices(Np+ind2) = i
        endif
      enddo
      do i=N+1,T
        if(config(sb+T+i).eq.1.0d0) then
          ind2             = ind2 + 1
          indices(Np+ind2) = i
        else
          ind              = ind + 1
          indices(ind )    = i
        endif
      enddo

      do i=1,T
        if(config(sb+T+indices(i)).eq.1.0d0) then
          tempBogo(sb+1:sb+2*T ,sb+T+i)  = Bogo(sb+1:sb+2*T,sb+T+indices(i))  
          tempconfig  (sb+T+i)           = config(sb+T+indices(i))
!      config(sb+T+1:sb+2*T) = 1.0d0
        elseif(config(sb+T+indices(i)).eq.0.5d0) then
          tempBogo(sb+1:sb+2*T ,sb+T+i)  = Bogo(sb+1:sb+2*T,sb+T+indices(i))  
          tempconfig  (sb+T+i)           = config(sb+T+indices(i))
        else
          tempBogo(sb+1  :sb+  T,sb+T+i) = Bogo(sb+T+1:sb+2*T,sb+T+indices(i))    
          tempBogo(sb+T+1:sb+2*T,sb+T+i) = Bogo(sb+  1:sb+  T,sb+T+indices(i))    
          tempconfig  (sb+T+i)           = 1-config(sb+T+indices(i))
        endif
      enddo    
      deallocate(indices)

      effBlocks(B)   = Np
      effBlocks(B+1) = T - Np

      sb = sb + 2*T
    enddo
    Bogo = tempbogo

    !---------------------------------------------------------------------------
    ! We have moved the Bogoliubov transformation around, such that we have 
    ! U and V matrices that are "selected" on the r.h.s.
    !
    ! This situation might have also produced itself by reading a .wf file 
    ! that had information from a calculation with a gradient solver. Hence, 
    ! we still need to count the block sizes (again) to be sure we pass the 
    ! gradient solver the right structure of U and V matrices.
    sb = 0
    do B=1,8,2
      N = HFBlocks(B)   ; if(N.eq.0) cycle
      N2= HFBlocks(B+1)
      T = N+N2
      
      Np = 0
      do i=1,N+N2
        !  compare = 1     compare = 0     (to numerical precision)
        !
        !   ^                 ^
        !   |                 | 
        !  
        ! ( U^+ )           ( 0  )
        ! ( 0   )           ( U^-)
        ! ( 0   )           ( V^-)
        ! ( V^+ )           ( 0  )
        !    
        ! 
        compare    = sum(Bogo(sb+1    :sb+N,sb+T+i)**2) &
        &          + sum(Bogo(sb+T+N+1:sb+2*T,sb+T+i)**2)
        if(abs(compare-1).lt.0.0001d0) Np = Np+ 1
      enddo
      
      effBlocks(B)   = Np
      effBlocks(B+1) = T - Np
      sb = sb + 2 * T
    enddo
    
    !---------------------------------------------------------------------------
    ! Finally, we change the configmatrix around.
    sb = 0
    config = 0.0d0
    do B=1,8,2
      N = HFBlocks(B)   ; if(N.eq.0) cycle
      N2= HFBlocks(B+1)
      T = N+N2

      do i=1, T
        config(sb  +i) = 1 - tempconfig(sb+T+i)
        config(sb+T+i) =     tempconfig(sb+T+i)
      enddo
      sb = sb + 2*T
    enddo
    
  end subroutine reorganise_Bogo_gradient

  pure function calc_dispersion_HFB(rho,kappa) result(dispersion)
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
    ! Input: 
    !      rho  : Normal and anomalous density matrix  
    !      kappa|
    ! Output:
    !      dispersion: particle number dispersion for both nucleon species
    !---------------------------------------------------------------------------
    !
    ! Note, that at T = 0, we have that (kappa * kappa^{\dagger}) = rho(1-rho).
    ! So in that case, we have 
    !  < Delta N^2 > = < N^2 > - <N>^2 = 2 * Tr(rho(1-rho))
    ! which is the old formula from EV8, CR8, etc...
    !
    ! We implement however the formula above, since this is the one that 
    ! correctly generalizes to T != 0.
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)    :: rho(:,:), kappa(:,:)
    real(KIND=dp), allocatable   :: chi(:,:)
$NTR real(KIND=dp), allocatable  :: k2(:,:)
    real(KIND=dp)                :: dispersion(2)
    integer                      :: N, B, si, i, it, N2, T

    si            = 0
    dispersion = 0.0
    allocate(chi(nwt, nwt))
    do B=1,8,2 
      N  = HFBlocks(B)  ;  if(N .eq. 0) cycle 
      N2 = HFBlocks(B+1) 
      
      T = N + N2

      it  = 1
      if(B .gt. 4) it = 2
      
      ! rho squared
      chi = matmul(rho(si+1:si+T, si+1:si+T), rho(si+1:si+T, si+1:si+T) ) 
      
      ! kappa * kappa^*
$NTR  k2  = matmul(           kappa(si+1:si+T, si+1:si+T), & 
$NTR      &         transpose(kappa(si+1:si+T, si+1:si+T))) 
         
      !                                       Tr rho - Tr rho^2
      do i=1,T
          dispersion(it) = dispersion(it) + rho(si+i, si+i)           &
          &                                     - chi(i,i)            &
$TR       &                                     + kappa(si+i, si+i)**2
$NTR      &                                     + k2(i,i)
      enddo
      si = si + T
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
            &     config(sb+  k)*bogo(sb+     i,column) * bogo(sb+     j,column)
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

  function ConstructHFBHamil(sphamil, gaps, N, N2) result(H)
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
    !
    ! 13/04/21, W.R.: this option is no longer supported.
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
    ! Output: 
    !       H  : the HFB Hamiltonian in a format that can be easily passed to
    !            diagonalization routines.
    !
    !  -  -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! REMOVED:
    !  r, k    : rho and kappa pairing matrices for the addition of the 
    !            generalized density matrix. Should also be in correct
    !            block structure.
    !   
    !  gauge   : real parameter alpha
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)   :: sphamil(:,:), gaps(:,:)
    real(KIND=dp), allocatable  :: H(:,:)
    integer, intent(in)         :: N, N2
    integer                     :: T

    T = N + N2
    allocate(H(2*T,2*T)) 
    H = 0

    !---------------------------------------------------------------------------
    ! SP hamil
    ! Block 1
    H(  1:N   ,  1:N)    =  sphamil(  1:N   ,   1:N)    
    H(N+1:N+N2,N+1:N+N2) = -sphamil(N+1:N+N2, N+1:N+N2) 

    ! Block 2
    H(N+  N2+1:  N+2*N2,N+  N2+1:  N+2*N2) =  sphamil(N+1:N+N2,N+1:N+N2)
    H(N+2*N2+1:2*N+2*N2,N+2*N2+1:2*N+2*N2) = -sphamil(  1:N   ,   1:N)  

    !---------------------------------------------------------------------------
    ! Gaps
    ! Block 1
    H(  1:N   ,N+1:2*N)  = gaps(  1:N   , N2+1:N2+N) 
    H(N+1:2*N ,  1:N  )  = gaps(  1:N   , N2+1:N2+N) 

    ! Block 2
    H(2*N   +1:2*N+N2  ,2*N+N2+1:2*N+2*N2) =  gaps(N+1:N+N2, 1:N2)
    H(2*N+N2+1:2*N+2*N2,2*N   +1:2*N+  N2) =  gaps(N+1:N+N2, 1:N2) 
    

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

    call start_timer(T_gaps)

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
        
$NTR        do wave2=1,N
$TR         do wave2=wave1,N
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
$TR       HFBgaps(indb,inda) = HFBgaps(inda,indb)
        enddo
      enddo
      !-------------------------------------------------------------------------      
      ! We have now calculated the gaps (without cutoffs) in the basis that 
      ! is currently in storage. This can either be the HF basis or not, but
      ! we need the gaps in the HF-basis to calculate the cutoffs. 
      !-------------------------------------------------------------------------      
      ! I could do this transformation with the transformation routine below, 
      ! but this was  historically first and block-wise, so I keep it. 
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
    
    call stop_timer(T_gaps)

  end subroutine calcHFBgaps

  subroutine calc_gaps_HF(gaps, transfo)
    !---------------------------------------------------------------------------
    ! Input:
    !    gaps   : pairing gaps in the basis in memory
    !    transfo: unitary transformation from the basis in memory to the one 
    !             we want.
    !
    ! Output:
    !    gaps   :  gaps in the basis defined by the transformation.
    !---------------------------------------------------------------------------
    real(KIND=dp), allocatable, intent(inout) :: gaps(:,:) 
    real(KIND=dp), allocatable, intent(in)    :: transfo(:,:)
    real(KIND=dp), allocatable                :: temp(:,:)
    integer :: si, N, B, T
    
    temp = gaps
    gaps = 0.0d0
    
    si = 0
    do B=1,8,2
      N = HFblocks(B)         ; if(N.eq.0) cycle
      T = HFBlocks(B+1) + N
        
      gaps(si+1:si+T, si+1:si+T) = &
      & matmul(transpose(transfo(si+1:si+T,si+1:si+T)),&
      &                                          temp(si+1:si+T,si+1:si+T))
      gaps(si+1:si+T, si+1:si+T) = &
      &    matmul(gaps(si+1:si+T,si+1:si+T),transfo(si+1:si+T,si+1:si+T))

      si = si + T
    enddo
  
  end subroutine calc_gaps_HF 

  subroutine PrintHFBconvergence(rho_pairing, kappa_pairing, Bogo)
    !---------------------------------------------------------------------------
    ! Prints out some convergence info on the HFB subproblem.
    !   a) sqrt(sum( (rho*rho - rho + kapppa * kappa^T)**2)
    !   b) sqrt(sum( (rho*kappa - kappa * rho)**2))
    ! 
    ! These should be (numerically) vanishingly small at convergence, when we
    ! have solved the HFB problem consistently.
    !
    ! Input : 
    !   rho_pairing  : density matrix
    !   kappa_pairing: anomalous density matrix
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Note that these numbers do not vanish when we do either
    !   (i)  finite-temperature calculations (they increase as T increases)
    !   (ii) Equal Filling-style blocking
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in)  :: rho_pairing(:,:),kappa_pairing(:,:),bogo(:,:)
    real(KIND=dp)              :: test1(8), test2(8), test3(8)
    real(KIND=dp), allocatable :: A(:,:) , r(:,:), k(:,:), check(:,:)
    integer                    :: N, B, si, N2, T, sb, i,j

    1 format (' HFB convergence:   (N,+)    (N,-)    (P,+)    (P,-)')
    2 format ('  r^2-r+k*k^T    = ',  4es9.2)
    3 format ('  r*k-k*r        = ',  4es9.2)
    4 format ('  unitarity      = ',  4es9.2)    
!    5 format (' Grad. norm      = ',  4es9.2)

    print *
    print 1
    
    si = 0
    do B=1,8,2
        N = HFBlocks(B) ; if (N.eq. 0) cycle
        N2= HFBlocks(B+1)

        T = N + N2
        allocate(A(T,T), r(T,T), k(T,T))
        r = rho_pairing(si+1:si+T  ,si+1:si+T)
        k = kappa_pairing(si+1:si+T,si+1:si+T)
        
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

    si = 0 ; sb = 0
    test3 = 0.0d0 
    do B=1,8,2
        N = HFBlocks(B) ; if (N.eq. 0) cycle
        N2= HFBlocks(B+1)

        T = N + N2
        check = matmul(transpose(Bogo(sb+1:sb+2*T,sb+1:sb+2*T)), &
        &                        Bogo(sb+1:sb+2*T,sb+1:sb+2*T))
        
        do i=1,2*T
          do j=1,2*T
            if(i.eq.j) then
              test3(B) = test3(B) + (check(i,j) - 1.0d0)**2
            else 
              test3(B) = test3(B) + (check(i,j))**2
            endif
          enddo
        enddo
        sb = sb + 2*T
    enddo
    test3 = sqrt(test3)
    print 4, test3(1),test3(3),test3(5),test3(7) 

  end subroutine PrintHFBconvergence

  subroutine Canonical(rho_pairing, kappa_pairing, rho_can, kappa_can,    &
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
    integer, intent(out)       :: ifail
    real(KIND=dp), allocatable :: tmp(:,:), work(:), temp_occ(:)
    
    integer, allocatable       :: indices(:)
    integer                    :: si, N, N2, B, i,j,  lwork, effN, ind
    integer                    :: ii, jj
    $NTR   real(KIND=dp)       :: mindiff, diff
 
    !---------------------------------------------------------------------------
    ! a) Diagonalize rho
    !
    ! Transforms as rho'  = D^dagger rho D
    !---------------------------------------------------------------------------
    si         = 0
    rhotransfo = 0 ; kappatransfo = 0
    do B=1,8
      N =HFBlocks(B) ;  if(N .eq. 0) cycle 
      
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Search for all the states in the Hartree-Fock basis that don't 
      ! participate in pairing; for a state i this means
      ! 
      !       rho_{ij} = 0 if i != j and rho_ii = 0 or 1.
      !       kappa_ij = 0 for all j      
      effN = 0 ; ind = 1
      allocate(indices(N)) ; indices = 0
      do i=1,N
        if (      any(abs(rho_pairing(i,i+1:N)) .gt. rho_cutoff)  &
        &  .or. (       abs(rho_pairing(i,i)  ).gt.rho_cutoff  &    
        &         .and. abs(rho_pairing(i,i)-1).gt.rho_cutoff))&
        then
            indices(ind) = i
            ind          = ind + 1
        endif        
      enddo
      effN = ind - 1
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Construct the submatrix of rho that we will diagonalize
      allocate(tmp(effN,effN), temp_occ(effN)) 

      do i=1, effN
        ii = indices(i)
        do j=1, effN
          jj = indices(j)
          tmp(i,j) = rho_pairing(si + ii, si + jj)
        enddo
      enddo
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! Diagonalize rho in the limited subspace
      if(effN .ne. 0) then         ! only do this if pairing has not collapsed
        lwork = -1 ; allocate(work(1))
        call DSYEV( 'V', 'U', effN, tmp, effN, temp_occ, work, lwork, ifail)
        lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
        call DSYEV( 'V', 'U', effN, tmp, effN, temp_occ, work, lwork, ifail)
        deallocate(work)
      endif
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! And then reconstruct (1) the occupations in the full space
      !                      (2) the transformation in the full space
      
      ! Make this whole transformation trivial in the HF basis
      ! And use the diagonal elements for rho_can
      rhotransfo(si+1:si+N,si+1:si+N) = 0  
      do i=1, N
        rhotransfo(si+i,si+i) = 1
        rho_can(si+i) = rho_pairing(si+i, si+i)
      enddo
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      ! And then overwrite these arrays with the non-trivial information
      do i=1, effN
        ii = indices(i)
        rho_can(si+ii) = temp_occ(i)
        do j=1, effN
          jj = indices(j)
          
          rhotransfo(si+ii,si+jj) = tmp(i,j)
        enddo
      enddo
      si = si + N
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
      deallocate(tmp, temp_occ, indices)
    enddo
    
    ! Time-reversal
$TR    rho_can = 2*rho_can

    do i=1,nwt
$TR      if(rho_can(i).gt.2.0d0) rho_can(i) = 2.0d0
$NTR     if(rho_can(i).gt.1.0d0) rho_can(i) = 1.0d0
      if(rho_can(i).lt.0.0) rho_can(i) = 0.0d0
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

!-------------------------------------------------------------------------------
!    This way of constructing the conjugate partners is close to that used
!    in Esperance. However, it depends on two numerical cutoffs, and failing 
!    those could result in small differences in the final energy.
!    W.R. 26/07/22
!$NTR  do i=1, N
!$NTR    do j=N+1, N+N2
!$NTR      if(abs(rho_can(si+i) - rho_can(si+j)).lt.1d-6) then
!$NTR        if( abs(tmp(i,j)).gt. 2e-5) then
!$NTR            conjugp(si+i) =  si+j
!$NTR            conjugp(si+j) =  si+i
!$NTR            kappa_can(si+i) = tmp(i,j)
!$NTR            kappa_can(si+j) = tmp(j,i)   
!$NTR        endif
!$NTR      endif    
!$NTR    enddo
!$NTR  enddo
!-------------------------------------------------------------------------------
!    This way of constructing conjugp is somewhat less arbitrary: it scans the
!    opposite signature block for the level with the occupation that is the 
!    closest in value. This might not necessarily be better than the other 
!    way commented above, but it at least does not depend on parameters.
!    W.R. 26/07/22
$NTR     do i=1, N
$NTR       mindiff = 10000
$NTR       conjugp(si+i) = 0
$NTR       do j=1, N2
$NTR          diff = abs(rho_can(si+i) - rho_can(si+N+j))
$NTR          if(diff .lt. mindiff ) then
$NTR            mindiff = diff
$NTR            conjugp(si+i) = si+N+j
$NTR          endif
$NTR       enddo
$NTR       !print *, i, mindiff, conjugp(si+i)
$NTR       conjugp(conjugp(si+i)) = si + i
$NTR       kappa_can(si+i)          = tmp(i,conjugp(si+i)-si)
$NTR       kappa_can(conjugp(si+i)) = tmp(conjugp(si+i)-si,i)
$NTR     enddo
!-------------------------------------------------------------------------------
      deallocate(tmp)
      si = si + N + N2
    enddo    
    
   end subroutine Canonical

   pure function construct_generalized_density(rho, kappa) result(R)
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

    allocate(R(4*N, 4*N)) ; R = 0.0d0

    R(  1:  N,   1:  N) = rho
    R(N+1:2*N, N+1:2*N) = rho
    
    R(2*N+1:3*N, 2*N+1:3*N) =-rho
    R(3*N+1:4*N, 3*N+1:4*N) =-rho

    do i=1, 2*N
        R(2*N+i, 2*N+i) = R(2*N+i, 2*N+i) + 1.0d0
    enddo
    
    R(  1:  N,3*N+1:4*N) = kappa
    R(N+1:2*N,2*N+1:3*N) = -transpose(kappa)

    R(2*N+1:3*N, N+1:2*N) = -kappa
    R(3*N+1:4*N,   1:  N) = transpose(kappa)


!$NTR    R(N+1:2*N,  1:  N) =-kappa
!$TR     R(N+1:2*N,  1:  N) =+kappa

   end function construct_generalized_density

   subroutine clean_HFB
    if(allocated(HFBgaps))  deallocate(HFBGaps)
   end subroutine clean_HFB

   subroutine update_qp_angmom(Bogo)
    !---------------------------------------------------------------------------
    ! Calculate the angular momentum expectation values for the HFB 
    ! quasiparticle operators. 
    !
    ! What we calculate here is 
    !
    !  <  qp=k | J_mu | qp=k >  = sum_i - |V_{i,k}|^2 < i | J_mu | i >
    !                                   + |U_{i,k}|^2 < i | J_mu | i >
    !
    ! where i is a single-particle index and k is a quasi-particle index.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Currently only <J_z>, hardcoded for CR8-like symmetries.
    !
    !---------------------------------------------------------------------------
    integer                   :: B, k, N, N2, si, sb, wave, i
    real(KIND=dp), intent(in) :: Bogo(:,:)
    
    if(.not.allocated(qp_J)) then
      allocate(qp_J(3,2*nwt), qp_JTR(3,2*nwt), qp_JTI(3,2*nwt))
    endif
    qp_J   = 0.0d0 ;  qp_JTR = 0.0d0 ; qp_JTI = 0.0d0
    
    si = 0 ; sb = 0
    do B=1,8,2
      N = HFBlocks(B)   ; if(N.eq.0) cycle
      N2= HFBlocks(B+1)
      do wave=1,2*N+2*N2
        do k=1,3
          qp_J(k,sb+wave)   = 0.0d0
          qp_JTR(k,sb+wave) = 0.0d0
          qp_JTI(k,sb+wave) = 0.0d0
          do i=1,N+N2
            qp_J(k,sb+wave) = qp_J(k,sb+wave) +                                &
            !          U^2                       V^2
            &   (Bogo(sb+i,sb+wave)**2 - Bogo(sb+N+N2+i,sb+wave)**2)           &
            &                                              * spwf_J(k,si+i,si+i)   
          enddo
        enddo
      enddo
      
      si = si +  N +   N2
      sb = sb +2*N + 2*N2
    enddo
    
   end subroutine update_qp_angmom
end module
