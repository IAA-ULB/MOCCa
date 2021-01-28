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
  !---------------------------------------------------------------------------
  ! History of the pairing matrices, for mixing purposes.
  real(KIND=dp), allocatable ::  rho_history(:,:), kappa_history(:,:)
  real(KIND=dp), allocatable ::  configmatrix_history(:) 
  real(KIND=dp), allocatable ::  Bogoliubov_history(:,:)
  real(KIND=dp), allocatable ::  overrho, overkap

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

  subroutine solvepairing_HFB(fermi, Bogoliubov, rho_pairing, kappa_pairing,   &
  &                           configmatrix, qpenergies,HFBmix, HFBmixtype,     &
  &                     BlockType,Blockindices, blocklowest, blocked_qps, ifail)

    !---------------------------------------------------------------------------
    ! Driver routine for the solving of the HFB equations, represented in the 
    ! Hartree-Fock basis.
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(inout) :: Fermi(2)
    real(KIND=dp), intent(inout) :: Bogoliubov(:,:)
    real(KIND=dp), intent(inout) :: kappa_pairing(:,:), rho_pairing(:,:)
    real(KIND=dp), intent(inout) :: configmatrix(:), qpenergies(:)
    
    ! Options for the mixing of the HFB configurations
    real(KIND=dp)                :: HFBmix
    integer, intent(in)          :: HFBmixtype

    ! Configuration for the blocking
    integer, intent(in)          :: Blockindices(:)
    integer, intent(in)          :: BlockType
    integer, intent(out)         :: ifail
    character(len=2), intent(in) :: BlockLowest(:)
    integer, allocatable         :: neutron_block(:), proton_block(:)
    integer, allocatable         :: blocked_qps(:), p_blocked(:), n_blocked(:)

    ! Quantities for the HFB hamiltonian
    real(KIND=dp)              :: sphamil(nwt,nwt), HFBHamil(2*nwt, 2*nwt)
    
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
    ! a) We construct the HFB-hamiltonian for every block. 
    !    We pass in everything to the routine by PAIRS of blocks
    sphamil = 0
    si      = 0 ; sb = 0
    do B=1,8,2
      N  = HFBlocks(B)    ! Size of the first partner block
      N2 = HFBlocks(B+1)  ! Size of the second partner block
      
      it = 1 ; if (B .gt. 4) it = 2

      do wave1=1,N+N2
        sphamil(si+wave1,si+wave1) = spenergies(si+wave1)
      enddo

      if(.not.allocated(HFBgaps)) then
        print *, 'HFB gaps are not allocated yet.'
        stop  
      endif

      HFBHamil(sb+1:sb+2*N+2*N2, sb+1:sb+2*N+2*N2) = ConstructHFBHamil(        &
      &                           sphamil(si+1:si+N+N2,si+1:si+N+N2),          &
      &                           HFBgaps(si+1:si+N+N2,si+1:si+N+N2), N, N2,   &
      &                           rho_history(si+1:si+N+N2,si+1:si+N+N2),      &
      &                           kappa_history(si+1:si+N+N2,si+1:si+N+N2),    &
      &                           HFBgauge(it))  

      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
    enddo

    !---------------------------------------------------------------------------
    ! b) We repeatedly diagonalize the matrix to find a suitable Fermi energy, 
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
         blocked_qps(NN+i) = p_blocked(i) + 2*sum(HFBlocks(1:4))
       enddo
    endif
    !---------------------------------------------------------------------------
    ! Naively writing down the HFB hamiltonian, we have 
    ! (if there is a linear, antihermitian conserved symmetry)  
    ! 
    !
    !       (  h+  0     0     d+- )
    !  H =  (  0   h-    d-+   0   )
    !       (  0  -d-+  -h+    0   )
    !       ( -d+- 0     0    -h-  ) 
    !
    !
    !  but we diagonalize in practice the two submatrices
    !
    !  H+ =  (  h+     d+-)     H- = ( h-     d-+ )
    !        (  -d+-  -h- )          ( d-+   -h+  )
    !
    !  as
    !         H = ( H+ 0 )
    !             ( 0  H-)
    !     
    ! Hence, the Bogoliubov transformation in memory is structured as
    !        
    !       (  V^*+  U+    0     0   )
    !  W =  (  U^*+  V+    0     0   )
    !       (  0     0     V^*-  U-  )
    !       (  0     0     U^*-  V-  )
    !   
    !  which we need to correct by moving things around to
    !
    !       (  V^*+   0     U+  0   )
    !  W =  (  0      V^*-  0   U-  )      (*)
    !       (  0      U^*-  0   V-  )
    !       (  U^*+   0     V+  0   )
    !
    ! Note that this reordering is necessary too for the
    !  (i)  QPenergies
    !  (ii) Configmatrix
    ! 
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -    
    ! Almost(!) all of these matrices are ordered by increasing value of 
    ! quasi-particle energy. So, in the Bogoliubov transformation at (*),
    ! the columns correspond to the following ordering of qp energies   
    !
    !  -E_-N ... -E_-1 , -E_+N ... -E_+1, E_+1, ... E_+N, E_-1 .... E_-N
    !
    !  where E_+/-1 .... E_+/-N are increasing sequences, and I draw attention
    !  to the fact that the first half of quasi-particles is organised in 
    !  reverse order from the second half. 
    !
    ! The EXCEPTION to this rule is the configuration matrix configmatrix. That
    ! is constructed in blocks of increasing ABSOLUTE quasi-particle energy
    ! meaning that in that matrix the ordering is  
    !   
    ! -E_+1 ... -E_+N , -E_-1 ... -E_-N, E_+1, ... E_+N, E_-1 .... E_-N
    !
    ! meaning that the order has been reversed for the first half. This is 
    ! convenient for the construction of the pairing matrices, as then the  
    ! full generalized density matrix in quasi-particle representation 
    !
    ! looks like
    !          ( f_+1  0                     ......                    0  )
    !    R =   ( 0    ....                                                ) 
    !          (            f_-1                                          ) 
    !          (                  .....                                   )
    !          (                         1-f_+1                           ) 
    !          (                                  ......                  ) 
    !          (                                          1-f_-1          )
    !          ( 0                                                 .......)
    !
    !
    ! These warnings about matrix ordering of the first half of the Bogoliubov
    ! transformation is a little bit academical: the code has been constructed
    ! such that only the right half of the Bogoliubov transformation in (*), 
    ! their quasi-particle energies and the FULL configuration matrix enter
    ! the generalized density matrix, which is the only quantity affecting the
    ! rest of the program. 
    !---------------------------------------------------------------------------
    temp = Bogoliubov ;  tempqe = QPenergies  ; tempc        = configmatrix
    Bogoliubov = 0    ;  QPenergies = 0.0d0   ; configmatrix = 0.0

    si = 0 ; sb = 0
    do B=1,8,2
      N = HFBlocks(B) ; N2 = HFBlocks(B+1)

      !-------------------------------------------------------------------------
      ! Moving the first block

      ! First N eigenvalues of  H+ 
      ! (negative qp energies generally, but not always)
      Bogoliubov(sb+       1:sb  +N   ,sb+1:sb+N) = &
      &                                        temp(sb+  1:sb+  N,sb+  1:sb+N)
      Bogoliubov(sb+N+2*N2+1:sb+2*N+2*N2,sb+1:sb+N) = &
      &                                        temp(sb+N+1:sb+2*N,sb+  1:sb+N)

      ! Second set of N eigenvalues of H+
      Bogoliubov(sb+       1:sb+  N     ,sb+N+  N2+1:sb+2*N+  N2) = &
      &                                        temp(sb+  1:sb+  N,sb+N+1:sb+2*N)
      Bogoliubov(sb+N+2*N2+1:sb+2*N+2*N2,sb+N+  N2+1:sb+2*N+  N2) = &
      &                                        temp(sb+N+1:sb+2*N,sb+N+1:sb+2*N)
      !-------------------------------------------------------------------------
      ! Moving the second block (which doesn't exist if T is conserved)

      ! First N2 eigenvalues of H-
      ! (negative qp energies generally, but not always)
      Bogoliubov(sb+N+1:sb+N+2*N2,sb+N+1:sb+N+N2) = &
      &                          temp(sb+2*N+1:sb+2*N+2*N2,sb+2*N+1:sb+2*N+N2)
      ! Second set of N2 eigenvalues of H-
      Bogoliubov(sb+N+1:sb+N+2*N2,sb+2*N+N2+1:sb+2*N+2*N2) = &
      &                       temp(sb+2*N+1:sb+2*N+2*N2,sb+2*N+N2+1:sb+2*N+2*N2)

      !-------------------------------------------------------------------------
      ! First N eigenvalues of H+
      configmatrix(sb+     1:sb+N    )  = tempc (sb  +1:sb+  N)
      ! Second set of N eigenvalues of H+
      configmatrix(sb+N+N2+1:sb+2*N+N2) = tempc (sb+N+1:sb+2*N) 
      !-------------------------------------------------------------------------
      ! First set of N2 eigenvalues of H-
      configmatrix(sb+N+     1:sb+N+N2  )   = tempc (sb+2*N+1   :sb+2*N+  N2)
      ! Second set of N2 eigenvalues of H-
      configmatrix(sb+2*N+N2+1:sb+2*N+2*N2) = tempc (sb+2*N+N2+1:sb+2*N+2*N2)

      !-------------------------------------------------------------------------
      qpenergies  (sb+     1:sb+N)      = tempqe(sb    +1:sb+N   ) 
      qpenergies  (sb+N+   1:sb+N+N2)   = tempqe(sb+2*N+1:sb+2*N+N2) 

      qpenergies  (sb+  N+N2+1:sb+2*N+  N2) = tempqe(sb+  N   +1:sb+2*N   )
      qpenergies  (sb+2*N+N2+1:sb+2*N+2*N2) = tempqe(sb+2*N+N2+1:sb+2*N+2*N2)

      !-------------------------------------------------------------------------
      si = si +   N +   N2
      sb = sb + 2*N + 2*N2
    enddo

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

    ! Calculating overlaps    
    overrho = sum(rho_pairing * rho_history)/sum(rho_pairing**2)
    overkap = sum(kappa_pairing * kappa_history)/sum(kappa_pairing**2)

    !---------------------------------------------------------------------------
    ! e) Optionally mix these densities
    if(.not.all(rho_history.eq.0.0)) then
      if(HFBmixtype .eq. 0) then
        rho_pairing   =  HFBmix * rho_pairing   + (1.0d0-HFBmix) * rho_history
        kappa_pairing =  HFBmix * kappa_pairing + (1.0d0-HFBmix) * kappa_history
      endif
    endif  
    
    !---------------------------------------------------------------------------
    ! f) Calculate the dispersion, while we are at it. 

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

    si            = 0
    HFBdispersion = 0.0
    allocate(chi(nwt, nwt))
    do B=1,8
      N = HFBlocks(B)  ;  if(N .eq. 0) cycle 

      it  = 1
      if(B .gt. 4) it = 2
      
      ! rho squared
      chi(si+1:si+N, si+1:si+N) = matmul(rho_pairing(si+1:si+N, si+1:si+N),    &
      &                                      rho_pairing(si+1:si+N, si+1:si+N) )    

      !                                       Tr rho - Tr rho^2
      do i=1,N
          HFBdispersion(it) = HFBdispersion(it) + rho_pairing(si+i, si+i)      &
          &                                     - chi(si+i, si+i)              &
          &                                     + kappa_pairing(si+i, si+i)**2
      enddo
      si = si + N
    enddo
    deallocate(chi)

    ! Time-reversal introduces a factor of two
$TR    HFBdispersion = 2 * HFBdispersion 

  end subroutine solvepairing_HFB

  function ConstructConfiguration(Bogo, Eqp, blocks, blocktype, blockconf,     &
  &                               blocked_qp)   result(R)
    !---------------------------------------------------------------------------
    ! Construct the configuration matrix, based on the various user options.
    !
    !---------------------------------------------------------------------------

    integer, intent(in)          :: BlockType
    integer, intent(in)          :: blocks(4) 
    integer, intent(in)          :: blockconf(:)
    integer, allocatable         :: blocked_qp(:)
    real(KIND=dp), intent(in)    :: Eqp(:), Bogo(:,:)
    real(KIND=dp), allocatable   :: R(:)
      
    integer                      :: N, N2,B, sb, i, NB, j, qblock, ind, si, bi
    real(KIND=dp)                :: compare, occ, qpmin
    integer                      :: toblock(4), qpb

    N = size(Eqp) 
    allocate(R(N)) ;  R = 0
    qpb = -1

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
        if(.not.allocated(blocked_qp)) then
          allocate(blocked_qp(NB)) ; blocked_qp = 0
        endif        

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
                N = blocks(B)
                if(qblock.eq.B) then
                    R(sb+ind-N)       = occ 
                    R(sb+ind)         = 1 - occ
                    ! Save which one we blocked
                    blocked_qp(j) = sb+ind
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
            do B=1,4,2
              N = blocks(B) ; if (N.eq.0) cycle
              N2= blocks(B+1) 
              if(Eqp(sb+N+toblock(B)+1) .lt. qpmin) then
                qpmin = Eqp(sb+N+toblock(B)+1)
                qpb   = B
              endif
              si = si +   N + N2
            enddo
            toblock(qpb) = toblock(qpb) + 1
          enddo
        endif

        NB = sum(toblock)
        if(allocated(blocked_qp)) then
          deallocate(blocked_qp)
        endif

        if(.not.allocated(blocked_qp)) then
          allocate(blocked_qp(NB)) ; blocked_qp = 0
        endif        

        !  For every block, we flip the required number of qps.
        sb = 0 ; si = 0 ; ind = 0
        do B=1,4
          N = blocks(B) ; if (N.eq.0) cycle
          do j=1,toblock(B)
            ! The qps are ordered in energy from the diagonalization
            !  So we simply flip the first ones
            R(sb + N + j ) = 1 - occ
            R(sb     + j ) =     occ

            ! Saving the one we flipped
            ind = ind + 1
            blocked_qp(ind) = sb + j
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

  subroutine FindFermi_secant(H, blocks, targetparticles, config, Bogo, Eqp,   & 
            &       lambda, maxhfbiter, blocktype, blockconf, blocked_qp, ifail)
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
      integer, intent(out)         :: ifail
      integer, allocatable         :: blocked_qp(:)

      real(KIND=dp)                :: df, dn(2), particles
      integer                      :: iter

      ! Initialization
      df = 0 ; dn = 0.0
      do iter=1, maxHFBiter

        particles = & 
        &   diagbyblock(H,blocks, config, Bogo,Eqp,lambda, blocktype,blockconf,&
        &               blocked_qp, ifail)
        ! Return if we do not want to readjust the Fermi energy
        if(MaxHFBiter.eq.1) return
        
        !-----------------------------------------------------------------------
        ! Readjust the Fermi energy based on the number of particles.
        ! We use the secant method.
        dn(2) = dn(1)
        dn(1) = particles - targetparticles

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
  end subroutine FindFermi_secant

  function diagbyblock(H, blocks, config, Bogo,Eqp,lambda, blocktype,blockconf,&
  &                    blocked_qp, ifail) result(particles)
      !-------------------------------------------------------------------------
      ! Routine that diagonalizes, block by block, a HFB Hamiltonian that is 
      ! passed in. It does the low-level work for all the high-level routines
      ! in this moodule. 
      !
      ! Input
      !   H        : HFB hamiltonian, without Fermi energy
      !   blocks   : Sizes of the symmetry blocks that can be used to simplify 
      !              the problem.
      !   lambda   : Fermi energy.
      !
      ! Output
      !   config   : Configuration matrix of the final solution
      !   Eqp      : Quasiparticle energies of the final solution
      !   Bogo     : Bogoliubov transformation that diagonalizes H
      !   Lambda   : Final fermi energy
      ! 
      !   particles: total number of particles for this input
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in)    :: H(:,:)
      real(KIND=dp), intent(out)   :: config(:), Bogo(:,:), Eqp(:)
      real(KIND=dp), intent(inout) :: lambda
      integer, intent(in)          :: blocks(4), blocktype
      integer, intent(in)          :: blockconf(:)
      integer, allocatable         :: blocked_qp(:)

      real(KIND=dp), allocatable   :: eigen(:), work(:), A(:,:)
      real(KIND=dp)                :: particles
      integer                      :: sb, si, N, B, i, ifail, lwork

      allocate(eigen(2*sum(blocks)))
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
        lwork = -1; allocate(work(1))
        call DSYEV( 'V', 'U', 2*N, A, 2*N, eigen(sb+1:sb+2*N),work,lwork,ifail)
        lwork = int(work(1)); deallocate(work) ; allocate(work(lwork))
        call DSYEV( 'V', 'U', 2*N, A, 2*N, eigen(sb+1:sb+2*N),work,lwork,ifail)
        deallocate(work)

        Bogo(sb+1:sb+2*N, sb+1:sb+2*N) = A

        if(ifail.ne.0) then
          print *, 'WARNING: diagon failed in subroutine DiagByBlock.'
          print *, '         Problematic block B = ', B
          deallocate(A, eigen)
          particles = 0.0
          return
        endif

        Eqp(sb+1:sb+N) = eigen(sb+1:sb+2*N)
        deallocate(A)
        ! Indices for the next block
        si = si +   N
        sb = sb + 2*N
      enddo
      !-----------------------------------------------------------------------
      ! b) We construct the configuration matrix that was asked for
      config = ConstructConfiguration(Bogo, Eqp, blocks, blocktype, blockconf, &
      &                               blocked_qp) 
      !-----------------------------------------------------------------------
      ! c) Count the total number of particles that we have.
      si = 0 ; sb = 0
      particles   = 0
      do B=1,4                
          N = Blocks(B) ;  if(N .eq. 0) cycle 
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
$TR   particles = 2 * particles                 

      deallocate(eigen)
  end function diagbyblock

  subroutine FindFermi_Brent(H, blocks, targetparticles, config, Bogo, Eqp,    & 
   &                lambda, maxhfbiter, blocktype, blockconf, blocked_qp, ifail)
      !-------------------------------------------------------------------------
      ! Subroutine that diagonalizes the HFB hamiltonian (repeatedly) to find  
      ! the correct Fermi energy that fixes the average number of particles.
      ! The routine only solves this for one particular isospin.
      !
      ! This particular subroutine employs Brents method to fix the Fermi 
      ! energy. See
      ! 
      ! https://en.wikipedia.org/wiki/Brent%27s_method
      ! 
      ! which combines bisection, secant method and inverse quadratic 
      ! interpolation.The original source is probably
      ! R. P. Brent (1973), "Chapter 4: An Algorithm with Guaranteed Convergence
      ! for Finding a Zero of a Function", Algorithms for Minimization without
      ! Derivatives, Englewood Cliffs, NJ: Prentice-Hall,  
      !
      ! Input
      !   H        : HFB hamiltonian, without Fermi energy
      !   blocks   : Sizes of the symmetry blocks that can be used to simplify 
      !              the problem.
      !   lambda   : Initial guess for the Fermi energy
      !   maxhfbiter: Maximum number of iterations to perform
      !   targetparticles: Average number of particles to target. 
      !
      ! Output
      !   config   : Configuration matrix of the final solution
      !   Eqp      : Quasiparticle energies of the final solution
      !   Bogo     : Bogoliubov transformation that diagonalizes H
      !   Lambda   : Final fermi energy
      !-------------------------------------------------------------------------
      ! This routine is very heavily inspired/copy-pasted by the routines 
      ! implemented in MOCCa by M. Bender. 
      !-------------------------------------------------------------------------
      real(KIND=dp), intent(in)    :: H(:,:), targetparticles
      real(KIND=dp), intent(out)   :: config(:), Bogo(:,:), Eqp(:)
      real(KIND=dp), intent(inout) :: lambda
      integer, intent(in)          :: blocks(4), maxhfbiter, blocktype
      integer, intent(in)          :: blockconf(:)
      integer, intent(out)         :: ifail  

      real(KIND=dp)                :: InitialBracket(2), FA, FB, N
      integer                      :: idir = 0 , idirsig = 1, FailCount
      logical                      :: Success
      integer, allocatable         :: blocked_qp(:)
 
      !-------------------------------------------------------------------------
      ! STEP 1: set up an initial bracket
      !-------------------------------------------------------------------------
      N = diagbyblock(H, blocks, config, Bogo,Eqp,lambda, blocktype,blockconf, &
      &               blocked_qp, ifail)
      N = N - targetparticles
      ! Check if this guess for lambda is good enough
      if(abs(N).lt.pairing_prec) return

      ! Use present Fermi energy as starting point and check the direction
      ! where the zero of <N>-N0 can be expected.
      ! If <N>-N0 <  0, search at higher values.
      ! If <N>-N0 >= 0, search at lower  values.
      ! Initialize InitialBracket(it,1) = A, InitialBracket(it,2) = B with 
      ! present  Fermi energy.
      ! "dir" is the label of the InitialBracket(it,idir) that has to be moved,
      ! "idirsig" is the sign of steps needed to go into that direction.

      InitialBracket(:) = lambda
      if ( N .lt. 0.0_dp ) then 
        idir   =  2 ;  idirsig =  1
      else 
        idir   =  1 ;  idirsig = -1
      endif

      ! Try to find a boundary that brackets the Fermi energy in the direction 
      ! into which the Fermi energy has to be changed.
      FailCount = -1 ;  Success = .false.

      do while(.not. Success)
        FailCount = FailCount + 1
        
        ! update moving boundary and recalculate particle numbers at both.
        InitialBracket(idir) = &
        &                 InitialBracket(idir) + idirsig * 0.01_dp*(FailCount+1)

        FA = diagbyblock(H,blocks,config,Bogo,Eqp,InitialBracket(1),blocktype, &
        &                                          blockconf, blocked_qp, ifail)
        FB = diagbyblock(H,blocks,config,Bogo,Eqp,InitialBracket(2),blocktype, &
        &                                          blockconf, blocked_qp, ifail)
        FA = FA - targetparticles ; FB = FB - targetparticles

        ! check if N(epsilon_F) is a monotonically growing function.
        ! It should be, but who knows, pigs may fly ...
        if ( FB .lt. FA ) then 
          print '(" : Warning N(eps_F) decreases ")'
          print '(" A = ",f13.8," FA = ",f14.8," B = ",f13.8," FB = ",f14.8)', &
               & InitialBracket(1),FA+N, InitialBracket(2),FB+N
        endif

        ! diagnostic printing for convergence analysis (usually commented out)
!        print '(" Bracketing ",i4,1l2,(2(f13.8,es16.7)))',        &
!              & FailCount,Success,InitialBracket(1),FA, InitialBracket(2),FB          
        ! code failure (Fermi energy has changed by 30 MeV)
        if (Failcount .gt. 76) then
          print '(/," A = ", f13.8, "FA = ",1es12.4,              &
               &    " B = ", f13.8, "FB = ",1es12.4)',            &
               &     InitialBracket(1),FA,InitialBracket(2),FB 
          ifail = 1
          return
          !stop 'FindFermiBrent: Search for InitialBracket failed.'
        endif
        ! check if root is bracketed for isospin it after the update
        if( FA*FB .lt. 0.0_dp ) then 
            ! Correct Bracket found!
            Success = .true.
        endif
      enddo
      !-------------------------------------------------------------------------
      ! STEP 2: call the routine for the actual bisection
      !-------------------------------------------------------------------------
      call BrentBisection(lambda,N,InitialBracket(1), InitialBracket(2),FA,FB, &
      &                   maxHFBIter,H,blocks, targetparticles, config, Bogo,  & 
      &                   Eqp, blocktype, blockconf, blocked_qp)
  
  end subroutine FindFermi_brent

  subroutine BrentBisection(lambda,particles, X1,X2,FX1,FX2,Depth, H, blocks,  & 
    &      targetparticles, config, Bogo, Eqp, blocktype, blockconf, blocked_qp)
    !---------------------------------------------------------------------------
    ! This routine searches for the Fermi energy
    ! by Brent's methods https://en.wikipedia.org/wiki/Brent%27s_method
    ! which combines bisection, secant method and inverse quadratic 
    ! interpolation. The original source is probably
    ! R. P. Brent (1973), "Chapter 4: An Algorithm with Guaranteed Convergence
    ! for Finding a Zero of a Function", Algorithms for Minimization without
    ! Derivatives, Englewood Cliffs, NJ: Prentice-Hall, 
    !---------------------------------------------------------------------------
    ! see pages 1188 - 1189 of http://apps.nrbook.com/fortran/index.html
    ! W. H. Press, S. A. Teukolsky, W. T. Vetterling and B. P. Flannery,
    ! Numerical Recipes in Fortran in Fortran 90, Second Edition (1996).
    !---------------------------------------------------------------------------
    ! Input
    !   H        : HFB hamiltonian, without Fermi energy
    !   blocks   : Sizes of the symmetry blocks that can be used to simplify 
    !              the problem.
    !   maxhfbiter: Maximum number of iterations to perform
    !   targetparticles: Average number of particles to target. 
    !
    ! Output
    !   config   : Configuration matrix of the final solution
    !   Eqp      : Quasiparticle energies of the final solution
    !   Bogo     : Bogoliubov transformation that diagonalizes H
    !   Lambda   : Final fermi energy
    !   Particles: Final number of particles
    !---------------------------------------------------------------------------
    
    real(KIND=dp), intent(in)    :: H(:,:), targetparticles
    real(KIND=dp), intent(out)   :: config(:), Bogo(:,:), Eqp(:), lambda
    real(KIND=dp), intent(out)   :: particles
    integer, intent(in)          :: blocks(4), blocktype
    integer, intent(in)          :: blockconf(:)
    integer, intent(in)          :: Depth
    real(KIND=dp), intent(in)    :: X1 , X2, FX1 , FX2 
    integer, allocatable         :: blocked_qp(:)

    real(KIND=dp)                :: A , B, C , FA, FB , FC
    real(KIND=dp)                :: D , E, S , P  , Q , R 
    real(KIND=dp)                :: Num , Tol , XM 
    real(KIND=dp)                :: eps = 1.d-9
    integer                      :: FailCount, ifail
    logical                      :: Found

    A  = X1 ; B  = X2 
    FA = FX1; FB = FX2
    Found = .false.    
    if (A .eq. B) then 
      !-------------------------------------------------------------------------
      ! This signals that FA = FB is zero within the tolerance.
      ! Either near-converged HFB or HF case of completely broken-down pairing
      ! which also satisfies FA = FB = 0 within an interval. The 
      ! latter case cannot be handled by the algorithm below.
      !-------------------------------------------------------------------------
      Found = .true.
    endif

    C = B ; FC = FB 
    E = -1000000 ; D = -1000000
  
    FailCount = -1

    do while(.not.Found) 
      FailCount = FailCount + 1
      if ( ( FB .gt. 0.0_dp .and. FC .gt. 0.0_dp ) .or. & 
         & ( FB .lt. 0.0_dp .and. FC .lt. 0.0_dp ) )  then
        C  = A     ;  FC = FA
        D  = B - A ;  E  = D
      endif
      if ( abs(FC) .lt. abs(FB) ) then
        A  = B ;  FA = FB
        B  = C ;  FB = FC
        C  = A ;  FC = FA
      endif
      !-------------------------------------------------------------------------
      ! Convergence check
      ! Note (W.R.): I have tightened convergence a bit compared to the values
      !              in MOCCa by M.B. 
      !-------------------------------------------------------------------------
      Tol  = 2.0_dp * eps * abs(B) + 0.05_dp * Pairing_prec
      XM   = 0.5_dp * (C-B)
      !----------------------------------------------------------------
      ! Note: the tolerance is on the precision of the Fermi energy,
      ! NOT the nearness of the particle number to the targeted value.
      !----------------------------------------------------------------
      if ( abs(XM) .le. Tol .or. FB .eq. 0.0_dp ) then
        Lambda =  B
        Found = .true. 
        cycle
      endif
      if ( abs(E) .ge. Tol .and. abs(FA) .gt. abs(FB) ) then
        S = FB/FA
        if ( A .eq. C ) then
          P = 2.0_dp * XM * S
          Q = 1.0_dp - S
        else
          Q = FA/FC
          R = FB/FC
          P = S * (2.0_dp * XM * Q * (Q-R) & 
                  &    - (B-A)*(R-1.0_dp))
          Q = (Q-1.0_dp)*(R-1.0_dp)*(S-1.0_dp)
        endif
        if ( P .gt. 0.0_dp ) Q = -Q
        P = abs(P)
        if (2.0_dp * P .lt. min(3.0_dp*XM*Q - abs(Tol*Q),abs(E*Q))) then
          E = D
          D = P / Q
        else
          D = XM
          E = D 
        endif
      else
        D = XM
        E = D 
      endif
      A  = B 
      FA = FB
      B  = B + merge(D,sign(Tol,XM),abs(D) .gt. Tol)    
  
      !-------------------------------------------------------------------------
      ! B is present best guess for the fermi energy, FB the corresponding 
      ! particle number.
      !-------------------------------------------------------------------------
      Num = diagbyblock(H, blocks, config, Bogo,Eqp,B, blocktype,blockconf,    &
      &                 blocked_qp, ifail)
      FB  = Num - targetparticles

      !-------------------------------------------------------------------------
      ! diagnostic printing for convergence analysis (usually commented out)
      !-------------------------------------------------------------------------
      ! NOTE: B is the the best guess for the zero of F. A has been the previous
      ! "closest" interval boundary that is not updated after Found
      ! is set to .true. The actual zero might therefore be outside the 
      ! interval [A,B]. If so, the true zero is typically closer to B than 
      ! A is to B.
      !-------------------------------------------------------------------------
      ! Note further: as A and B are swapped from time to time, B might be 
      ! smaller than A when printed here
      !-------------------------------------------------------------------------
!      print '(" BrentBisection ",i4,(1l2,2(f13.8,es16.7),f14.8))',    &
!           & FailCount, Found,A,FA,B,FB,Num
!     
      if ( FailCount .gt. Depth ) then
        print '(/," Warning: BrentBisection did not converge after ",i4," iterations")', & 
        &      FailCount
      endif
    enddo
    ! Output
    Lambda    = B ; particles = FB
  end subroutine BrentBisection

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
    integer                   :: wave1, wave2, iso, si,  B, N, inda, indb, N2
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

          HFBgaps(inda,indb) =     sum(hfpsi(:,:,indb)*deltapsi)*dv *          &
          &                          Pcutoffs(inda)*Pcutoffs(indb)*stabfac

          ! The full matrix Delta is antisymmetric...
$NTR      HFBgaps(indb,inda) =  - HFBgaps(inda,indb)
          ! ... but the stored matrix is symmetric when time-reversal is 
          ! conserved; but this is not exploited at the moment!
!$TR       HFBgaps(indb,inda) = HFBgaps(inda,indb)

        enddo
      enddo

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
    4 format (' overlap rho,kap = ',  2es9.2)

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
    print 4, overrho, overkap

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


!      do i=1, N+N2
!        print ('(99f10.3)'), tmp(i, 1:N+N2)
!      enddo
!      print *

!      do i=1, N
!        print *, 'CONJUG', i, conjugp(si+i), si+i+N, conjugp(conjugp(si+i)), si+i, kappa_can(si+i)
!      enddo
!      print *

      deallocate(tmp)
      si = si + N + N2
    enddo    
    
   end subroutine Canonical
   
   subroutine ConstructCanonicalBasis(Transfo)
    !---------------------------------------------------------------------------
    ! Transform the spwf wavefunctions in the HFBasis into Canbasis, with the 
    ! passed in Transfo. 
    !---------------------------------------------------------------------------

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
      !-------------------------------------------------------------------------
      ! Apply the transformation in this symmetry block
      do wave1=1, N 
        do wave2=1,N
          CanPsi(:,:,si+wave1)  = CanPsi(:,:,si+wave1) +                       &
          &                     Transfo(si+wave2,si+wave1) * HFPsi(:,:,si+wave2) 
          
          canenergies(si+wave1) = canenergies(si+wave1) +                      &
          &             abs(Transfo(si+wave2,si+wave1)**2) *spenergies(si+wave2) 
        enddo 
      enddo
      
      si = si +  N
      !-------------------------------------------------------------------------
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

   subroutine clean_HFB
    if(allocated(HFBgaps))  deallocate(HFBGaps)
   end subroutine clean_HFB

end module
