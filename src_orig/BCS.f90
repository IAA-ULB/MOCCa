module BCS
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
 ! Module implementing the routines for the solution of the BCS equations.
 !
 !==============================================================================

 use compilation
 use geninfo
 use wavefunctions
 use pairingcutoffs

 implicit none
 
 !------------------------------------------------------------------------------
 real(KIND=dp) :: FermiPrec = 1d-9
 integer       :: maxBCSiter = 500
 logical       :: ConstantGap = .false.
 
 real(KIND=dp) :: BCSdispersion(2)
 !------------------------------------------------------------------------------
 ! Values of the BCS gaps Delta.
 real(KIND=dp), allocatable :: BCSGaps(:)
 ! BCS quasiparticle energies
 real(KIND=dp), allocatable :: BCSqps(:)
 !------------------------------------------------------------------------------
 ! BCS occupations, i.e. 2 * v_i^2 (factor 2 due to time-reversal)
 ! - - - - - - - - - - - - - - - --  - - - - - - - - - - - - - - - - - - - - - -
 ! IMPORTANT NOTE:
 !  These are NOT ALWAYS equal to the diagonal matrix elements of the 
 !  density matrix. This is ONLY true at zero temperature. They are ALWAYS
 !  equal to 2 * v_i**2, but NOT always equal to rho_ii !
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: BCSoccupations(:)
 !------------------------------------------------------------------------------
 ! BCS occupation factors f, typically defined in the context of 
 ! finite-temperature calculations, but also useful in the context of 
 ! Equal Filling Approximation blocking calculations.
 real(KIND=dp), allocatable :: BCSf(:)
 !------------------------------------------------------------------------------
 procedure(delta_action_dummy), pointer :: delta_action_BCS

  interface
   function delta_action_dummy(psi,&
$N1DELTA                    &      dpsi, &
$N2DELTA                    &            ddpsi, &
$N3DELTA                    &                   dddpsi, &
$SYMDELTA                   &                          sx,sy,sz, &
&                                                               iso, onthefly) &
                                                                result(deltapsi)
      !-------------------------------------------------------------------------
      ! Dummy function to allow this module to acces the functional.f90 module 
      ! to acces the information on the acces of deltas.
      ! Note that the actual delta_action routine's interface is decided by 
      ! Hephaestos at compiletime, and as such this dummy interface has to also
      ! be decided at that time.
      !-------------------------------------------------------------------------
      real*8, intent(in)    :: psi(:,:)  
$N1DELTA      real*8, intent(inout) ::   dpsi(:,:,:)
$N2DELTA      real*8, intent(inout) ::  ddpsi(:,:,:)
$N3DELTA      real*8, intent(inout) :: dddpsi(:,:,:)
$SYMDELTA     integer, intent(in)   :: sx(:),sy(:),sz(:)
      integer, intent(in)   :: iso
      real*8, allocatable   :: deltapsi(:,:)
      logical, intent(in)   :: onthefly
   end function
  end interface

contains
 
 subroutine uv_from_occupation(occ, u, v)
  !-----------------------------------------------------------------------------
  ! Simple function that calculates BSC u and v factors in a numerically
  ! safe way from the BCS occupation v^2.
  !
  ! Attention: only VALID if u and v can be assumed to be real.
  !
  ! Input:
  !      occ  : 2 v^2, i.e. what is stored in BCSoccupations
  ! Output:
  !      u, v : bcs factors
  !-----------------------------------------------------------------------------

  real(KIND=dp), intent(in) :: occ
  real(KIND=dp), intent(out) :: u,v

  if(occ/2 .gt. 0.0d0) then
    v = sqrt(occ/2)
  else
    v = 0.0d0
  endif
  
  u = 1 - occ/2
  if(u .gt. 0.0d0) then
    u = sqrt(u)
  else
    u = 0.0d0
  endif
 end subroutine uv_from_occupation
 
 subroutine solvepairing_BCS(fermi, rho_can, kappa_can, qpenergies, gas,       &
 &                             BlockType,Blockindices, blocklowest, blocked_qps)
  !-----------------------------------------------------------------------------
  ! Driver routine for the solving of the BCS equations.
  !
  ! |----
  ! | Do until the Fermi energy is stationary
  ! |   1) Calculate the BCS quasiparticle energies
  ! |        using the matrix elements of h and delta already calculated!
  ! |   2) Construct the correct configuration with either
  ! |       *) ordinary, ground-state BCS
  ! |       *) finite-temperature BCS
  ! |       *) equal-filling blocking BCS
  ! |   3) Calculate the Fermi energy with analytic formula
  ! |----
  !     3) Calculate the entries in rho_pairing and kappa_pairing
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(inout) :: fermi(2)
  real(KIND=dp), intent(inout) :: rho_can(:)
  real(KIND=dp), intent(inout) :: kappa_can(:), qpenergies(:)
  real(KIND=dp)                :: oldfermi(2), fac
  integer                      :: iter, wave
  integer, intent(in)          :: gas

  ! Options for the determination of a blocking configuration
  integer, intent(in)          :: Blockindices(:)
  integer, intent(in)          :: BlockType
  character(len=2), intent(in) :: BlockLowest(:)
  integer, allocatable         :: blocked_qps(:)
  
  1 format('---------------------------------------',/,     &
    &        ' Warning! ',/,                                &
    &        ' BCS calculations did not converge.',/,       &
    &        ' Pairing iterations:  ', i4,/,                &
    &        ' Old Fermi:    ', 2f12.7,/,                   &
    &        ' New Fermi:    ', 2f12.7,/,                   &
    &        '---------------------------------------')

  if(.not.allocated(Bcsf)) then
    allocate(bcsf(nwt)) ; bcsf = 0
  endif

  !-----------------------------------------------------------------------------
  ! Iteration start
  do iter =1, maxBCSiter
    oldfermi = fermi
    ! Calculate the quasiparticle energies
    call BCSQPEnergies(Fermi) 
    BCSf = BCSconstructconfiguration(gas, BlockType,Blockindices, blocklowest, & 
                                                                    blocked_qps)
    call BCSFindFermiEnergy(Fermi)
    ! Check for convergence
    if( all(abs(fermi - oldfermi).lt.FermiPrec)) then
      exit
    elseif(iter.eq.maxBCSiter) then
      print 1, iter, oldfermi, fermi
    endif          
  enddo
  ! Iteration end  
  call calcBCSoccupations(Fermi)
  !-----------------------------------------------------------------------------
  ! Rho_pairing is diagonal for a BCS calculation
  rho_can = 0.0
  do wave=1,nwt
    ! Occupations are 
    !   n_a = f_i + v_i^2 (1 - 2 * f_i)
    fac = BCSf(wave) 
    ! Note that there is already a factor two due to timereversal in    
    ! BCSoccupations, but not in the first term in the formula above.
    rho_can(wave) = 2.0*fac + BCSoccupations(wave) * (1 - 2.0*fac)

    ! Failsafe
    if(rho_can(wave).lt.0.0) then
       rho_can(wave) = 0.0
    endif
  enddo

  !-----------------------------------------------------------------------------
  ! Kappa is in its canonical form for a BCS calculation.
  ! However, we are not storing the time-reversed partners in this case, so
  ! put the matrix elements (i,ibar) on the diagonal anyway.
  !
  ! u * v  = 0.5 * Delta/(sqrt(epsilon**2 + Delta**2))
  !-----------------------------------------------------------------------------
  kappa_can = 0.0
  do wave=1,nwt
     ! At finite temperature, the elements of kappa are
     ! kappa_i\bar{i} = u_i v_i ( 1 - 2 * f_i )
     fac             = BCSf(wave)
     kappa_can(wave) = 0.5 * BCSgaps(wave)/(BCSqps(wave)) * (1 - 2.0*fac)
  enddo

  ! Qpenergies in this case are the BCSqpenergies
  Qpenergies = BCSqps

  ! Side effect: calculate the dispersion 
  call calcBCSdispersion(rho_can, kappa_can)

 end subroutine solvepairing_BCS
 
 subroutine CalcBCSGaps(fermi, stabfactor)
    !---------------------------------------------------------------------------
    ! Calculate the BCS pairing gaps.
    !---------------------------------------------------------------------------
    integer                      :: wave, iso, wave_global
    real(KIND=dp)                :: deltapsi(mv,4), trash(2)
    real(KIND=dp), intent(in)    :: fermi(2), stabfactor(2)
#if(USE_MPI>0)
    integer                      :: mpi_err
#endif
    ! trash statement to stop the compiler complaining about unused dummy 
    ! variables
    trash = fermi

    if(ConstantGap) then  
      ! Constantgap pairing
      do wave=1,nwt
          BCSGaps(wave) = 2.0 * PCutoffs(wave)**2
      enddo
    else
       ! Use the delta_action function to calculate the matrix elements
       BCSgaps = 0.0d0 ! zeroing to be able to call MPI_ALLREDUCE later
       do wave=1,nwt_local
            wave_global = spwf_map(wave)
            if(wave_global .le. nwn) then
                iso = -1
            else
                iso = +1
            endif

            deltapsi = delta_action_BCS(  hfpsi(:,:,wave) ,             &
$N1DELTA    &                            hfdpsi(:,:,:,wave),            &
$N2DELTA    &                           hfddpsi(:,:,:,wave),            &
$N3DELTA    &                          hfdddpsi(:,:,:,wave),            &
$SYMDELTA   &              sx(:,wave), sy(:,wave), sz(:,wave),          &
            &                                                iso,.false.)
 
 
            BCSgaps(wave_global) =  sum(hfpsi(:,:,wave)*deltapsi)*dv*          &
            &             Pcutoffs(wave_global)**2 * (1 + stabfactor((iso+3)/2))
       enddo
#if(USE_MPI>0)
       call MPI_ALLREDUCE(MPI_IN_PLACE, BCSgaps,nwt,MPI_REAL8,MPI_SUM,         &
       &                                                 MPI_COMM_WORLD,mpi_err)
#endif
    endif
  end subroutine CalcBCSGaps

  subroutine BCSFindFermiEnergy (Fermi)
   !----------------------------------------------------------------------------
   ! Implements the analytical formula for the Fermi energy, for fixed 
   ! quasiparticle energies and occupation factors f.
   !  
   !
   ! Starting by demanding that
   !     N     = 2 sum_(i>0) [ f_i  + v_i^2 ( 1 - 2 f_i) ]
   ! and  
   !     v_i^2 = 1/2 ( 1 - (eps_i - mu)/E^qp_i ) 
   !
   ! we obtain 
   ! 
   !  mu = (N - a)/b
   !  
   ! with
   !   a = sum_(i>0) (1  - eps_i/E^qp_i * (1-2f_i))  
   !   b = sum_(i>0) ((1 - 2 f_i)/E^qp_i)  
   !----------------------------------------------------------------------------
    real(KIND=dp),intent(inout) :: Fermi(2)
    real(KIND=dp)               :: LambdaSums(2,2), eqp, nom, Particles(2), fac
    integer                     :: wave, it

    Particles(1) = Neutrons; Particles(2) = Protons

    LambdaSums=0.0_dp 
    !Sum the quasiparticle energies to get the correct lambda.
    do wave=1,nwt
      it = 1
      if(wave .gt. nwn) it = 2 
     
      eqp = BCSqps(wave)
      nom = spenergies(wave)

      fac =  BCSf(wave) 
      lambdasums(1,it)= lambdasums(1,it) +(1.0_dp - nom/eqp*(1-2*fac))
      lambdasums(2,it)= lambdasums(2,it) +       1.0_dp/eqp*(1-2*fac)         
    enddo
    do it=1,2
        Fermi(it)=(particles(it)  - lambdasums(1,it))/lambdasums(2,it)
    enddo
  end subroutine BCSFindFermiEnergy

  function BCSConstructConfiguration(gas, blocktype,Blockindices, blocklowest, & 
  &                                  blocked_qps) result(f)

    !---------------------------------------------------------------------------
    ! Construct the correct pairing configuration as a function of various 
    ! user options. 
    !  (a) Zero-temperature, no blocking  : all the BCS-f factors are zero.
    !  (b) Zero-temperature, EFA blocking : f_k = 1/2 for the blocked state
    !  (c) Finite-temperature, no blocking: the f_k are given by a fermi
    !                                       function (modulo different options 
    !                                       for the BCS gas)
    !  (d) Finite-temperature, blocking   : nothing implemented
    !---------------------------------------------------------------------------
    real*8              :: f(nwt), occ, qpmin
    integer             :: wave, NB, i, ind, si, N, B,  qpb, c

    integer, intent(in)          :: Blockindices(:)
    integer, intent(in)          :: BlockType, gas
    integer, allocatable         :: proton_block(:), neutron_block(:)
    integer, allocatable         :: blocked_qps(:), indices(:), toblock(:)
    character(len=2), intent(in) :: BlockLowest(:)
  
    f = 0 ; qpb = 0
    if(allocated(blocked_qps)) deallocate(blocked_qps)
    !---------------------------------------------------------------------------
    ! Zero-temperature
    if(inversetemp .lt. 0) then
      occ = 0.0
      select case(Blocktype)
      case(0)
        ! No blocking, all the f are zero
        return
      case(1,2,5)
        ! Time-reversal breaking blocking asked for, impossible to do in BCS
        call stp('The code cannot perform true blocking in BCS.')
      case(3,4,6)
        ! Equal filling blocking
        occ = 0.5d0
      end select

      select case(Blocktype)
      case(3)
        !-----------------------------------------------------------------------
        ! We search for a specific configuration, i.e. a quasiparticle with
        ! a specific sp index. In a BCS calculation, this is trivial. 
        NB = size(blockindices)      
        allocate(blocked_qps(NB))
        do i=1, NB
          ind    = blockindices(i)
          f(ind) = occ
          blocked_qps(i) = ind
        enddo
      case(4)
        !-----------------------------------------------------------------------
        ! We search for the lowest qp with specific quantum numbers.
        
        allocate(proton_block(5))  ; proton_block  = 0 
        allocate(neutron_block(5)) ; neutron_block = 0
        allocate(toblock(8))       ; toblock       = 0

        ! Count all qps that were asked for
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
        toblock(1:4) = neutron_block(1:4)
        toblock(5:8) = proton_block(1:4)
        ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        ! First see where the lowest qp energies are for those qps in block(5)
        if(neutron_block(5).ne.0) then
          do i = 1, neutron_block(5)
            qpmin = 10000000
            si    = 0
            qpb   = 0
            do B=1,4
              N = HFBlocks_global(B) ; if (N.eq.0) cycle
              if(bcsqps(si+toblock(B)+1) .lt. qpmin) then
                qpmin = bcsqps(si+toblock(B)+1)
                qpb   = B
              endif
              si = si +   N
            enddo
            toblock(qpb) = toblock(qpb) + 1
          enddo
        endif

        if(proton_block(5).ne.0) then
          do i = 1, proton_block(5)
            qpmin = 10000000
            si    = sum(HFBlocks_global(1:4)) 
            qpb   = 0
            do B=5,8
              N = HFBlocks_global(B) ; if (N.eq.0) cycle

              if(bcsqps(si+toblock(B)+1) .lt. qpmin) then
                qpmin = bcsqps(si+toblock(B)+1)
                qpb   = B
              endif
              si = si +   N
            enddo
            toblock(qpb) = toblock(qpb) + 1
          enddo
        endif

        allocate(blocked_qps(sum(toblock))) ; blocked_qps = 0
        si = 0
        c  = 0
        do B=1,8
            N = HFBlocks_global(B) ; if(N.eq.0) cycle
            indices = Order(BCSqps(si+1:si+N))
            do i=1, toblock(B)
              f(si+indices(i)) = occ
              c = c+1
              blocked_qps(c) = si+indices(i)
            enddo            
            si = si + N
        enddo
        deallocate(proton_block, neutron_block) 
        !-----------------------------------------------------------------------
      ! functionality woth modelspwfs removed since it was not useful
      !case(6)
      !  !-----------------------------------------------------------------------
      !  ! We search for the spwf with the largest overlap with the model
      !  ! wavefunction
      !  allocate(blocked_qps(1)) ; blocked_qps = 0!!

      !  B = modelblock
      !  if(B.gt.1) then
      !    si = sum(HFBlocks_global(1:B-1))
      !  else
      !    si = 0
      !  endif
      !  N       = HFBlocks_global(B)
      !  maxover = -10
      !  indover =   0
      !  do i=1, N
!     !       overlap = abs(sum(HFpsi(:,:,si+i) * modelspwf)) * dv
      !      if(overlap .gt. maxover) then
      !        maxover = overlap
      !        indover = i
      !      endif
      !  enddo
      !  f(si + indover) = occ
      !  blockoverlap    = maxover
      !  blocked_qps(1)  = si + indover
      end select
    !---------------------------------------------------------------------------
    ! Finite-temperature
    else
      if(blocktype.ne. 0) then
        call stp('Cannot do finite-temperature BCS with blocking.')
      endif
      
      ! Select occupations based on the type of treatment of the gas
      select case(gas)
      case(0)
        ! No special treatment of the gas
        f = 1./(1. + exp(inversetemp * BCSqps))
      case(1)
        ! Not implemented!
        call stp('gastype = 1 is not implemented in the BCS module.')
      case(2)
        ! Only take into account the bound states
        do wave=1,nwt
          if(spenergies(wave) .lt. 0) then
            f(wave) = 1./(1. + exp(inversetemp * BCSqps(wave)))
          else 
            f(wave) = 0
          endif
        enddo
      case DEFAULT
        call stp('Unknown value for "gas" in the BCS module.')
      end select
    endif
    return
  end function BCSConstructConfiguration
  
  subroutine BCSQPEnergies(Fermi)
   !----------------------------------------------------------------------------
   ! This subroutine calculates the BCS quasiparticle energies as a function of
   ! the Fermi energies and the pairing gaps.
   !
   ! eqp_k = sqrt[(epsilon_k - \lambda)**2 + Delta_{k, kbar}**2]
   !
   ! which is formula (6.72) on page 235 in Ring & Shuck.
   !----------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Fermi(2)
    integer                   :: wave, it
    real(KIND=dp)             :: epsilon, lambda

    if(.not.allocated(BCSqps)) then
      allocate(BCSqps(nwt)) ; BCSqps=0
    endif
   
    do wave=1,nwt
        it = 1
        if(wave .gt. nwn) it = 2

        epsilon = spenergies(wave)
        lambda  = Fermi(it)
        BCSqps(wave)  = sqrt((epsilon - lambda)**2 + abs(BCSgaps(wave))**2)
    enddo

  end subroutine BCSQPEnergies

  subroutine calcBCSOccupations(Fermi)
   !----------------------------------------------------------------------------
   ! Find the occupation numbers of the HFBasis from the quasiparticle energies
   ! and the Fermi energies.
   !
   ! v^2_k = 0.5 * (1 - (Epsilon - Lambda)/(E_{qp}))
   !
   ! or formula (6.51) on page 231 in Ring & Schuck. Note that these are only
   ! entries of the density matrix in the case of zero-temperature calculations.
   !----------------------------------------------------------------------------
    real(KIND=dp), intent(in)   :: Fermi(2)
    integer                     :: wave,it
    real(KIND=dp)               :: eqp
    
    if(.not.allocated(BCSoccupations)) then
      allocate(BCSoccupations(nwt)) ; BCSoccupations = 0.0
    endif
    
    do wave=1,nwt
      it  = 1
      if(wave .gt. nwn) it = 2
      
      eqp = BCSqps(wave)

      ! Zero-temperature BCS
      BCSOccupations(wave) = 0.5*(1 - (spenergies(wave) - Fermi(it))/eqp)
    enddo
    ! Time reversal symmetry
    BCSOccupations = 2 * BCSOccupations

   end subroutine calcBCSOccupations

   subroutine calcBCSdispersion(rho_can, kappa_can)
      !-------------------------------------------------------------------------
      ! Calculate the dispersion 
      !  Tr rho - rho^2
      !-------------------------------------------------------------------------
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
      real(KIND=dp), intent(in):: rho_can(:), kappa_can(:)
      integer :: wave

      ! BCS dispersion
      BCSdispersion = 0.0    
      do wave=1,nwn
        BCSdispersion(1) = BCSdispersion(1)                                    &
        &                             + 0.5*rho_can(wave)*(1-0.5*rho_can(wave))&
        &                             + kappa_can(wave)**2
      enddo

      do wave=nwp+1,nwt
        BCSdispersion(2) = BCSdispersion(2)                                    &
        &                             + 0.5*rho_can(wave)*(1-0.5*rho_can(wave))&
        &                             + kappa_can(wave)**2
      enddo
      BCSdispersion  = 2 * BCSdispersion

   end subroutine calcBCSdispersion

   function average_gap_BCS() result(gap)
      !-------------------------------------------------------------------------
      ! Calculation of two types of "average gap", based on 
      !
      ! M. Bender et al., EPJA 8, 59-75 (2000).
      !
      ! <v2 Delta > = sum f_k v^2_k   Delta_k / sum f_k v^2_k
      ! <uv Delta > = sum f_k u_k v_k Delta_k / sum f_k u_k v_k
      !
      ! Note that the f_k in the reference is the square of our cutoff!
      !
      ! I calculate the average gap here WITHOUT cutoffs, as our definition 
      ! of the gaps already includes all of the cutoff factors already, in
      ! contrast to the EPJA paper. 
      !-------------------------------------------------------------------------

      real(KIND=dp) :: gap(2,2), norm(2,2), v2, uv
      integer       :: it, wave

      gap = 0 ; norm = 0
      if(.not.allocated(Pcutoffs)) return      
      do wave=1,nwt
        it = 1
        if(wave.gt.nwn) it = 2
      
        uv = 0.5 * BCSgaps(wave)/(BCSqps(wave))
        v2 =       BCSoccupations(wave)

        ! Note that the definition of the gaps include the cutoff factors. 
        !  v^2 weighted 
        gap(1,it) = gap(1,it)    + v2 * BCSgaps(wave) !*  Pcutoffs(wave)**2
        norm(1,it)= norm(1,it)   + v2                 !*  Pcutoffs(wave)**2 
        ! uv weighted
        gap(2,it) = gap(2,it)    + uv * BCSgaps(wave) !*  Pcutoffs(wave)**2
        norm(2,it)= norm(2,it)   + uv                 !*  Pcutoffs(wave)**2
      enddo
      gap = gap/norm
      
   end function average_gap_BCS

   subroutine clean_BCS()

     if(allocated(BCSgaps))        deallocate(BCSgaps) 
     if(allocated(BCSqps))         deallocate(BCSqps) 
     if(allocated(BCSoccupations)) deallocate(BCSoccupations)
     if(allocated(BCSf))           deallocate(BCSf)
 
   end subroutine clean_BCS

   function Order(energies) result(Indices)
    !---------------------------------------------------------------------------
    ! Returns the indices for an ordered traversal of the input array.
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Input:
    !  energies :  real*8, a set of energies to be ordered
    ! Output:
    !  indices  :  integer, the indices to get the energies in ascending order  
    !---------------------------------------------------------------------------
    integer, allocatable       :: Indices(:)
    real(Kind=dp),intent(in)   :: Energies(:)
    real(Kind=dp),allocatable  :: Eswap(:)
    integer                    :: i, nwf,  HolePos, ToInsertIndex
    real(Kind=dp)              :: ToInsert

    nwf = size(energies)
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Filling Energies & Indices
    if(allocated(indices))  deallocate(indices)
    allocate(Indices(nwf), Eswap(nwf))
    do i=1,nwf
       Indices(i) = i 
    enddo

    Eswap = Energies
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !Sort the energies
    do i=2,nwf
      !Make a hole at index i
      ToInsert = Eswap(i)
      HolePos  = i
      ToInsertIndex = Indices(i)
      do while(ToInsert.lt.Eswap(HolePos-1))
        !Move the hole one place down
        Eswap(HolePos) = Eswap(HolePos-1)
        Indices(HolePos) = Indices(HolePos-1)
        HolePos = HolePos - 1
        if(HolePos.eq.1.0_dp) exit
      enddo
      !Insert the energy at the correct place
      Eswap(HolePos)    = ToInsert
      Indices(HolePos)  = ToInsertIndex
    enddo

    deallocate(Eswap)
  end function Order
!===============================================================================
!  Never to be used function to define an interface for delta_action
!===============================================================================   
! function delta_action_dummy(psi, dpsi, ddpsi, dddpsi, sx,sy,sz,iso, onthefly) &
!                                                                result(deltapsi)
!    !---------------------------------------------------------------------------
!    ! Dummy function to allow this module to acces the functional.f90 module 
!    ! to acces the information on the acces of deltas.
!    !---------------------------------------------------------------------------
!    
!    real(KIND=dp), intent(in)    :: psi(mv,4)  
!    real(KIND=dp), intent(inout) :: dpsi(mv,3,4),ddpsi(mv,6,4), dddpsi(mv,10,4)
!    integer, intent(in)       :: sx(4),sy(4),sz(4),   iso
!    logical, intent(in)       :: onthefly
! end function
  
end module
