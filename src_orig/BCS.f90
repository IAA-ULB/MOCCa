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
 ! BCS occupations, i.e. 2 * v_i^2 
 ! (factor 2 due to time-reversal)
 ! - - - - - - - - - - - - - - - --  - - - - - - - - - - - - - - - - - - - - - -
 ! IMPORTANT NOTE:
 !  These are NOT ALWAYS equal to the diagonal matrix elements of the 
 !  density matrix. This is ONLY true at zero temperature. They are ALWAYS
 !  equal to 2 * v_i**2, but NOT always equal to rho_ii !
 !------------------------------------------------------------------------------
 real(KIND=dp), allocatable :: BCSoccupations(:)
 !------------------------------------------------------------------------------
 ! BCS occupation factors f
 !  f_i = 1/(1 + exp(beta * E_qp))
 real(KIND=dp), allocatable :: BCSf(:)
 !------------------------------------------------------------------------------
 procedure(delta_action_dummy), pointer :: delta_action_BCS

 interface

  function delta_action_dummy(psi, dpsi, ddpsi, dddpsi, sx,sy,sz,iso, onthefly)&
                                                                result(deltapsi)
    !---------------------------------------------------------------------------
    ! Dummy function to allow this module to acces the functional.f90 module 
    ! to acces the information on the acces of deltas.
    !---------------------------------------------------------------------------
    
    real*8, intent(in)    :: psi(:,:)  
    real*8, intent(inout) :: dpsi(:,:,:),ddpsi(:,:,:), dddpsi(:,:,:)
    integer, intent(in)   :: sx(:),sy(:),sz(:),iso
    logical, intent(in)   :: onthefly
    real*8, allocatable   :: deltapsi(:,:)
   end function
 end interface

contains
 
 subroutine solvepairing_BCS(fermi, rho_can, kappa_can, qpenergies,gas)
  !-----------------------------------------------------------------------------
  ! Driver routine for the solving of the BCS equations.
  !
  ! |----
  ! | Do until the Fermi energy is stationary
  ! |   1) Calculate the BCS quasiparticle energies
  ! |        using the matrix elements of h and delta already calculated!
  ! |   2) Calculate the Fermi energy with closed formula
  ! |----
  !     3) Calculate the entries in rho_pairing and kappa_pairing
  !-----------------------------------------------------------------------------
  real(KIND=dp), intent(inout) :: fermi(2)
  real(KIND=dp), intent(inout) :: rho_can(:)
  real(KIND=dp), intent(inout) :: kappa_can(:), qpenergies(:)
  real(KIND=dp)                :: oldfermi(2), fac
  integer                      :: iter, wave
  integer, intent(in)          :: gas
  
  1 format('---------------------------------------',/,     &
    &        ' Warning! ',/,                                &
    &        ' BCS calculations did not converge.',/,       &
    &        ' Pairing iterations:  ', i4,/,                &
    &        ' Old Fermi:    ', 2f12.7,/,                   &
    &        ' New Fermi:    ', 2f12.7,/,                   &
    &        '---------------------------------------')

  do iter =1, maxBCSiter
    oldfermi = fermi
    call BCSQPEnergies(Fermi, gas)
    call BCSFindFermiEnergy(Fermi, gas)
    
    ! Check for convergence
    if( all(abs(fermi - oldfermi).lt.FermiPrec)) then
      exit
    elseif(iter.eq.maxBCSiter) then
      print 1, iter, oldfermi, fermi
    endif          
  enddo
  
  call calcBCSoccupations(Fermi, gas)
  
  !-----------------------------------------------------------------------------
  ! Rho_pairing is diagonal for a BCS calculation
  rho_can = 0.0
  do wave=1,nwt

    if(inversetemp.eq.-1) then
      rho_can(wave) = BCSoccupations(wave)
      if(rho_can(wave).lt.0.0) then
         rho_can(wave) = 0.0
      endif
    else
      ! Occupations are 
      !   n_a = f_i + v_i^2 (1 - 2 * f_i)
      fac = BCSf(wave) 
      ! Note that there is already a factor two due to timereversal in    
      ! BCSoccupations, but not in the first term in the formula above.
      rho_can(wave) = 2.0*fac + BCSoccupations(wave) * (1 - 2.0*fac)
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
    if(inversetemp.eq.-1) then
      kappa_can(wave) = 0.5 * BCSgaps(wave)/(BCSqps(wave))
    else
      ! At finite temperature, the elements of kappa are
      ! kappa_i\bar{i} = u_i v_i ( 1 - 2 * f_i )
      fac             = BCSf(wave)
      kappa_can(wave) = 0.5 * BCSgaps(wave)/(BCSqps(wave)) * (1 - 2.0*fac)
    endif
  enddo

  ! Qpenergies in this case are the BCSqpenergies
  Qpenergies = BCSqps

  ! Side effects, calculate the dispersion 
  call calcBCSdispersion(rho_can, kappa_can)

 end subroutine solvepairing_BCS
 
 subroutine CalcBCSGaps(fermi, stabfactor)
    !---------------------------------------------------------------------------
    ! Calculate the BCS pairing gaps.
    !---------------------------------------------------------------------------
    integer                      :: wave, iso
    real(KIND=dp)                :: deltapsi(mv,4)
    real(KIND=dp), intent(in)    :: fermi(2), stabfactor(2)
    
    if(ConstantGap) then  
      ! Constantgap pairing
      do wave=1,nwt
          BCSGaps(wave) = 2.0 * PCutoffs(wave)**2
      enddo
    else
      ! Use the delta_action to calculate the elements in the gaps
       do wave=1,nwt
            if(wave .le. nwn) then
                iso = -1
            else
                iso = +1
            endif
            
            if(.not.associated(Delta_action_BCS)) stop
             
            deltapsi = delta_action_BCS(  hfpsi(:,:,wave)  ,                   &
            &                            hfdpsi(:,:,:,wave),                   &
            &                           hfddpsi(:,:,:,wave),                   &
            &                          hfdddpsi(:,:,:,wave),                   &
            &              sx(:,wave), sy(:,wave), sz(:,wave),iso,.false.)
 
            BCSgaps(wave) =     sum(hfpsi(:,:,wave)*deltapsi)*dv*              &
            &               Pcutoffs(wave)**2 * (1 + stabfactor((iso+3)/2))
       enddo
    endif
  end subroutine CalcBCSGaps

  subroutine BCSFindFermiEnergy (Fermi, gas)
   !----------------------------------------------------------------------------
   ! Function that finds the correct Fermi energy to satisfy the particle 
   ! number constraints 
   !----------------------------------------------------------------------------
    real(KIND=dp),intent(inout) :: Fermi(2)
    real(KIND=dp)               :: LambdaSums(2,2), eqp, nom, Particles(2), fac
    integer                     :: wave, it
    integer, intent(in)         :: gas

    Particles(1) = Neutrons; Particles(2) = Protons

    LambdaSums=0.0_dp 
    !Sum the quasiparticle energies to get the correct lambda.
    do wave=1,nwt
      it = 1
      if(wave .gt. nwn) it = 2 
     
      eqp = BCSqps(wave)
      nom = spenergies(wave)

      if(inversetemp.eq.-1) then  
          lambdasums(1,it)= lambdasums(1,it) +(1.0_dp - nom/eqp)
          lambdasums(2,it)= lambdasums(2,it) +       1.0_dp/eqp
      else
        select case(gas)
        case(0)
          ! Ordinary finite-temperature BCS
          fac =  BCSf(wave) !1.0/(1 + exp(inversetemp * eqp))
          lambdasums(1,it)= lambdasums(1,it) +(1.0_dp - nom/eqp*(1-2*fac))
          lambdasums(2,it)= lambdasums(2,it) +       1.0_dp/eqp*(1-2*fac)  
        case(1)
          ! Subtraction method for gas degrees of freedom          

        case(2)
          ! Taking into account only bound states
          fac =  BCSf(wave) !1.0/(1 + exp(inversetemp * eqp))

          if(nom .lt. 0) then
            ! Only the bound states contribute
            lambdasums(1,it)= lambdasums(1,it) +(1.0_dp - nom/eqp*(1-2*fac))
            lambdasums(2,it)= lambdasums(2,it) +       1.0_dp/eqp*(1-2*fac)  
          endif
        end select

      endif

    enddo
    ! There is a nice analytical formula for the Fermi energy in the BCS case
    do it=1,2
        Fermi(it) =  (particles(it) - lambdasums(1,it))/lambdasums(2,it)
    enddo
  end subroutine BCSFindFermiEnergy
  
  subroutine BCSQPEnergies(Fermi, gas)
   !----------------------------------------------------------------------------
   ! This subroutine calculates the BCS quasiparticle energies as a function of
   ! the Fermi energies and the pairing gaps.
   !
   ! eqp_k = sqrt[(epsilon_k - \lambda)**2 + Delta_{k, kbar}**2]
   !
   ! which is formula (6.72) on page 235 in Ring & Shuck.
   !----------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Fermi(2)
    integer, intent(in)       :: gas
    integer                   :: wave, it
    real(KIND=dp)             :: epsilon, lambda

    if(.not.allocated(BCSqps)) then
      allocate(BCSqps(nwt)) ; BCSqps=0
      allocate(BCSf(nwt)) ; BCSf = 0.0
    endif
   
    do wave=1,nwt
        it = 1
        if(wave .gt. nwn) it = 2

        epsilon = spenergies(wave)
        lambda  = Fermi(it)
        BCSqps(wave)  = sqrt((epsilon - lambda)**2 + abs(BCSgaps(wave))**2)
    enddo
    
    if(inversetemp.ne.-1) then
      select case(gas)
      case(0)
        BCSf = 1./(1. + exp(inversetemp * BCSqps))
      case(1)
  
      case(2)
        do wave=1,nwt
          if(spenergies(wave) .lt. 0) then
            BCSf(wave) = 1./(1. + exp(inversetemp * BCSqps(wave)))
          else 
            BCSf(wave) = 0
          endif
        enddo
      end select
    else
      BCSf = 0
    endif

  end subroutine BCSQPEnergies

  subroutine calcBCSOccupations(Fermi, gas)
   !----------------------------------------------------------------------------
   ! Find the occupation numbers of the HFBasis from the quasiparticle energies
   ! and the Fermi energies.
   !
   ! v^2_k = 0.5 * (1 - (Epsilon - Lambda)/(E_{qp}))
   !
   ! or formula (6.51) on page 231 in Ring & Schuck.
   !----------------------------------------------------------------------------
    real(KIND=dp), intent(in)   :: Fermi(2)
    integer, intent(in)         :: gas
    integer                     :: wave,it
    real(KIND=dp)               :: eqp
    
    if(.not.allocated(BCSoccupations)) then
      allocate(BCSoccupations(nwt)) ; BCSoccupations = 0.0
    endif
    
    do wave=1,nwt
        it  = 1
        if(wave .gt. nwn) it = 2
        
        eqp = BCSqps(wave)

        if(inversetemp.eq.-1) then
          ! Zero-temperature BCS
          BCSOccupations(wave) = 0.5*(1 - (spenergies(wave) - Fermi(it))/eqp)
        else
          ! Finite temperature BCS
          select case(gas)
          case(0,1)
            BCSOccupations(wave) = 0.5*(1 - (spenergies(wave) - Fermi(it))/eqp)
          case(2)
            if(spenergies(wave).gt.0) then
              BCSoccupations(wave) = 0
            else
              BCSoccupations(wave) = 0.5*(1 - (spenergies(wave) - Fermi(it))/eqp)
            endif
          end select
        endif
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
        gap(1,it) = gap(1,it)    + v2 * BCSgaps(wave) *  Pcutoffs(wave)**2
        norm(1,it)= norm(1,it)   + v2                 *  Pcutoffs(wave)**2 
        ! uv weighted
        gap(2,it) = gap(2,it)    + uv * BCSgaps(wave) *  Pcutoffs(wave)**2
        norm(2,it)= norm(2,it)   + uv                 *  Pcutoffs(wave)**2
      enddo
      gap = gap/norm
      
   end function average_gap_BCS

   subroutine clean_BCS()

     if(allocated(BCSgaps))        deallocate(BCSgaps) 
     if(allocated(BCSqps))         deallocate(BCSqps) 
     if(allocated(BCSoccupations)) deallocate(BCSoccupations)
     if(allocated(BCSf))           deallocate(BCSf)
 
   end subroutine clean_BCS
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
