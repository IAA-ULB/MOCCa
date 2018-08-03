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
 integer       :: maxBCSiter = 100
 logical       :: ConstantGap = .false.
 
 !------------------------------------------------------------------------------
 ! Values of the BCS gaps Delta.
 real(KIND=dp), allocatable :: BCSGaps(:)
 ! BCS quasiparticle energies
 real(KIND=dp), allocatable :: BCSqps(:)
  
 real(KIND=dp), allocatable :: BCSoccupations(:)
 
 !------------------------------------------------------------------------------
 procedure(delta_action_dummy), pointer :: delta_action_BCS

contains
 
 subroutine solvepairing_BCS(fermi, rho_pairing, kappa_pairing)
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
  real(KIND=dp), intent(inout) :: rho_pairing(nwt,nwt), kappa_pairing(nwt,nwt)
  real(KIND=dp)                :: oldfermi(2)
  integer                      :: iter, wave
  
  1 format('---------------------------------------',/,     &
    &        ' Warning! ',/,                                &
    &        ' BCS calculations did not converge.',/,       &
    &        ' Pairing iterations:  ', i4,/,                &
    &        ' Old Fermi:    ', 2f12.7,/,                   &
    &        ' New Fermi:    ', 2f12.7,/,                   &
    &        '---------------------------------------')

  do iter =1, maxBCSiter
    oldfermi = fermi
    call BCSQPEnergies(Fermi)
    call BCSFindFermiEnergy(Fermi)

    ! Check for convergence
    if( all(abs(fermi - oldfermi).lt.FermiPrec)) then
      exit
    elseif(iter.eq.maxBCSiter) then
      print 1, iter, oldfermi, fermi
    endif          
  enddo
  
  ! Compute the cutoffs: this is not needed at this point, but it is needed 
  ! for the calculation of the pairdensities afterwards
  call ComputePairingCutoffs(fermi)
  
  call calcBCSoccupations(Fermi)
  !-----------------------------------------------------------------------------
  ! Rho_pairing is diagonal for a BCS calculation
  rho_pairing = 0.0
  do wave=1,nwt
    rho_pairing(wave,wave) = BCSoccupations(wave)
  enddo
 
  !-----------------------------------------------------------------------------
  ! Kappa is in its canonical form for a BCS calculation.
  ! However, we are not storing the time-reversed partners in this case, so
  ! put the matrix elements (i,ibar) on the diagonal anyway.
  !
  ! u * v  = 0.5 * Delta/(sqrt(epsilon**2 + Delta**2))
  !
  kappa_pairing = 0.0
  do wave=1,nwt
    kappa_pairing(wave,wave) = 0.5 * BCSgaps(wave)/(BCSqps(wave))
  enddo

 end subroutine solvepairing_BCS
 
 subroutine CalcBCSGaps(fermi)
    !---------------------------------------------------------------------------
    ! Calculate the BCS pairing gaps.
    ! Currently only does constant gap. 
    !---------------------------------------------------------------------------
    integer                      :: wave, iso
    real(KIND=dp)                :: deltapsi(mv,4)
    real(KIND=dp), intent(in)    :: fermi(2)
    
    call ComputePairingCutoffs(fermi)

    if(ConstantGap) then  
      ! Constantgap pairing
      do wave=1,nwt
          BCSGaps(wave) = 2.0 * PCutoffs(wave)**2
      enddo
    else
      ! Use the delta_action to calculate the elements in the gaps
       do wave=1,nwt
            if(wave .lt. nwn) then
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
 
            BCSgaps(wave) = sum(hfpsi(:,:,wave)*deltapsi)*dv * Pcutoffs(wave)**2
            
            print *, 'GAPS', BCSgaps(wave), BCSgaps(wave)/2
       enddo
    endif

  end subroutine CalcBCSGaps

  subroutine BCSFindFermiEnergy (Fermi)
   !----------------------------------------------------------------------------
   ! Function that finds the correct Fermi energy to satisfy the particle 
   ! number constraints 
   !----------------------------------------------------------------------------
    real(KIND=dp),intent(inout) :: Fermi(2)
    real(KIND=dp)               :: LambdaSums(2,2), eqp, nom, Particles(2)
    integer                     :: wave, it

    Particles(1) = Neutrons; Particles(2) = Protons

    LambdaSums=0.0_dp
    !Sum the quasiparticle energies to get the correct lambda.
    do wave=1,nwt
      it = 1
      if(wave .gt. nwn) it = 2 
     
      eqp = BCSqps(wave)
      nom = spenergies(wave)

      lambdasums(1,it)= lambdasums(1,it) +(1.0_dp - nom/eqp)
      lambdasums(2,it)= lambdasums(2,it) + 1.0_dp/eqp
    enddo
    ! There is a nice analytical formula for the Fermi energy in the BCS case
    do it=1,2
        Fermi(it) =  (Particles(it) - lambdasums(1,it))/lambdasums(2,it)
    enddo
  end subroutine BCSFindFermiEnergy
  
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
   ! v^2_k = 0.5 * 1 - (Epsilon - Lambda)/(E_{qp})
   !
   ! or formula (6.51) on page 231 in Ring & Schuck.
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
        !There is already an intrinsic factor two here due to Time-reversal
        BCSOccupations(wave) = (1 - (spenergies(wave) - Fermi(it))/eqp)
    enddo

   end subroutine calcBCSOccupations

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
  
end module
