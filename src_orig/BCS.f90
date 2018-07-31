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
 !
 integer       :: maxBCSiter = 100
  
 !------------------------------------------------------------------------------
 ! Values of the BCS gaps Delta.
 real(KIND=dp), allocatable :: BCSGaps(:)
 ! BCS quasiparticle energies
 real(KIND=dp), allocatable :: BCSqps(:)
  
 real(KIND=dp), allocatable :: BCSoccupations(:)

contains
 
 subroutine solvepairing_BCS(fermi, rho_pairing, kappa_pairing)
  !-----------------------------------------------------------------------------
  ! Driver routine for the solving of the BCS equations.
  !
  ! |----
  ! | Do until the Fermi energy is stationary
  ! |   1) Determine the pairing gaps Delta_{ij} using the pairing fields.
  ! |      from the fields module.
  ! |   2) Calculate the BCS quasiparticle energies with those gaps
  ! |   3) Calculate the Fermi energy with closed formula
  ! |----
  !     4) Calculate the entries in rho_pairing and kappa_pairing
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
    
    call ComputePairingCutoffs(fermi)
    
    oldfermi = fermi

    call calcBCSGaps()
    call BCSQPEnergies(Fermi)
    call BCSFindFermiEnergy(Fermi)
    ! Check for convergence
    
    if( all(abs(fermi - oldfermi).lt.FermiPrec)) then
      exit
    elseif(iter.eq.maxBCSiter) then
      print 1, iter, oldfermi, fermi
    endif          
  enddo
  
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
  kappa_pairing = 0.0
  do wave=1,nwt
    kappa_pairing(wave,wave) = BCSgaps(wave)/(BCSqps(wave))
  enddo

 end subroutine solvepairing_BCS
 
 subroutine CalcBCSGaps()
    !---------------------------------------------------------------------------
    ! Calculate the BCS pairing gaps.
    ! Currently only does constant gap. 
    !---------------------------------------------------------------------------
    integer                      :: wave

    if(.not.allocated(BCSGaps)) then
        allocate(BCSGaps(nwt)) ; BCSGaps = 0.0
    endif

    do wave=1,nwt
        BCSGaps(wave) = 2.0 * PCutoffs(wave)**2
    enddo

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
  
end module
