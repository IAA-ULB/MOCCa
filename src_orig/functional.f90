module functional
 !==============================================================================
 !  #######   ##   #    # #####   ##   #      #    #  ####
 !     #     #  #  ##   #   #    #  #  #      #    # #
 !     #    #    # # #  #   #   #    # #      #    #  ####
 !     #    ###### #  # #   #   ###### #      #    #      #
 !     #    #    # #   ##   #   #    # #      #    # #    #
 !     #    #    # #    #   #   #    # ######  ####   ####
 !
 !  Copyright W. Ryssens & M. Bender
 !
 !==============================================================================
 !
 ! Module containing the means to calculate (and print) the mean-field energy.
 ! Note that the actual coupling constants are contained in the constants.f90
 ! file. 
 !==============================================================================
 
 use compilation
 use geninfo
 use densities
 
 implicit none
 
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! All the different Energies
 real(KIND=dp) :: Kinetic(2)     ,CoulombEnergy
 real(KIND=dp) :: CoulombExchange,CoMCorrection(2,2)
 real(KIND=dp) :: TotalEnergy    ,OldEnergy(7)
 real(KIND=dp) :: SpEnergy       ,LNENergy(2)
 real(KIND=dp) :: Routhian       ,OldRouthian(7)
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! Different forms of the Skyrme energy
 ! Note the indexing
 !    first  index = indicates term
 !    second index = indicates if the term is LO, NLO, N2LO or N3LO
 !
 ! Energies are calculated in the BFH representation using the B coupling 
 ! coefficients, and then recombined into the isospin representation using 
 ! the C coefficients. By default only the latter is printed, but the BFH
 ! representation can be asked for for debugging purposes. 
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 real(KIND=dp) :: SkyrmeTerms_C(11,3), SkyrmeTerms_B(11,3)
 
 
 procedure(PrintEnergy_interface), pointer :: PrintEnergy
 
 abstract interface
    subroutine PrintEnergy_interface(Lagrange)
      logical, intent(in), optional :: Lagrange
    end subroutine
  end interface
 
 contains
 
 subroutine CompEnergy()
    !---------------------------------------------------------------------------
    ! Subroutine that computes all the different energies of the main program
    ! state.
    !---------------------------------------------------------------------------

    integer :: i
    
    ! Making sure the PrintEnergy routine is associated
    if(.not.associated(PrintEnergy)) then
        !PrintEnergy => PrintEnergy_termbyterm
    endif

    !Initialise the contributions to the energy
    CoulombEnergy=0.0_dp; CoMCorrection=0.0_dp ; Kinetic      =0.0_dp
    SkyrmeTerms_C=0.0_dp; SkyrmeTerms_B=0.0_dp

    !Shift the entries in OldEnergy by one place and put in
    !the previous value of the Energy.
    do i=0,5
        OldEnergy(7-i) = OldEnergy(6-i)
    enddo
    OldEnergy(1)=TotalEnergy

    !call CompKinetic
    
!    SkyrmeTerms = compSkyrme(Density)
!    TotalEnergy = sum(SkyrmeTerms)

!    ! Calculate the N2LO terms
!    N2LOterms = N2LO(Density)

!    ! Calculate the N3LO term
!    N3LOterms = N3LO(Density)

    !Pairing Energy
    !PairingEnergy = CompPairingEnergy(Delta)
    !Lipkin-Nogami Energy
    !if(Lipkin) then
    !  LNEnergy      = - LNLambda * PairingDisp
    !endif

    !Calculating the CoulombEnergy
    !CoulombEnergy   = CompCoulombEnergy(Density)
    !CoulombExchange = CompCoulombExchange(Density)

    !COM Correction
    !call CompCOMCorrection

!    !Sum of all energies (Note that the Skyrme sum already was included)
!    TotalEnergy = TotalEnergy        + sum(Kinetic)  + CoulombEnergy  +        &
!    &             sum(CoMCorrection) + CoulombExchange + sum(PairingEnergy) +  &
!    &             sum(LNEnergy) + sum(N2LOterms)

!    !Calculate the Routhian too
!    do i=0,5
!        OldRouthian(7-i) = OldRouthian(6-i)
!    enddo
!    OldRouthian(1) = Routhian
!    Routhian = TotalEnergy                                                     &
!    !                         Contribution of multipole moment constraints
!    &                      +    sum(ConstraintEnergy*Density%Rho)*dv           &
!    !                         Contribution of the cranking constraints (Omega*J)
!    &                      +    sum(CrankEnergy)

!    !Energy due to the single particle states.
!    SpEnergy = SpwfEnergy()
    return
  end subroutine CompEnergy
 
end module functional
