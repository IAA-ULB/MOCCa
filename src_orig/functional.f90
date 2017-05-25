module functional
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

  function CompKinetic() result(kinetic)
    !---------------------------------------------------------------------------
    ! This subroutine computes the total kinetic energy,
    ! according to the following formula:
    !    E_k = -\hbar/2m \int d^3x \sum_{k} v_{k} \Psi_k^* \Delta \Psi_k
    !---------------------------------------------------------------------------
    ! Note that the 1-body c.o.m. correction is not taken into account here!
    !---------------------------------------------------------------------------
    use Constants

    integer          :: wave, it,k,i
    real(KIND=dp)    ::  Inproduct
    real(KIND=dp)    :: Kinetic(2)

    ! Kinetic Energy
    Kinetic = 0.0_dp
    do wave=1,nwt
        ! Isospin is neutron in the first half of blocks, proton in the rest
        it = 2
        if(wave.le.sum(HFBlocks(1:Blocks/2))) it = 1

        Inproduct = 0.0_dp
        do k=1,4          
                do i=1,mv
                       Inproduct = Inproduct + HFPsi(i,1,1,k,wave) *  & 
                       &  ( HFddPsi(i,1,1,k,1,1,wave) + &
                       &    HFddPsi(i,1,1,k,2,2,wave) + &
                       &    HFddPsi(i,1,1,k,3,3,wave))
                enddo
        enddo
        Kinetic(it)= Kinetic(it) + Occupations(wave)*Inproduct
    enddo

    Kinetic=-Kinetic * hbm * dv
    return
  end function CompKinetic
  
  function Skyrme_LO() result(LOTerms)
    !---------------------------------------------------------------------------
    ! Leading order terms in the Skyrme functional. 
    ! Everything is calculated in the BFH representation (neutron-proton) 
    ! 
    ! LOTerms
    ! (1,2)   rho^2
    ! (3,4)   s^2
    ! (5,6)   rho^(2+alpha)
    ! (7,8)   s^(2+alpha)
    !---------------------------------------------------------------------------
    use Constants
    
    real(KIND=dp) :: LOTerms(8)
    integer       :: it
    
    LOTerms = 0.0_dp
    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! rho^2 term                                                Time-even
    LOTerms(1)   =              B1 * sum(sum(Rho,2)**2)
    do it=1,2
      LOTerms(2) = LOTerms(2) + B2 * sum(Rho(:,it)**2)
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! rho^(2+alpha) term                                        Time-even
    LOTerms(5)   =              B7*sum(sum(Rho,2)**(2 + byt3))
    do it=1,2
      LOTerms(6) = LOTerms(6) + B8*sum(Rho(:,it)**2*sum(Rho,2)**(byt3))
    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! s^(2) term                                                Time-odd
!    LOTerms(3)   =              B10*sum(sum(vecs,3)**2)
!    do it=1,2
!      LOTerms(4) = LOTerms(4) + B11*sum(vecs(:,:,it)**2)
!    enddo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! s^(2 + alpha) term                                        Time-odd
!    LOTerms(7)   =              B12*sum(sum(vecs,3)**(2 + byt3))
!    do it=1,2
!      LOTerms(8) = LOTerms(8) + B13*sum(vecs(:,:,it)**2*sum(vecs,3)**(byt3))
!    enddo
!    
    ! Don't forget the volume element
    LOTerms = LOTerms * dv
  end function Skyrme_LO
  
  function Skyrme_NLO() result(NLOTerms)
    !---------------------------------------------------------------------------
    ! Next-to-leading order terms in the Skyrme functional. 
    ! Everything is calculated in the BFH representation (neutron-proton) 
    !---------------------------------------------------------------------------
    use Constants
    
    real(KIND=dp) :: NLOTerms(24)
    integer       :: it
    
    NLOTerms = 0.0_dp
        
    NLOTerms(5) =                   B3 * sum(sum(Rho,2) * sum(tau(:,1,1,:) &
    &                                                   + tau(:,2,2,:)     &
    &                                                   + tau(:,3,3,:),2))    
    do it=1,2
        NLOTerms(6) = NLOTerms(6) + B4 * sum(Rho(:,it)*(tau(:,1,1,it)  &
    &                                                  +tau(:,2,2,it)  &
    &                                                  +tau(:,3,3,it)))    
    enddo
    
    NLOterms(7)   = B5*sum(sum(Rho,2)*sum(lap_Rho,2))

    !B6 Terms
    do it=1,2
      NLOterms(8) = NLOterms(8) + B6*sum(Rho(:,it)*lap_Rho(:,it))
    enddo

    ! Don't forget the volume element
    NLOTerms = NLOTerms * dv
  end function Skyrme_NLO
  
  function Skyrme_N2LO() result(N2LOTerms)
    !---------------------------------------------------------------------------
    ! Next-to-newt-to-leading order terms in the Skyrme functional. 
    ! Everything is calculated in the BFH representation (neutron-proton) 
    !---------------------------------------------------------------------------
    use Constants
    
    real(KIND=dp) :: N2LOTerms(32)
    integer       :: it
    
    N2LOTerms = 0.0_dp

    print *, 'N2RHOQ',N2RHOQ

    !---------------------------------------------------------------------------
    ! Delta rho Delta rho                              T-even
    N2LOterms(1) = sum(sum(lap_rho,2)**2)                            * N2D2rho(1)
    do it=1,2
        N2LOterms(2) =  N2LOterms(2) + sum(Lap_Rho(:,it)**2)         * N2D2rho(2) 
    enddo
    !----------------------------------------------------------------------------
    !  rho Q                                           T-even 
    N2LOterms(3) = sum(sum(rho,2) * sum(QN2LO,2))                    * N2rhoQ(1)
    do it=1,2
        N2LOterms(4) = N2LOterms(4) + sum(rho(:,it) * QN2LO(:,it)) 
    enddo
    N2LOterms(4) = N2LOterms(4)                                      * N2rhoQ(2)
    !----------------------------------------------------------------------------
    !  tau^2                                           T-even
    N2LOterms(5) = sum(sum(tau(:,1,1,:) + tau(:,2,2,:)                  &
                 &                      + tau(:,3,3,:),2)**2)         * N2tau(1)
    do it=1,2
        N2LOterms(6) = N2LOterms(6) + sum((tau(:,1,1,it) + tau(:,2,2,it) &
                    &               + tau(:,3,3,it))**2)              * N2tau(2)
    enddo
    
    ! Don't forget the volume element
    N2LOTerms = N2LOTerms * dv
  
  end function Skyrme_N2LO
 
end module functional
