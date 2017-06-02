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
 
  
        real(KIND=dp) :: t0=-2488.913 
        real(KIND=dp) :: x0=0.834
        real(KIND=dp) :: t1=486.818
        real(KIND=dp) :: x1=-0.344
        real(KIND=dp) :: t2=-546.395
        real(KIND=dp) :: x2=-1.0
        real(KIND=dp) :: t3=13777.0
        real(KIND=dp) :: x3=1.354
        real(KIND=dp) :: yt3a=0.166666666666666666667 
        real(KIND=dp) :: te=0.0
        real(KIND=dp) :: to=0.0
        real(KIND=dp) :: wso=123.0
        real(KIND=dp) :: wsoq=123.0

        real(KIND=dp) :: t1n2=24.3409
        real(KIND=dp) :: t2n2=-27.31975
        real(KIND=dp) :: x1n2=-0.344
        real(KIND=dp) :: x2n2=-1.0   
         
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 ! Energies are calculated in the BFH representation using the B coupling 
 ! coefficients, and then recombined into the isospin representation using 
 ! the C coefficients. By default only the latter is printed, but the BFH
 ! representation can be asked for for debugging purposes. 
 ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
 
 real(KIND=dp) :: Kinetic(2), Skyrme(2)
 
 ! Declaration of the energy terms and the coupling coefficients 
$DECLARATION
 

 contains
 
 subroutine calcedfcoefs()
     
$CALCCOEF
    
 end subroutine calcedfcoefs
 
 subroutine printedfcoefs
    !---------------------------------------------------------------------------
    ! Print the values of the EFD coefs used.
   
    1 format ('-----------------------------------------------------------------')
    2 format (' Skyrme coupling constants ')
    3 format (26x, 'Isospin representation BFH representation')
    4 format (26x, 'scalar      vector      total       (n-p)')
    
    print 1
    print 2
    print 3
    print 4
    print 1
$PRINTCOEF    
    print 1
 end subroutine printedfcoefs
 
 subroutine CompSkyrme()
    !---------------------------------------------------------------------------
    !
    ! 
    !
    !
    !
    integer       :: it,m,n,k
    real(KIND=dp) :: Edensity(mv,3)
    
$CALCULATION
    
    Skyrme = $TOTAL
    
 end subroutine CompSkyrme
 
 subroutine PrintSkyrme()
    !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ! Print all contributions to the energy.
    !
    
    1 format ('-----------------------------------------------------------------')
    2 format (' Skyrme Energy ')
    3 format (30x, 'isoscalar   isovector   total')
    4 format ('Total :' 3f15.6)
    
    
    print 1
    print 2
    print 3
    print 1
    
    $PRINT 

    print 1
    print 4, Skyrme, sum(Skyrme)
    print 1    
 end subroutine PrintSkyrme

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

  subroutine calcFields()
        
        integer :: m, n, k, it

        do it=1,2
$CALCFIELDS
        enddo
  end subroutine calcFields 

!  
!  function Skyrme_LO() result(LOTerms)
!    !---------------------------------------------------------------------------
!    ! Leading order terms in the Skyrme functional. 
!    ! Everything is calculated in the BFH representation (neutron-proton) 
!    ! 
!    ! LOTerms
!    ! (1,2)   rho^2
!    ! (3,4)   s^2
!    ! (5,6)   rho^(2+alpha)
!    ! (7,8)   s^(2+alpha)
!    !---------------------------------------------------------------------------
!    use Constants
!    
!    real(KIND=dp) :: LOTerms(8)
!    integer       :: it
!    
!    LOTerms = 0.0_dp
!    
!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!    ! rho^2 term                                                Time-even
!    LOTerms(1)   =              B1 * sum(sum(Rho,2)**2)
!    do it=1,2
!      LOTerms(2) = LOTerms(2) + B2 * sum(Rho(:,it)**2)
!    enddo
!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!    ! rho^(2+alpha) term                                        Time-even
!    LOTerms(5)   =              B7*sum(sum(Rho,2)**(2 + byt3))
!    do it=1,2
!      LOTerms(6) = LOTerms(6) + B8*sum(Rho(:,it)**2*sum(Rho,2)**(byt3))
!    enddo
!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!    ! s^(2) term                                                Time-odd
!!    LOTerms(3)   =              B10*sum(sum(vecs,3)**2)
!!    do it=1,2
!!      LOTerms(4) = LOTerms(4) + B11*sum(vecs(:,:,it)**2)
!!    enddo
!    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
!    ! s^(2 + alpha) term                                        Time-odd
!!    LOTerms(7)   =              B12*sum(sum(vecs,3)**(2 + byt3))
!!    do it=1,2
!!      LOTerms(8) = LOTerms(8) + B13*sum(vecs(:,:,it)**2*sum(vecs,3)**(byt3))
!!    enddo
!!    
!    ! Don't forget the volume element
!    LOTerms = LOTerms * dv
!  end function Skyrme_LO
!  
!  function Skyrme_NLO() result(NLOTerms)
!    !---------------------------------------------------------------------------
!    ! Next-to-leading order terms in the Skyrme functional. 
!    ! Everything is calculated in the BFH representation (neutron-proton) 
!    !---------------------------------------------------------------------------
!    use Constants
!    
!    real(KIND=dp) :: NLOTerms(24)
!    integer       :: it
!    
!    NLOTerms = 0.0_dp
!        
!    NLOTerms(5) =                   B3 * sum(sum(Rho,2) * sum(tau(:,1,1,:) &
!    &                                                   + tau(:,2,2,:)     &
!    &                                                   + tau(:,3,3,:),2))    
!    do it=1,2
!        NLOTerms(6) = NLOTerms(6) + B4 * sum(Rho(:,it)*(tau(:,1,1,it)  &
!    &                                                  +tau(:,2,2,it)  &
!    &                                                  +tau(:,3,3,it)))    
!    enddo
!    
!    NLOterms(7)   = B5*sum(sum(Rho,2)*sum(lap_Rho,2))

!    !B6 Terms
!    do it=1,2
!      NLOterms(8) = NLOterms(8) + B6*sum(Rho(:,it)*lap_Rho(:,it))
!    enddo

!    ! Don't forget the volume element
!    NLOTerms = NLOTerms * dv
!  end function Skyrme_NLO
!  
!  function Skyrme_N2LO() result(N2LOTerms)
!    !---------------------------------------------------------------------------
!    ! Next-to-newt-to-leading order terms in the Skyrme functional. 
!    ! Everything is calculated in the BFH representation (neutron-proton) 
!    !---------------------------------------------------------------------------
!    use Constants
!    
!    real(KIND=dp) :: N2LOTerms(32)
!    integer       :: it
!    
!    N2LOTerms = 0.0_dp

!    !---------------------------------------------------------------------------
!    ! Delta rho Delta rho                              T-even
!    N2LOterms(1) = sum(sum(lap_rho,2)**2)                            * N2D2rho(1)
!    do it=1,2
!        N2LOterms(2) =  N2LOterms(2) + sum(Lap_Rho(:,it)**2)         * N2D2rho(2) 
!    enddo
!    !----------------------------------------------------------------------------
!    !  rho Q                                           T-even 
!    N2LOterms(3) = sum(sum(rho,2) * sum(QN2LO,2))                    * N2rhoQ(1)
!    do it=1,2
!        N2LOterms(4) = N2LOterms(4) + sum(rho(:,it) * QN2LO(:,it)) 
!    enddo
!    N2LOterms(4) = N2LOterms(4)                                      * N2rhoQ(2)
!    !----------------------------------------------------------------------------
!    !  tau^2                                           T-even
!    N2LOterms(5) = sum(sum(tau(:,1,1,:) + tau(:,2,2,:)                  &
!                 &                      + tau(:,3,3,:),2)**2)         * N2tau(1)
!    do it=1,2
!        N2LOterms(6) = N2LOterms(6) + sum((tau(:,1,1,it) + tau(:,2,2,it) &
!                    &               + tau(:,3,3,it))**2)              * N2tau(2)
!    enddo
!    
!    ! Don't forget the volume element
!    N2LOTerms = N2LOTerms * dv
!  
!  end function Skyrme_N2LO
 
end module functional
