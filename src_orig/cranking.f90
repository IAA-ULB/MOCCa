module cranking
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

 use compilation
 use derivatives
 use wavefunctions
 use geninfo
 use nil8
 use pairing

 implicit none

 !------------------------------------------------------------------------------
 ! Omega:
 !   Cranking frequency or Lagrange multiplier of the angular momentum in 
 !   the three Cartesian directions.
 !------------------------------------------------------------------------------
 real(KIND=dp) :: Omega(3)      = 0.0_dp
 !------------------------------------------------------------------------------
 ! TotalAngMom:
 !    Total angular momentum in the three Cartesian directions, calculated
 !    by summation of the single-particle contributions.
 ! AngMomOld: 
 !    Values of the total angular momentum at the previous iteration, used for
 !    readjustment of the cranking constraints.
 ! J2:
 !    Values of the total angular momentum squared, <J_i^2>, for the three
 !    Cartesian directions. 
 ! AMBlock:
 !    Values of the total angular momentum, split by quantum number block.
 !------------------------------------------------------------------------------
 real(KIND=dp), public :: TotalAngMom(3)= 0.0_dp, AngMomOld(3)  = 0.0_dp
 real(KIND=dp), public :: J2(3)         = 0.0_dp, AMBlock(8,3)  = 0.0_dp

contains

  subroutine updateAM
    !---------------------------------------------------------------------------
    ! Calculate the total angular momentum and related observables.
    !---------------------------------------------------------------------------  
    integer :: B, N, wave, si

    si = 0    
    do B=1,8
      N = HFBlocks(B) ; if(N .eq. 0) cycle
      do wave = 1, N  
        TotalAngMom = TotalAngMom + rho_can(si+wave) * spwf_J(:,si+wave)
      enddo
      si = si + N
    enddo
    
  end subroutine updateAM

  subroutine PrintCranking
    !---------------------------------------------------------------------------
    ! Prints all kinds of information about the expectation value of the 
    ! angular momentum operator and all kinds of angles.
    !---------------------------------------------------------------------------

    1 format (2x,74('_') )
   10 format (2x,74('-'))
    2 format (25('-'), ' Angular Momentum (hbar) ',26('-') )
    3 format (15x, 'Spwfs(*)  ',2x, 'Desired', 5x, 'Omega', 7x, 'E (MeV)' 6x,'Densit. ')

    4 format (3x,'J_',a1,'   ','|', 5f12.5 )
   31 format (3x,'Size  |', 3f12.5,12x,1f12.5)
!   32 format (1x,'ReJT',a1,'   ','|', 5f12.5 )
!   33 format (1x,'ImJT',a1,'   ','|', 5f12.5 )
!   34 format (2x,'|J|' ,a1,'   ','|', 5f12.5 )

    print 2
    print *
    print 3
    print 1

    print 4, 'x', TotalAngMom(1), 0.0, Omega(1), 0.0, 0.0 
    print 4, 'y', TotalAngMom(2), 0.0, Omega(2), 0.0, 0.0 
    print 4, 'z', TotalAngMom(3), 0.0, Omega(3), 0.0, 0.0
    print 1
    print 31, sqrt(sum(totalangmom(1:3)**2)), 0.0, &
    &         sqrt(sum(omega(1:3)**2))      , 0.0
    print 10
  end subroutine PrintCranking


end module 
