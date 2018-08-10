module hartreefock
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
 ! Module containing all of the necessary routines to correctly decide on the 
 ! occupation of single-particle wave-functions in the case of a Hartree-Fock
 ! calculation.
 !
 !
 !==============================================================================
 
 use wavefunctions
 
 implicit none
 
contains
 
 subroutine NaiveFill(occupations)
    !---------------------------------------------------------------------------
    ! This subroutine finds the orbitals with the lowest single particle
    ! energy and fills them, after sorting all the levels.
    !---------------------------------------------------------------------------
    
    integer :: i,j,n,p, ProtonUpperBound, NeutronUpperBound
    integer :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp), intent(out) :: occupations(2*nwt)
    
    n=0; p=0
    !---------------------------------------------------------------------------
    !Setting all occupation numbers to Zero
    Occupations = 0.0_dp
    !---------------------------------------------------------------------------
    ProtonUpperBound  = floor(Protons)
    NeutronUpperBound = floor(Neutrons)
    !---------------------------------------------------------------------------
    !Finding the order of the spwfs, in terms of energy
    ProtonOrder = OrderSpwfsISO(+1)
    NeutronOrder= OrderSpwfsISO(-1)
    !---------------------------------------------------------------------------
    !Filling in the lowest Proton orbitals. This is easy, since we know
    !the order of Spwfs.
    i=1 
    do while(p.lt.ProtonUpperBound .and. i.le.nwp)
      j              = ProtonOrder(i)
      Occupations(j) = 2.0_dp
      p = p + 2
      i = i + 1
    enddo
    i=1
    do while(n.lt.NeutronUpperBound .and. i.le.nwn)
      j              = NeutronOrder(i)
      Occupations(j) = 2.0_dp
      n = n + 2
      i = i + 1
    enddo
    return
  end subroutine NaiveFill
  
  subroutine CalcHFgaps(Fermi)
    !---------------------------------------------------------------------------
    ! Dummy routine.
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Fermi(2)  
  end subroutine calcHFgaps

end module hartreefock
