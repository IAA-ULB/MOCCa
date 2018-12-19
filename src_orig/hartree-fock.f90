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
    !
    ! At zero temperature, the occupations are simply set to be one or zero.
    ! Note that this does not necessarily play well near spherical symmetry for
    ! non-magic numbers of nucleons.
    !---------------------------------------------------------------------------
    
    integer :: i,j,n,p, ProtonUpperBound, NeutronUpperBound
    integer :: ProtonOrder(nwp), NeutronOrder(nwn)
    real(KIND=dp), intent(out) :: occupations(nwt)
    
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

  subroutine FiniteTemperatureHF(occupations,  Fermi)
    !---------------------------------------------------------------------------
    ! Simple bisection routine to find the correct Fermi energy for a finite
    ! temperature HF calculation. A Newton type method might be more efficient, 
    ! but I fear for rounding errors and instability with exp() of either large
    ! or small numbers.
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(out) :: occupations(nwt), energies(nwt,2)
    real(KIND=dp)              :: Fermi(2), Fmin, Fplus, betaE, Nmin, Nmax
    integer                    :: Order(nwt,2), nw
    integer                    :: it, maxiter=100, i, N

    
    !---------------------------------------------------------------------------
    !Finding the order of the spwfs, in terms of energy
    Order = 0    
    Order(1:nwn,1) = OrderSpwfsISO(-1)
    Order(1:nwp,2) = OrderSpwfsISO(+1)
    
    do i=1,nwn
        energies(i,2) = spenergies(Order(i,1))
    enddo
    do i=1,nwp
        energies(i,2) = spenergies(Order(i,2)) 
    enddo      

    do it=1,2
        
        if(it .eq. 1) then
            N = neutrons ; nw = nwn
        else then
            N = protons  ; nw = nwp
        endif

        ! Establish a search interval
        Fmin = spenergies(Order( 1,it))
        Fmax = spenergies(Order(nw,it))
    
        Nmin = sumoccupations(Fermi(it), energies(1:nw,it))
        Nmax = 
    



    enddo    

  end subroutine FiniteTemperatureHF

  function sumoccupations(mu, energies) result (N)
    !---------------------------------------------------------------------------
    ! Simple function that sums the occupations for given mu and sp-energies.
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: energies(:), mu
    real(KIND=dp)             :: N
    integer                   :: i

    N = 0
    do i=1, size(energies)
        betaE = inversetemp * (energies(i) - mu)
        N = N + 1.0/(1 + exp(betaE))
    enddo

  end function sumoccupations
  

  recursive function FermiBisection()
  

  end function FermiBisection

  subroutine CalcHFgaps(Fermi)
    !---------------------------------------------------------------------------
    ! Dummy routine.
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Fermi(2)  

  end subroutine calcHFgaps

end module hartreefock
