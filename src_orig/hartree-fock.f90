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

 real(KIND=dp) :: HFdispersion(2) = 0.0
 
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

    real(KIND=dp), intent(out) :: occupations(nwt)      
    real(KIND=dp)              :: energies(nwt,2), N
    real(KIND=dp)              :: Fermi(2), Fmin, Fmax, betaE, Nmin, Nmax
    integer                    :: Order(nwt,2), nw
    integer                    :: it, maxiter=100, i

    if( fixfermi ) then
        !-----------------------------------------------------------------------
        ! We perform a calculation at fixed chemical potential
        Fermi(1) = mun ; Fermi(2) = mup
    else
        !-----------------------------------------------------------------------
        ! We fix the particle number to <N> = neutrons, <Z> = protons.
    
        !Finding the order of the spwfs, in terms of energy
        Order = 0    
        Order(1:nwn,1) = OrderSpwfsISO(-1)
        Order(1:nwp,2) = OrderSpwfsISO(+1)
        
        do i=1,nwn
            energies(i,1) = spenergies(Order(i,1))
        enddo
        do i=1,nwp
            energies(i,2) = spenergies(Order(i,2)) 
        enddo      

        do it=1,2
            
            if(it .eq. 1) then
                N = neutrons ; nw = nwn
            else 
                N = protons  ; nw = nwp
            endif
            
            !-------------------------------------------------------------------
            ! Establish a search interval
            Fmin = -1000
            Fmax = 100

            Nmin = FToccupations(Fmin, energies(1:nw,it)) - N
            Nmax = FToccupations(Fmax, energies(1:nw,it)) - N
        
            if(Nmin .gt. 0 .or. Nmax .lt. 0) then
                print *, 'Bracketing of Fermi energy is wrong.'
                print *, Fmin, Nmin
                print *, Fmax, Nmax
                print *, spenergies(Order(1,it)),spenergies(Order(nw,it)) 
                stop
            endif

            Fermi(it) = FermiBisection(Fmin,Nmin,Fmax,Nmax, energies(1:nw,it),N)
        enddo       
    endif

    !---------------------------------------------------------------------------
    ! Actually calculate the occupations for the fixed chemical potential
    do i = 1,nwt
        if (i .gt. nwn) then
            it    = 2 
        else
            it    = 1
        endif        
        betaE = inversetemp * (spenergies(i) - Fermi(it))
        occupations(i) = 2.0/(1 + exp(betaE))
    enddo

    !---------------------------------------------------------------------------
    ! Calculate the HF dispersion
    !  DN = 2  Tr (rho (1-rho))
    !---------------------------------------------------------------------------
    HFdispersion = 0.0
    do i=1,nwt
        if (i .gt. nwn) then
            it    = 2 
        else
            it    = 1
        endif   
        ! Note the time-reversal factors of two
        HFdispersion(it) = HFdispersion(it) +                                  &
        &                            occupations(i)/2.0 * (1-occupations(i)/2.0)
    enddo
    ! Factor of two from the formula
    HFdispersion = 2 * HFDispersion
    ! Factor two for time-reversal
    HFdispersion = 2 * HFDispersion
  end subroutine FiniteTemperatureHF

  function FToccupations(mu, energies) result (N)
    !---------------------------------------------------------------------------
    ! Simple function that sums the occupations for given mu and sp-energies.
    !---------------------------------------------------------------------------

    real(KIND=dp), intent(in) :: energies(:), mu
    real(KIND=dp)             :: N, betaE
    integer                   :: i

    N = 0
    do i=1, size(energies)
        betaE = inversetemp * (energies(i) - mu)
        N     = N + 1.0/(1 + exp(betaE))
    enddo

    ! Time-reversal gives the factor 2
    N = 2 * N

  end function FToccupations
  
  recursive function FermiBisection(xmin, Nmin, xmax, Nmax, energies, N) result(x)
    !---------------------------------------------------------------------------
    ! Bisection search for an appropriate chemical potential for the FT HF.
    !---------------------------------------------------------------------------

    real(KIND=dp)             :: xmin, xmax, x, Nmin, Nmax, xnew, Nnew
    real(KIND=dp), intent(in) :: energies(:), N

    if(abs(Nmin) .lt. pairing_prec) then
        x = xmin        
        return 
    endif
    if(abs(Nmax) .lt. pairing_prec) then
        x = xmax
        return
    endif   
    
    xnew = 0.5 * (xmin + xmax)
    Nnew = FToccupations(xnew, energies) - N

    if(Nnew .gt. 0) then
        x = FermiBisection(xmin, Nmin, xnew, Nnew, energies,N)
    else
        x = FermiBisection(xnew, Nnew, xmax, Nmax, energies,N)            
    endif

  end function FermiBisection

  subroutine CalcHFgaps(Fermi)
    !---------------------------------------------------------------------------
    ! Dummy routine.
    !
    !---------------------------------------------------------------------------
    real(KIND=dp), intent(in) :: Fermi(2)  

  end subroutine calcHFgaps

end module hartreefock
