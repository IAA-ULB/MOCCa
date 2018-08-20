program Tantalus
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
 use geninfo
 use wavefunctions
 use IO
 
 implicit none

 100 format &
     &  (/,8x,' ___________________________________________________________', &
     &   /,8x,'|                                                          |', &
     &   /,8x,'| MOCCa v2.0 =                                             |', &
     &   /,8x,'|                                                          |', &
     &   /,8x,'|  #######   ##   #    # #####   ##   #      #    #  ####  |', &
     &   /,8x,'|     #     #  #  ##   #   #    #  #  #      #    # #      |', &
     &   /,8x,'|     #    #    # # #  #   #   #    # #      #    #  ####  |', &
     &   /,8x,'|     #    ###### #  # #   #   ###### #      #    #      # |', &
     &   /,8x,'|     #    #    # #   ##   #   #    # #      #    # #    # |', &
     &   /,8x,'|     #    #    # #    #   #   #    # ######  ####   ####  |', &
     &   /,8x,'|                                                          |', &
     &   /,8x,'|  Copyright  P.-H. Heenen, M.Bender & W. Ryssens          |', &
     &   /,8x,'|                                                          |')
 
 299 format ( 8x,'|-------------- Version Information -----------------------|')
 300 format ( 8x,'| VERSION1',  8x, '|') ! Git commit
 301 format ( 8x,'| VERSION2',  7x, '|') ! Author of commit
 302 format ( 8x,'| VERSION3', 17x, '|') ! Date
 303 format ( 8x,'|__________________________________________________________|')


 print *
 print 100
 print 299
 print 300
 print 301
 print 302
 print 303
 
 !------------------------------------------------------------------------------
 ! Read input from STDIN
 call ReadInput
 !------------------------------------------------------------------------------
 ! Initalize relevant matrices throughout the code.
 call inilag() ! Derivative matrices. 
 !------------------------------------------------------------------------------
 ! Read wavefunctions
 call ReadWavefunction
 !------------------------------------------------------------------------------
 ! Print all relevant input gleaned from STDIN and the wf file.
 call PrintInput
 !------------------------------------------------------------------------------
 ! Go out and try to reach convergence, only to fail time and time again....
 call ReachForWaterAndFood
 
end program Tantalus

subroutine Converged(C) 
    !---------------------------------------------------------------------------
    ! Checks if the code has converged using the following convergence 
    ! criteria. 
    !
    !   energy_prec     1d-9     abs((E^(i) - E^(i-1))/E^(i))     < energy_prec 
    !                            This needs to be true across 5 iterations. 
    !
    !   moment_prec     1d-3     abs((Qlm^(i) - Qlm^(i))/Qlm^(i)) < moment_prec
    !                                if Qlm^(i) is large enough
    !   disp_prec       1d-5     abs(sum_i v^2_i <psi|h^2|psi> - epsilon^2)
    !                                     < disp_prec
    !---------------------------------------------------------------------------

    use Moments
    use functional

    logical       :: C
    integer       :: i
    real(KIND=dp) :: dE(5), dQ

    type(Moment), pointer  :: Current 

    C = .true.

    !---------------------------------------------------------------------------
    ! Checking the evolution of the energy
    do i =1,4
            dE(i) = abs(Ehistory(i) - Ehistory(i+1))/abs(totalE)
    enddo
    dE(5) = abs(TotalE - Ehistory(1))/abs(totalE)
    
    if(.not. all(dE .lt. energy_prec)) then
     C = .false.
    endif
    !---------------------------------------------------------------------------
    ! Checking the weighted dispersion
    if(d2H .gt. disp_prec) C = .false.

    !---------------------------------------------------------------------------
    ! Check all of the multipole moments
    Current => Root

    do while(associated(Current%next)) 
        Current => Current%next
        if(sum(Current%Value)>1) then
            dQ = abs(sum(Current%history)-sum(Current%value))
            dQ = dQ/abs(sum(Current%value))
            if(dQ > moment_prec) C = .false.
        endif
    enddo   

end subroutine Converged

subroutine ReachForWaterAndFood
    !---------------------------------------------------------------------------
    ! Evolve the single-particle wavefunctions and densities.
    !
    ! 
    ! The overall iterative scheme is as explained in
    !   W. Ryssens et al. [HEAVY-BALL PAPER]
    !
    !
    ! Which is
    !
    !   Initialization
    ! 
    !   Until convergence do
    !   |  1. Calculate matrix elements of h and Delta
    !   |  2. Evolve the HF-basis with the heavy-ball method
    !   |  3. Solve the pairing equations with matrix elements from 1.
    !   |     HF : fill the lowest levels
    !   |     BCS: solve the BCS equations to obtain the occupations
    !   |     HFB: a. solve the HFB equations in the HF-basis
    !   |          b. construct the canonical basis
    !   |  4. Perform feasible projection if asked for
    !   |  5. Construct the densities
    !   |    5b. Update the Lagrange multipliers of the constraints
    !   |  6. Construct the fields
    !   |     (including Coulomb and potential constraint contribution)
    !   |  7. Print iteration info
    !   |_____________________________
    !---------------------------------------------------------------------------
    use compilation
    use derivatives
    use wavefunctions
    use constants
    use densities
    use functional
    use evolution
    use IO
    use moments
    use coulombmod
    use pairing 
    use printing
    
    implicit none

    1 format('----------------------------------')
    2 format('| Convergence criteria satisfied.|')
    3 format('| dE < ', e10.3, ' | ')
    4 format('| dQ < ', e10.3, ' | ')
    5 format('| dH < ', e10.3, ' | ')        
    6 format('| Ending the iterative proces.   |')
   
    integer :: iter
    logical :: ConvergenceAchieved

    ConvergenceAchieved = .false.   
    
    !---------------------------------------------------------------------------
    ! Initial calculations
    !---------------------------------------------------------------------------
    ! Solve the pairing, with the current values of <h> and the pairing gaps.
    call SolvePairing()
    
    ! Derive all the single-particle wavefunctions
    call deriveHF()
   
    ! Calculate the initial densities.
    call densit(SaveRho=.false.)

    call CalculateMoments()
    call calcFields()
    call CalcGaps(FermiEnergy)
    call CalcEnergy()

    ! Initial printout
    call printSpwfs
    call printQps
    call printallmoments
    call printpairing
    call PrintEnergy 
  
    !---------------------------------------------------------------------------
    ! Start of the iterations
    !---------------------------------------------------------------------------
    do iter=1,maxiter
    
        ! Calculate the gaps Delta with the current 
        ! a) fields 
        ! b) density matrix and anomalous density matrix 
        ! c) Fermi-energy
        call CalcGaps(FermiEnergy)
        
        ! One heavy-ball step.
        ! Note that the (diagonal) matrix elements of <h> get calculated here
        call Evolve(iter)
       
        ! Solve pairing problem
        call SolvePairing()
        
        if(projectpresent) then
          ! Update the densities
          ! Note that this update is incorrect, as we do not want to perform a 
          ! set of derivatives
          call densit(SaveRho=.true.)
          
          call CalculateMoments()
          ! Readjust the projection constraints here, to not take into account
          ! the update from the projection
          call ReadjustAllMoments(2)
          
          ! Do an approximate projection on the feasible set
          call feasibleproject()

          ! Solve the pairing problem.
          call SolvePairing()
        endif
 
        ! Restore all the different derivatives.
        call deriveHF()
 
        ! Update the densities
        if(projectpresent) then 
          call densit(SaveRho=.false.)
        else
          call densit(SaveRho=.true.)
        endif
        ! Calculate a) moments values, b) readjustment and c) finally the 
        ! contribution to the sphamiltonian.
        call CalculateMoments()
        call ReadjustAllMoments(1)
        call Sphamilcontribution()
        
        ! Update the fields
        call calcFields()
        ! Recalculate the energy
        call CalcEnergy()

        ! Check for convergence
        call Converged(ConvergenceAchieved)
        
        !-----------------------------------------------------------------------
        ! Decide between full or partial printout.
        if(mod(iter,PrintIter).eq.0 .or. ConvergenceAchieved) then
            call PrintSpwfs
            call PrintQps
            call printallmoments
            call printpairing
            call PrintEnergy           
        else
            call printsummary(iter)
        endif
        !-----------------------------------------------------------------------
        if(ConvergenceAchieved) then
            print 1
            print 2
            print 3, energy_prec
            print 4, moment_prec
            print 5, disp_prec
            print 6
            print 1
            exit
        endif
    enddo
      
    !---------------------------------------------------------------------------
    ! Write output to the outputfile.
    call WriteTantalus(12, outputfilename)
    !---------------------------------------------------------------------------
end subroutine ReachForWaterAndFood


subroutine printsummary(iter)
    !---------------------------------------------------------------------------
    ! Not very advanced printing of a summary of the iteration.
    ! 
    !---------------------------------------------------------------------------
    use functional
    use evolution 
    integer, intent(in) :: iter

    print *,  '*************************************'
    print *,  ' Iteration ', iter
    print *,  ' dt=', dt, ' mu=', momentum
    print *,  ' Energy =  ', totalE, ' Spwfs  =  ', spwfenergy
    print *,  ' GradNorm =', gradientnorm, ' Dispersion=', d2h
    print *,  '*************************************'
        
end subroutine printsummary
