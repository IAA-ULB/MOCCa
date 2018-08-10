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

 100 format (/,' ___________________________________________________________', &
     &       /,'|                                                          |', &
     &       /,'|                                                          |', &
     &       /,'|                                                          |', &
     &       /,'|  #######   ##   #    # #####   ##   #      #    #  ####  |', &
     &       /,'|     #     #  #  ##   #   #    #  #  #      #    # #      |', &
     &       /,'|     #    #    # # #  #   #   #    # #      #    #  ####  |', &
     &       /,'|     #    ###### #  # #   #   ###### #      #    #      # |', &
     &       /,'|     #    #    # #   ##   #   #    # #      #    # #    # |', &
     &       /,'|     #    #    # #    #   #   #    # ######  ####   ####  |', &
     &       /,'|                                                          |', &
     &       /,'|  Copyright  P.-H. Heenen, M.Bender & W. Ryssens          |', &
     &       /,'|__________________________________________________________|')
 
 299 format (  ' -------------- Version Information ------------------------')
 300 format (  '  VERSION1') ! Git commit
 301 format (  '  VERSION2') ! Author of commit
 302 format (  '  VERSION3') ! Date


 print 100
 print 299
 print 300
 print 301
 print 302
 
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
   
    integer :: iter
   
    !---------------------------------------------------------------------------
    ! Initial calculations
    !---------------------------------------------------------------------------
    ! Solve the pairing, with the current values of <h> and the pairing gaps.
    call SolvePairing()
    
    ! Derive all the single-particle wavefunctions
    call deriveall()
   
    ! Calculate the initial densities.
    call densit(SaveRho=.false.)

    call CalculateMoments()
    call calcFields()
    
    stop
    call CalcGaps(FermiEnergy)
    
    call CalcEnergy()

    ! Initial printout
    call printSpwfs
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
        call deriveall()
 
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
        
        !-----------------------------------------------------------------------
        ! Decide between full or partial printout.
        if(mod(iter,PrintIter).eq.0) then
            call PrintSpwfs
            call printallmoments
            call printpairing
            call PrintEnergy           
        else
            call printsummary(iter)
        endif
        !-----------------------------------------------------------------------
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
