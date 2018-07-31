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

 print 100
 
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
    !
    ! Evolve the single-particle wavefunctions and densities.
    !
    !
    !
    !
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
    
    implicit none
   
    integer :: iter
   
    ! Derive all the single-particle wavefunctions
    call deriveall()
   
    ! Solve the pairing.
    call SolvePairing()
    
    ! Calculate the initial densities.
    call densit(SaveRho=.false.)

    call CalculateMoments()
    
    call calcFields()
    call CalcEnergy()

    ! Initial printout
    call printSpwfs
    call printallmoments
    call printpairing
    call PrintEnergy 

    !---------------------------------------------------------------------------
    ! Start of the iterations
    do iter=1,maxiter
        ! One step in the evolution.
        call Evolve(iter)
       
        ! Solve pairing problem the first time
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
        endif
        
        ! Restore all the different derivatives.
        call deriveall()
        ! Solve the pairing problem.
        call SolvePairing()
        
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
