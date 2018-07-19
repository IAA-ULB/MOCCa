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
     &       /,'|                                            V 0.1 (Pan)   |', &
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
    use hartreefock
    use functional
    use evolution
    use IO
    use moments
    use coulombmod
    
    implicit none
   
    integer :: iter
   
    ! Derive all the single-particle wavefunctions
    call deriveall()
   
    ! Solve the pairing. For now just HF.
    call NaiveFill(occupations)
   
    ! Calculate the initial densities.
    call densit()

    call DealWithMoments()
    
    call calcFields()
    call CalcEnergy()

    ! Initial printout
    call printSpwfs
    call printallmoments
    call PrintEnergy 

    !---------------------------------------------------------------------------
    ! Start of the iterations
    do iter=1,maxiter
        ! One step in the evolution.
        call Evolve(iter)
        ! Solve pairing problem the first time
        call NaiveFill(occupations)
        
        if(projectpresent) then
          ! Do an approximate projection on the feasible set
          call feasibleproject()
          ! Solve pairing problem the second time
          call NaiveFill(occupations)
        endif
        
        ! Restore all the different derivatives.
        call deriveall()
        ! Solve the pairing problem.
        call NaiveFill(occupations)
        
        ! Update the densities
        call densit()
        
        call DealWithMoments()
        
        ! Update the fields
        call calcFields()
        ! Recalculate the energy
        call CalcEnergy()
        
        !-----------------------------------------------------------------------
        ! Decide between full or partial printout.
        if(mod(iter,PrintIter).eq.0) then
            call PrintSpwfs
            call printallmoments
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
