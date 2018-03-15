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
    
    implicit none
   
    integer :: iter, i
   
    call deriveall()
    !---------------------------------------------------------------------------
    ! Solve the pairing. For now just HF.
    call NaiveFill(occupations)
    !---------------------------------------------------------------------------
    ! Calculate the initial densities.
    call densit(0)
    ! And the initial fields.
    call calcFields()
    ! And even the initial energy.
    call CalcEnergy()

    call printSpwfs
    call PrintEnergy 
    
    do iter=1,maxiter
        ! One step in the evolution.
        call Evolve(iter)
        ! Restore all the different derivatives.
        call deriveall()
        ! Solve the pairing problem.
        call NaiveFill(occupations)
        ! Update the densities
        call densit(iter)
        ! Update the fields
        call calcFields()
        ! Recalculate the energy
        call CalcEnergy()
        ! Decide between full or partial printout.
        if(mod(iter,PrintIter).eq.0) then
            call PrintSpwfs
            call PrintEnergy            
        else
            call printsummary
        endif
    enddo
      
    !---------------------------------------------------------------------------
    ! Write output to the outputfile.
    call WriteTantalus(12, outputfilename)
    !---------------------------------------------------------------------------
end subroutine ReachForWaterAndFood


subroutine printsummary
    !---------------------------------------------------------------------------
    ! Not very advanced printing of a summary of the iteration.
    ! 
    !---------------------------------------------------------------------------
        print *,  '*************************************'
        print *,  ' Iteration ', iter
        print *,  ' Energy =  ', totalE
        print *,  ' Spwfs  =  ', spwfenergy
        print *,  ' GradNorm =  ', gradientnorm
        print *,  '*************************************'
        
end subroutine printsummary
