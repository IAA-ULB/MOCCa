module Tantalus

 implicit none

contains

subroutine Run_Tantalus(run_mode, file_number,input_file)
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
 !------------------------------------------------------------------------------
 ! While I (W.R.) like to think about Tantalus as a standalone code, this 
 ! 'main' routine is now written as a subroutine (with inputs!) to accomodate
 ! meta-codes that want to run Tantalus multiple times. 
 !==============================================================================

 use compilation
 use geninfo
 use wavefunctions
 use IO
 use temperature_projection 

 implicit none

 ! These inputs control where the code will look for its input. Leaving them 
 ! empty will have the code rely on STDIN for input.
 integer(dp), intent(in), optional   :: file_number   
 character(11), intent(in), optional :: input_file 
 character(len=*), intent(in)        :: run_mode 
 character(len=43)                   :: mode_print
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
     &   /,8x,'|  Copyright  P.-H. Heenen, M. Bender & W. Ryssens         |', &
     &   /,8x,'|                                                          |')
 
 200 format ( 8x, '|', 58('-'), '|'  ,/,8x, '| Runtype = ', a43, 4x, '|')

 299 format ( 8x,'|-------------- Version Information -----------------------|')
 300 format ( 8x,'| VERSION1',  8x, '|') ! Git commit
 301 format ( 8x,'| VERSION2',  7x, '|') ! Author of commit
 302 format ( 8x,'| VERSION3', 17x, '|') ! Date
 303 format ( 8x,'|__________________________________________________________|')


 print *
 print 100
 write(mode_print, '(a43)') run_mode
 print 200, adjustl(mode_print)
 print 299
 print 300
 print 301
 print 302
 print 303

 !------------------------------------------------------------------------------
 ! Read input from STDIN
 call ReadInput(file_number, input_file)
 !------------------------------------------------------------------------------
 ! Initalize relevant matrices throughout the code.
 call inilag() ! Derivative matrices. 
 !------------------------------------------------------------------------------
 ! Read wavefunctions
 call ReadWavefunction()
 !------------------------------------------------------------------------------
 ! Print all relevant input gleaned from STDIN and the wf file.
 call PrintInput(file_number, input_file)
 !------------------------------------------------------------------------------
 ! Go out and try to reach convergence, only to fail time and time again....
 call ReachForWaterAndFood
 !------------------------------------------------------------------------------
 ! Clean up after running, just in case we need to run again.
 call Cleanupthemess()
 
end subroutine Run_Tantalus

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
    use evolution

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
    ! Check all of the multipole moments that are large enough
    Current => Root

    do while(associated(Current%next)) 
        Current => Current%next
        if(Current%Beta(3).gt.0.05) then
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
    use temperature_projection
    use momentsofinertia    

    implicit none

    1 format('----------------------------------')
    2 format('| Convergence criteria satisfied.|')
    3 format('| Needed ', i4, ' iterations.', 8x,'|')
    4 format('| dE < ', e10.3, 15x, ' | ')
    5 format('| dQ < ', e10.3, 15x, ' | ')
    6 format('| dH < ', e10.3, 15x, ' | ')        
    7 format('| Ending the iterative proces.   |')
   
    integer :: iter
    logical :: ConvergenceAchieved
    ! Logical to see if any moments with projection are necessary
    logical :: projectpresent = .false.
    ! Message for the output of the code, useful for the Brussels group.
    character(len=99) :: iomsg = 'START'

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
    ! Construct the charge density on the (nx/ny/nz)-sized mesh.
    ! This was previously included in the Coulomb routines, but now needs to be 
    ! called separatedly, since the chargedensity is used for the calculation 
    ! of the rms radius.
    call ConstructChargeDensity(ChargeDensity)
    
    call CalculateMoments()

    ! Only calculate the fields that have not been initialized from file.
    call calcFields(calcall=.false.)
    
    call CalcGaps(FermiEnergy)
    call setBelyaevProcedure()
    call CalcEnergy()

    ! Initial printout
    call printSpwfs
    call printQps
    call printallmoments
    call PrintMomentsofInertia

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
       
        ! Save      
        FermiHistory = FermiEnergy
        ! Solve pairing problem
        call SolvePairing()

        projectpresent = checkconstraints()        

        if(projectpresent) then
          ! Update the densities
          ! Note that this update is incorrect, as we do not want to perform a 
          ! set of derivatives
          call densit(SaveRho=.true.)
          call ConstructChargeDensity(ChargeDensity)
          
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
        call ConstructChargeDensity(ChargeDensity)

        !See if some moments were temporary
        call TurnOffConstraints(iter)

        ! Calculate a) moments values, b) readjustment and c) finally their
        ! contribution to the sphamiltonian.
        call CalculateMoments()
        call ReadjustAllMoments(1)
        call Sphamilcontribution()
        
        ! Update all of the fields
        call calcFields(calcall=.true.)
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
            call PrintMomentsofInertia
            call printpairing
            call PrintEnergy           
        else
            call printsummary(iter)
        endif
        !-----------------------------------------------------------------------
        if(ConvergenceAchieved) then
            print 1
            print 2
            print 3, iter
            print 4, energy_prec
            print 5, moment_prec
            print 6, disp_prec
            print 7
            print 1

            iomsg='CONVERGED'
            exit
        endif
    enddo
    if(inversetemp .ne. -1) then
        call projectThermal
    endif    

    if(iter.eq.maxiter+1) then
      iomsg='MAXITER'  
    endif
    !---------------------------------------------------------------------------
    ! Write output to the outputfile.
    call WriteTantalus(12, outputfilename, iter-1, iomsg)     
    !---------------------------------------------------------------------------
end subroutine ReachForWaterAndFood

subroutine printsummary(iter)
    !---------------------------------------------------------------------------
    ! Short printout after an iteration
    ! 
    !---------------------------------------------------------------------------
    use functional
    use evolution 
    use moments    
    use pairing
    
    implicit none

    integer, intent(in)   :: iter
    type(Moment), pointer :: Q20, Q22, part
    real(KIND=dp)         :: dQ20, dQ22, dF(2), DN(2)
     

    1 format (80('-'))
    2 format (' Iteration = ',i4)
    3 format (' dt  = ', f8.4, 4x, '  mu  = ', f8.4, '  D2H = ', e8.1)
    4 format (' E   = ', f10.3,2x, '  DE  = ', e8.1)
    5 format (' Q20 = ', f12.4,    '  Q22 = ', f12.4, &
    &         ' dQ20= ', e8.1, 4x, '  dQ22= ', e8.1)
    6 format (' dmun= ', e8.1, 4x, '  dmup= ', e8.1)
    7 format (' dN  = ', e8.1, 4x, '  dZ  = ', e8.1)  

    part=>FindMoment(0,0,.false.)
    Q20 =>FindMoment(2,0,.false.     )
    Q22 =>FindMoment(2,2,.false., Q20)

    print 1
    print 2, iter
    print 3, dt, momentum, d2h
    print 4, totalE,  abs(totalE - Ehistory(1))/abs(totalE)

    if(fixfermi) then
        dN = part%value - part%history 
        print 7, dN
    else
        dF   = FermiEnergy - FermiHistory
        print 6, dF
    endif
    
    dQ20 = abs(sum(Q20%value) - sum(Q20%history))
    dQ22 = abs(sum(Q22%value) - sum(Q22%history))
    print 5, sum(Q20%value), sum(Q22%value), dQ20,dQ22
        
end subroutine printsummary

subroutine cleanupthemess()
  !-----------------------------------------------------------------------------
  !
  ! 
  !-----------------------------------------------------------------------------
  use geninfo
  use derivatives
  use wavefunctions
  use pairingcutoffs
  use BCS
  use HFB
  use pairing
  use densities
  use moments
  use coulombmod
  use evolution

  call clean_geninfo
  call clean_derivatives
  call clean_wavefunctions
  call clean_pairingcutoffs
  call clean_BCS
  call clean_HFB
  call clean_pairing
  call clean_densities
  call clean_moments
  call clean_coulomb
  call clean_evolution

end subroutine cleanupthemess
end module Tantalus
