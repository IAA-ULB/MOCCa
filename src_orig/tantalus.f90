module Tantalus

 use geninfo

 implicit none

contains

subroutine Run_Tantalus(run_mode, file_number,input_file)
 !==============================================================================
 !_________ _______  _       _________ _______  _                 _______ 
 !\__   __/(  ___  )( (    /|\__   __/(  ___  )( \      |\     /|(  ____ \
 !   ) (   | (   ) ||  \  ( |   ) (   | (   ) || (      | )   ( || (    \/
 !   | |   | (___) ||   \ | |   | |   | (___) || |      | |   | || (_____ 
 !   | |   |  ___  || (\ \) |   | |   |  ___  || |      | |   | |(_____  )
 !   | |   | (   ) || | \   |   | |   | (   ) || |      | |   | |      ) |
 !   | |   | )   ( || )  \  |   | |   | )   ( || (____/\| (___) |/\____) |
 !   )_(   |/     \||/    )_)   )_(   |/     \|(_______/(_______)\_______)
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
 use timing

 implicit none
 !------------------------------------------------------------------------------
 ! These inputs control where the code will look for its input. Leaving them 
 ! empty will have the code rely on STDIN for input.
 integer(dp), intent(in), optional   :: file_number   
 character(*), intent(in), optional  :: input_file 
 character(len=*), intent(in)        :: run_mode 
 character(len=43)                   :: mode_print
 character(len=26)                   :: symprint
 logical                             :: printed
 !------------------------------------------------------------------------------
 ! Information gleaned from git and the Makefile, to be used to identify the 
 ! executable
 character(len=57), parameter        :: version1 =VERSION1
 character(len=57), parameter        :: version2 =VERSION2
 character(len=57), parameter        :: version3 =VERSION3
 character(len=57), parameter        :: version4 =VERSION4
 character(len=57), parameter        :: compiler =COMPCOMP
 character(len=57), parameter        :: cflags   =CFLAGS
 character(len=57), parameter        :: optflags =OPTFLAGS

 !------------------------------------------------------------------------------
 ! MPI error code
#if(USE_MPI > 0)
 integer :: mpi_err
#endif

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
 300 format ( 8x,'| ', a57, '|') ! Git commit
 301 format ( 8x,'| ', a57, '|') ! Author of commit
 302 format ( 8x,'| ', a57, '|') ! Date
 303 format ( 8x,'| Branch: ', a49, '|') ! Branch
 304 format ( 8x,'|                                                          |')
 305 format ( 8x,'|-------------- Symmetry Information ----------------------|')
 306 format ( 8x,'| S.p. generators        = ', a26, 6x, '|')
 307 format ( 8x,'| Axis reduction  X Y Z  = ', 3i2, 26x, '|')
 308 format ( 8x,'| SYM_CODE               = ', a26, 6x, '|')
 309 format ( 8x,'| TRANS_CODE             = ', a26, 6x, '|')
 310 format ( 8x,'|-------------- Environment Information -------------------|')
 311 format ( 8x,'|  Number of MPI_ranks   = ', i6, 26x, '|')
 313 format ( 8x,'|-------------- Compilation Information -------------------|')
 314 format ( 8x,'| Compiled with:                                           |')
 315 format ( 8x,'| ', a57, '|')
 316 format ( 8x,'| Compilation flags reported:                              |')
 317 format ( 8x,'| ', a57, '|')
 318 format ( 8x,'| Optimisation flags reported:                             |')
 319 format ( 8x,'| ', a57, '|')
 320 format ( 8x,'|__________________________________________________________|')
 
 !------------------------------------------------------------------------------
 ! Start the different processes across MPI ranks and do MPI bookkeeping
#if(USE_MPI > 0) 
  call mpi_init(mpi_err)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, NCORES  , mpi_err)
  call MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, mpi_err)
  
  ! Set MPI errors to be fatal. This is the default setting, but it doesn't
  ! hurt to be verbose, precise and future-flexible.
  CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL,mpi_err)
#endif
 !------------------------------------------------------------------------------
 ! starting all timers 
 ! (disabled for now as I'm not sure how this interacts with MPI)
 call initialize_all_timers
 call start_timer(T_tantalus)

 !------------------------------------------------------------------------------
 ! Printing information to STDOUT on the run
 if(MPI_RANK .eq. 0) then
   print *
   print 100
   write(mode_print, '(a43)') run_mode
   print 200, adjustl(mode_print)
   print 299
   print 304
   print 300, version1
   print 301, version2
   print 302, version3
   print 303, version4
   print 304
   print 305
   print 304
   symprint = adjustl(SYMSTRING)
   print 306, symprint
   print 307, reduX, reduY, reduZ
   print 308, SYM_CODE
   print 309, TRANS_CODE
   print 310
   print 311, NCORES
   printed = .false.
   print 313
   print 314
   print 315, compiler
   print 316
   print 317, cflags
   print 318
   print 319, optflags
   print 320
 endif

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
 call ReachForWaterAndFood()
 !------------------------------------------------------------------------------
 ! Clean up after running, just in case we need to run again.
 call Cleanupthemess()
 !------------------------------------------------------------------------------
 ! end the processes across MPI ranks
#if(USE_MPI > 0) 
  call mpi_finalize(mpi_err)
#endif
 !------------------------------------------------------------------------------
 ! Print all timing info
 call stop_timer(T_tantalus)
! call print_all_timers()

 ! end of one mean-field calculation..;
end subroutine Run_Tantalus

subroutine ReachForWaterAndFood()
    !---------------------------------------------------------------------------
    ! Evolve the single-particle wavefunctions and densities.
    !
    ! 
    ! The overall iterative scheme is as explained in
    !   W. Ryssens, M. Bender, M. and P.-H. Heenen,
    !   Eur. Phys. J. A, 55, 93. https://doi.org/10.1140/epja/i2019-12766-6
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
    use cranking  
    use convergence
    use scfiteration
    use timing
    use fission_MOI

    implicit none

    1 format('----------------------------------')
    2 format('| Convergence criteria satisfied.|')
    3 format('| Needed ', i4, ' iterations.', 8x,'|')
    4 format('| dE    < ', es10.3, 12x, ' | ')
    5 format('| dQ2   < ', es10.3, 12x, ' | ')
    6 format('| d2H   < ', es10.3, 12x, ' | ')
   61 format('| |spg| < ', es10.3, 12x, ' | ')
    7 format('| dmu   < ', es10.3, 12x, ' | ')
   71 format('| dJz   < ', es10.3, 12x, ' | ')    
    8 format('| Ending the iterative proces.   |')

    9 format(' Iter =', i5, '; writing checkpoint to file ', a20, '.')
   10 format(86('-'))  
   11 format(30x, 'Iteration = ', i5, /)   
   12 format(24x, 'FINAL Iteration = ', i5, /)

    integer :: iter, iprint, scheme, ifail
    logical :: ConvergenceAchieved, calc_expensive
    ! Logical to see if any moments with projection are necessary
    logical :: projectpresent = .false.
    ! Message for the output of the code, useful for the Brussels group.
    character(len=99) :: iomsg = 'START'

    ifail = 0
    ConvergenceAchieved = .false.
    !---------------------------------------------------------------------------
    ! Initial calculations
    !---------------------------------------------------------------------------

    if( (Bogofromfile.and.readHFBinfofile) .and. pairingscheme.eq.1) then
      ! If using a gradient strategy and we want to continue from file.  
      ! Only allowed of course if we have actually read a Bogoliubov transfo.
      scheme = -1
    else
      scheme =  0
    endif
    call SolvePairing(scheme, ifail)
    if(ifail.ne.0) then
       ! Solve the pairing, with the current values of <h> and the pairing gaps.
       ! Note that this is ALWAYS a direct solve, i.e. we diagonalise the HFB 
       ! Hamiltonian with a LAPACK call. We do this if the code did not receive
       ! explicit instructions to start from the Bogoliubov transformation on 
       ! file. 
       print *, 'WARNING! Pairing solver failed.'
    endif
    if(bogofromfile .and. guessgaps .and. pairingscheme.eq.1) then
      ! We perform a few extra calls to solvepairing to take a few gradient
      ! steps, with finite values for Delta.
      call SolvePairing(pairingscheme, ifail)
    endif

    ! Construct the canonical basis    
    if(pairingtype.eq. 2) call ConstructCanonicalBasis()

    ! Derive all the single-particle wavefunctions in the HFPsi array
    call deriveHF()

    ! Calculate the initial densities and the charge density (separately)
    call densit(SaveRho=.false.)
    call ConstructChargeDensity(ChargeDensity)

    ! Adopt the relevant quantities to the centre-of-mass of the nucleus
    call adapt_com()

    call CalculateMoments()   !=> vital to be called here, 
                              !    (a) before the calculation of the fields
                              !    (b) after construction of the charge density
                              ! as
                              !  (a) the multipole cutoff is allocated in this
                              !      process, and is needed for the calculation
                              !      of the cranking fields
                              !  (b) the calculations of the charge rms radius
                              !      requires the charge density to be 
                              !      constructed

    ! Only calculate the fields that have not been read from either a 
    ! wavefunction file or a potential file. 
    call calcFields(calcall=.false.,precon= .false.)

    ! Update all spwf properties
    call update_spwf_properties( .true. ) ! expensive version

    call setBelyaevProcedure()
    call CalcEnergy(.true.)      ! Calculate the energy WITH all the expensive
                                 !   parts included. 
    call calc_avg_gap()

    ! Initial printout
    if(MPI_RANK .eq. 0) then
      ! only the very first MPI RANK prints all of this output
      call printSpwfs
      call printQps
      call printallmoments
      call print_boxsize_check
      call PrintMomentsofInertia
      call printcranking  
      call printpairing(pairstabfactor)
      call PrintEnergy 
    endif

    !---------------------------------------------------------------------------
    ! Start of the iterations
    !---------------------------------------------------------------------------
    do iter=1,maxiter
        call update_E_history()

        projectpresent   = checkconstraints() .or. check_cranking()    
        if(projectpresent) call feasibleproject()

        ! One evolution step for the spwfs
        call Evolve(iter)

        ! Calculate the gaps Delta with the current 
        ! a) fields 
        ! b) density matrix and anomalous density matrix 
        ! c) Fermi-energy
        PairStabfactor = CompStabilisingFactor(PairDenEnergy)
        call CalcGaps(FermiEnergy, PairStabFactor)

        ! Save Fermi energy
        FermiHistory   = FermiEnergy

        call SolvePairing(pairingscheme,ifail)
        if(pairingtype.eq. 2)  call ConstructCanonicalBasis()

        ! Derive all spwfs in the HF-basis
        call deriveHF()
        call densit(SaveRho=.true.)
        call ConstructChargeDensity(ChargeDensity)
        if(follow_com) call adapt_com()
        ! Calculate a) moments values, b) readjustment and c) finally their
        ! contribution to the sphamiltonian.
        call CalculateMoments()
        call ReadjustAllMoments(1)
        call ReadjustAllMoments(2)
        call Sphamilcontribution()

        ! Recalculate the fields, but only if MaxIter > FreezeIter
        if(iter .gt. freezeiter) then
          call calcFields(calcall=.true.,precon=.true.)
        endif

        ! Update all spwf properties
        call update_spwf_properties( .false. ) ! nonexpensive version

        call updateAM
        call ReadjustCranking
        !-----------------------------------------------------------------------
        ! Above: actual evolution of physical quantities
        ! Below: administration/bookkeeping
        !-----------------------------------------------------------------------
        !See if some moments were temporary
        call TurnOffConstraints(iter)

        ! Recalculate the energy
        if((mod(iter,PrintIter).eq.0) .or. (iter.eq.maxiter)) then
          iprint = 1
          calc_expensive = .true.
        else
          iprint = 0
          calc_expensive = .false.
        endif

        call CalcEnergy(calc_expensive)
        call calc_avg_gap()

        ! Check for convergence or a failed calculation
        if (ifail .ne. 0) then  
          iomsg               = 'FERMI'
          ConvergenceAchieved = .false.
          exit
        else  
          call Converged(ConvergenceAchieved, iter)  
        end if 

        if(convergenceAchieved) then
          iprint = 1
          ! Recalculate the energy with all parts included at the end 
          call CalcEnergy(.true.)
        endif
        !-----------------------------------------------------------------------
        ! Decide between full or partial printout.
        if(iprint .eq.1) then
            ! ... but update all spwf properties first to ensure correct prints
            call update_spwf_properties( .true. ) ! expensive version
            call updateAM 
            call ReadjustCranking

            if(MPI_RANK.eq.0) then
              print 10
              if((iter .eq. maxiter) .or. ConvergenceAchieved) then
                ! Add a clear indication this is the FINAL iteration
                print 12, iter  
              else
                print 11, iter
              endif
              call PrintSpwfs
              call PrintQps
              call printallmoments
              call print_boxsize_check
              call PrintMomentsofInertia
              call printcranking
              call printpairing(PairStabfactor)
              call printEnergy()
            endif
        elseif(MPI_RANK.eq.0) then
             ! ..... else print a summary
            call printsummary(iter)
        endif

        !-----------------------------------------------------------------------
        ! Write a wavefunction file according to checkpointiter
        if(checkpointiter.ne.0) then
          if(mod(iter,checkpointiter) .eq. 0) then
            if(MPI_RANK.eq.0) print 9, iter, outputfilename
            iomsg='CHECKPOINT'
            call WriteTantalus(12, outputfilename)     
          endif
        endif
        !-----------------------------------------------------------------------
        if(ConvergenceAchieved) then
          if(MPI_RANK .eq. 0) then
            print 1
            print 2
            print 3, iter
            print 4, energy_prec
            print 5, moment_prec
            print 6, disp_prec
            print 61, gradient_prec
            print 7, fermi_prec
            print 71, angmom_prec
            print 8
            print 1
          endif
          iomsg='CONVERGED'
          exit
        endif
    enddo
!    if(inversetemp .ne. -1) then
!        call projectThermal
!    endif    
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    ! Calculate and print the collective moment of inertias    
    if(N_inertia .gt. 0) then
      call calc_collective_inertia
      if(MPI_RANK.eq.0) then
        call print_collective_inertia
      endif
    endif
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
    if(iter.eq.maxiter+1) then
      iomsg='MAXITER'  
    endif
    !---------------------------------------------------------------------------
    ! Write output to the outputfile, i.e. the full wavefunction file
    call WriteTantalus(12, outputfilename)     
    !---------------------------------------------------------------------------
    ! Write other, advanced, output
    call write_advanced_output(iter-1,iomsg)
    !---------------------------------------------------------------------------
end subroutine ReachForWaterAndFood

subroutine printsummary(iter)
    !---------------------------------------------------------------------------
    ! Short printout after an iteration.
    ! 
    !---------------------------------------------------------------------------
    use functional
    use evolution 
    use moments    
    use pairing
    
    implicit none

    integer, intent(in)   :: iter
    type(Moment), pointer :: current, part
    real(KIND=dp)         :: dF(2), DN(2), dQ, dL, dev, val, devJ
    character(len=1)      :: t, spec

    1 format (86('-'))
    2 format (' Iteration = ',i4)
   21 format (' Potentials frozen.')
    3 format (' dt    = ', f8.4, 4x, '  mu   = ', f8.4, ' gradn = ', es12.3, ' D2H  = ', es12.3)
   31 format (' dtg   = ', f8.4, 4x, '  mug  = ', f8.4, ' gradn = ', es12.3)
    4 format (' E     = ', f10.3,2x, '  DE   = ', e12.5)
   41 format (' R     = ', f10.3,2x, '  DR   = ', e12.5)  
   42 format (' R-E   = ', f10.3,2x, 'D(R-E) = ', e12.5)  

    5 format (' ',a1, 'Q', 2i1,a1,' = ',f12.4, 3x, 'dQ = ', es8.1, 2x,         &
    &          'L = ',f12.4,2x,' dL = ', es8.1, 2x, 'dev = ', es8.1)

    6 format (' dmun  = ', es8.1, 4x, '  dmup= ', es8.1)
    7 format (' dN    = ', es8.1, 4x, '  dZ  = ', es8.1)  
    8 format (' Jz    = ', f12.4,  2x, ' dJZ= ', es8.1,  &
    &         ' Om = ' , f12.4,  2x ' dO = ', e8.1, 2x, 'dev = ', es8.1)

    part=>FindMoment(0,0,.false.)

    if(iter.eq.1) print 1
    print 2, iter
    if(freezeiter .gt. iter) print 21
    print 3, dt, momentum, gradientnorm, d2h
    if(pairingscheme.eq.1) then
      print 31, gradient_stepsize, gradient_mu, sqrt(sum(HFBGradnorm**2))
    endif
    print 4, totalE,     (totalE - Ehistory(1))/abs(totalE)
    print 41, Routhian,  (Routhian - Rhistory(1))/abs(Routhian)
    print 42, Routhian-totalE, &
    &  ((Routhian - Rhistory(1)) - (totalE - Ehistory(1)))/abs(totalE) 
    if(fixfermi) then
        dN = part%value - part%history 
        print 7, dN
    else
        dF   = FermiEnergy - FermiHistory
        print 6, dF
    endif
    
    Current => Root
    do while (associated(Current%next))
      Current => Current%next
      
      if((Current%constrainttype .ne. 0) .or. (Current%l .eq. 2)) then
        dL = Current%multiplier - Current%mult_hist
        if(Current%Impart) then
          t = 'I'
        else
          t = 'R'
        endif
        
        if(Current%constrainttype.ne.0) then
          dev = Current%deviation
        else
          dev = 0
        endif
        
        if(current%l.eq.2 .and. current%isoswitch .ne. 0) then
          dQ = sum(Current%value) - sum(Current%history)
          print 5, t, Current%l, Current%m, 't', sum(Current%value), dQ, &
          &              0.0d0, 0.0d0, 0.0d0
        endif
                
        select case(Current%isoswitch)
        case(0)
          dQ   = sum(Current%value) - sum(Current%history)
          spec = 't'
          val  = sum(Current%value)
        case(1,2)
          dQ = Current%value(Current%isoswitch)&
           & - Current%history(Current%isoswitch)     
          val =  Current%value(Current%isoswitch)   
          select case(Current%isoswitch)
          case(1)
            spec = 'n'
          case(2)
            spec = 'p'
          end select
        end select        
        print 5, t, Current%l, Current%m, spec, val, &
          &         dQ, Current%multiplier, dL, dev
      endif
    enddo    
    
    if(cranktype(3) .eq. 1) then
      devJ = TotalAngMom(3) - CrankValues(3)
    else
      devJ = 0.0
    endif
    print 8, totalangmom(3), totalangmom(3) - angmomold(3), &
    &        omega(3), omega(3)-omega_prev(3), devJ
    print 1
        
end subroutine printsummary

subroutine initialize_all_timers()
   !----------------------------------------------------------------------------
   ! Initialize all the timers that have been defined.
   ! Input: 
   !       NONE
   ! Output:
   !       NONE
   !----------------------------------------------------------------------------
   use timing
   
   call add_timer('Tantalus'                   , T_tantalus)  
   call add_timer('HF-basis Derivatives'       , T_derivatives)  
   call add_timer('Canonical basis Derivatives', T_derivatives_can)  
   call add_timer('Spwf evolution'             , T_evolution)  
   call add_timer('Orthonormalization'         , T_ortho)  
   call add_timer('Density calculations'       , T_densities)  
   call add_timer('Density: pp'                , T_den_pp)
   call add_timer('Density: ph'                , T_den_ph)  
   call add_timer('Density: derivatives'       , T_den_der)  
   call add_timer('Field calculations'         , T_fields)  
   call add_timer('Energy calculations'        , T_energy)  
   call add_timer('Pairing solver '            , T_pairing)  
   call add_timer('Sp. Hamiltonian '           , T_sphamil)  
   call add_timer('Coulomb solver'             , T_coulomb)  
   call add_timer('Can. basis construction'    , T_den_can)  
   call add_timer('Moments of inertia '        , T_MOI)  
   call add_timer('Centre-of-mass correction ' , T_COM)  
   call add_timer('COM one-body '              , T_COM1)  
   call add_timer('COM two-body '              , T_COM2)  
   call add_timer('Pairing gaps '              , T_gaps)  
   call add_timer('Multipole moments '         , T_moments)  
   call add_timer('Feas. Proj. step '          , T_feasible)  
   call add_timer('Spwf angular momentum '     , T_spwfangmom)  
   call add_timer('Charge density folding'     , T_chargedensity)  
   call add_timer('Collective MOIs'            , T_collective_moi)  
   call add_timer('Microscopic pairing'        , T_microscopic_pairing)  
   call add_timer('Orthogonalisation of h|psi>', T_Hortho)  
   call add_timer('Construction HF transfo'    , T_HFDiag)  

end subroutine initialize_all_timers

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
  use functional

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
  call clean_potentials

end subroutine cleanupthemess
end module Tantalus
