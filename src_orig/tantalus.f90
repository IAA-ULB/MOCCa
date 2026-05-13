module Tantalus

   use geninfo

   implicit none

contains

#if( $FAM == 0)
   subroutine Run_Tantalus(file_number, input_file)
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
      use timing
      use fission_MOI
      use version, only: print_header

      implicit none
      !------------------------------------------------------------------------------
      ! Convergence signals
      ! Iteration counter
      integer           :: iteration
      ! Message for the output of the code, useful for the Brussels group.
      character(len=99) :: iomsg = 'START'
      !------------------------------------------------------------------------------
      ! These inputs control where the code will look for its input. Leaving them
      ! empty will have the code rely on STDIN for input.
      integer(dp), intent(in), optional   :: file_number
      character(*), intent(in), optional  :: input_file
      !------------------------------------------------------------------------------
      ! MPI error code
#if(USE_MPI > 0)
      integer :: mpi_err
#endif

      !------------------------------------------------------------------------------
      ! Start the different processes across MPI ranks and do MPI bookkeeping
#if(USE_MPI > 0)
      call mpi_init(mpi_err)
      call MPI_COMM_SIZE(MPI_COMM_WORLD, NPROCS, mpi_err)
      call MPI_COMM_RANK(MPI_COMM_WORLD, MPI_RANK, mpi_err)

      ! Set MPI errors to be fatal. This is the default setting, but it doesn't
      ! hurt to be verbose, precise and future-flexible.
      CALL MPI_Comm_set_errhandler(MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL, mpi_err)
#endif
      !------------------------------------------------------------------------------
      ! starting all timers
      ! (disabled for now as I'm not sure how this interacts with MPI)
      call initialize_all_timers
      call start_timer(T_tantalus)

      !------------------------------------------------------------------------------
      ! Printing information to STDOUT on the run
      call print_header(.false.)
      !------------------------------------------------------------------------------
      ! Read input from STDIN
      call ReadInput(file_number, input_file)
      !------------------------------------------------------------------------------
      ! Initalize relevant matrices throughout the code.
      call inilag()
      !------------------------------------------------------------------------------
      ! Read all information from a wf file
      call ReadWavefunction()
      !------------------------------------------------------------------------------
      ! Print all relevant input gleaned from STDIN and the wf file.
      call PrintInput(file_number, input_file)
      !------------------------------------------------------------------------------
      ! Go out and try to reach convergence, only to fail time and time again....
      call ReachForWaterAndFood(iteration, iomsg)
      !------------------------------------------------------------------------------
      ! Perform analysis on the final many-body state
      ! (i) calculate and print the collective moment of inertias
      if (N_inertia .gt. 0) then
         call calc_collective_inertia
         if (MPI_RANK .eq. 0) then
            call print_collective_inertia
         end if
      end if
      !---------------------------------------------------------------------------
      ! Write other (optional) output files
      call write_advanced_output(iteration - 1, iomsg)
      !---------------------------------------------------------------------------
      ! Write output to the outputfile, i.e. the full wavefunction file
      call writewavefunction(12, outputfilename)
      !------------------------------------------------------------------------------
      ! Clean up after running, just in case we need to run again.
      call Cleanupthemess()
      !------------------------------------------------------------------------------
      ! Print all timing info
#if(USE_MPI > 0)
      ! guarantee that timing info is only at the very end
      call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
#endif
      call stop_timer(T_tantalus)
      call print_all_timers()
      !------------------------------------------------------------------------------
      ! end the processes across MPI ranks
#if(USE_MPI > 0)
      call mpi_finalize(mpi_err)
#endif

      ! end of one mean-field calculation..;
   end subroutine Run_Tantalus

   subroutine ReachForWaterAndFood(iter, iomsg)
      !---------------------------------------------------------------------------
      ! Evolve the single-particle wavefunctions and densities.
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
      !   |  6. Construct the potentials
      !   |     (including the contributions to F_I_I by Coulomb interaction
      !   |      and any multipole constraints)
      !   |  7. Print iteration info
      !   |_____________________________
      !
      ! TODO: correct this documentation
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Input:
      !   None.
      ! Output:
      !   iter  : number of iterations executed by this routine.
      !   iomsg : message about convergence that can be included in output files.
      !---------------------------------------------------------------------------
      use compilation
      use derivatives
      use wavefunctions
      use constants
      use densities
      use functional
      use evolution
      use IO
      use IO_wf, only: readHFBinfofile, write_tantalus_wf
#if(USE_HDF5>0)
      use IO_wf, only: write_tantalus_hdf5
#endif
      use moments
      use coulombmod
      use pairing
      use printing
      use momentsofinertia
      use cranking
      use convergence
      use scfiteration
      use timing

      implicit none

9     format(' Iter =', i5, '; writing checkpoint to file ', a20, '.')

      integer, intent(out)           :: iter
      character(len=99), intent(out) :: iomsg

      integer :: iprint, scheme, ifail
      logical :: ConvergenceAchieved, calc_expensive, print_all_spwf_properties
      logical :: potentials_frozen = .true.
      ! Logical to see if any moments with feasible set projection are necessary
      logical :: projectpresent = .false.

#if(USE_MPI > 0)
      integer :: mpi_err
#endif
#if(DEBUG_LEVEL == 1)
      character(len=40) :: denfile_iter, potfile_iter
#endif

      ifail = 0
      ConvergenceAchieved = .false.

      ! Provide memory for the derivatives of the spwfs
      call allocate_memory_derivatives(PairingType)

      !---------------------------------------------------------------------------
      ! Initial calculations
      !---------------------------------------------------------------------------
      if ((Bogofromfile .and. readHFBinfofile) .and. pairingscheme .eq. 1) then
         ! If using a gradient strategy and we want to continue from file.
         ! Only allowed of course if we have actually read a Bogoliubov transfo.
         scheme = -1
      else
         scheme = 0
      end if
      call SolvePairing(scheme, ifail)
      if (ifail .ne. 0) then
         ! Solve the pairing, with the current values of <h> and the pairing gaps.
         ! Note that this is ALWAYS a direct solve, i.e. we diagonalise the HFB
         ! Hamiltonian with a LAPACK call. We do this if the code did not receive
         ! explicit instructions to start from the Bogoliubov transformation on
         ! file.
         print *, 'WARNING! Pairing solver failed.'
      end if
      if (bogofromfile .and. guessgaps .and. pairingscheme .eq. 1) then
         ! We perform a few extra calls to solvepairing to take a few gradient
         ! steps, with finite values for Delta.
         call SolvePairing(pairingscheme, ifail)
      end if

      ! Derive all the single-particle wavefunctions in the HFPsi array
      if (store_derivatives) call deriveHF()
      ! Calculate the initial density vector
      call construct_canonical_basis(rho_pairing, kappa_pairing, rho_can, kappa_can)
      Density = densit(rho_can, kappa_pairing)

      call CalculateMoments(Density, .true.)
      !=> vital to be called here,
      !    (a) before the calculation of the potentials
      !    (b) after construction of the charge density
      ! as
      !  (a) the multipole cutoff is allocated in this
      !      process, and is needed for the calculation
      !      of the cranking potentials
      !  (b) the calculations of the charge rms radius
      !      requires the charge density to be
      !      constructed
      ! Adopt the relevant quantities to the centre-of-mass of the nucleus
      call adapt_com(Density)        !
      call CalculateMoments(Density, .false.) ! Recalculate because the COM might have changed.

      ! Update all spwf properties
#if(PASTA == 0)
      ! the memory and CPU time requirements of these routine scale very badly...
      call update_spwf_properties_HF() !
      if (PairingType .eq. 2) call update_spwf_properties_CAN()
      print_adv_spwf_properties = .true.
#else
      print_adv_spwf_properties = .false.
#endif

      call updateAM(Density, .true.) ! TODO: adapt the calculation of angular momentum
      !       to only ever use densities; this will avoid
      !       having to recalculate angular momentum matrix
      !       elements at every iteration

      ! Only calculate the fields that have not been read from either a
      ! wavefunction file or a potential file.
      if(allocated(potentials_read%F_I_I) .and. potentials_from_file) then
        potentials = calcPotentials(Density, potentials_read)
      else
        potentials = calcPotentials(Density)
      endif

      call setBelyaevProcedure()
      !---------------------------------------------------------------------------
      ! Calculate the energy WITH all the expensive parts included.
      call CalcEnergy(Density, Potentials, .true.)
      call calc_avg_gap()
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Initial printout
      call full_printout(0, .false., print_adv_spwf_properties)
      !---------------------------------------------------------------------------
      ! Start of the iterations
      !---------------------------------------------------------------------------
      do iter = 1, maxiter
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! First, do some bookkeeping
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! (1) update the arrays containing stuff at the last iteration
         call update_E_history()
         !     Save Fermi energy
         FermiHistory = FermiEnergy
         ! (2) check if some constraints should not be turned off
         call TurnOffConstraints(iter)
         ! (3) and decide whether we are constraining stuff or not
         projectpresent = checkconstraints() .or. check_cranking()

         !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Then, we update the reduced subspace spanned by our spwfs
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! TODO: include feasibleproject in the evolve_subspace code
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         if (projectpresent) call feasibleproject(Density)
         call Evolve_subspace(potentials, iter)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Calculate the single-particle hamiltonian ...
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         if (store_derivatives) call deriveHF() ! and update derivatives
         !  ... optionally perform a subspace rotation...
         if (subspace_rotation) then
            ! Full recalculation of h with the spwfs AFTER evolution
            sphamil = Calc_Sphamil(potentials, .true.)
            ! Mystery: performing this recalculation when subspace_rotation
            !          is NOT active will spoil convergence in many cases..
            call apply_subspace_rotation(sphamil, HFTransfo, spenergies)
            if (store_derivatives) call deriveHF() ! and update derivatives
         end if
         ! ..... and then calculate the pairing gaps
         call CalcGaps(FermiEnergy, PairStabFactor, Potentials)
         ! ... and use these matrices to build a new many-body state!
         call SolvePairing(pairingscheme, ifail)
         call construct_canonical_basis(rho_pairing, kappa_pairing,&
         &                              rho_can, kappa_can)
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! From the many-body state, we start calculating observables
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         Density = densit(rho_can, kappa_pairing)
         if (follow_com) call adapt_com(Density)
         ! Calculate the value of all multipole moments
         call CalculateMoments(Density, .true.)
         ! ...and readjust any constraints on them
         call ReadjustAllMoments(1) ! TODO: remove the input dependence here...
         call ReadjustAllMoments(2)
         ! Update value of the average angular momentum
         if (check_cranking() .and. .not. crank_smooth) then
            call update_spwf_properties_HF()
            if (PairingType .eq. 2) call update_spwf_properties_CAN()
         end if
         call updateAM(Density, .true.)
         ! .... and readjust any constraints on it
         call ReadjustCranking

         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Do a double take when constraints are present: use the updated
         ! Lagrange multipliers to correct our many-body state
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         !if(projectpresent) then
         !    ! Update the single-particle hamiltonian
         !    call update_sphamil_constraints(sphamil)
         !    if(subspace_rotation) then
         !        call apply_subspace_rotation(sphamil, HFTransfo, spenergies)
         !        if(store_derivatives) call deriveHF() ! and update derivatives
         !    endif
         !    ! .... and recalculate the gaps .....
         !    call CalcGaps(FermiEnergy, PairStabFactor, Potentials)
         !    ! ..... reconstruct a many-body state .....
         !    call construct_canonical_basis(rho_pairing,kappa_pairing,rho_can,kappa_can)
         !    Density = densit(rho_can, kappa_pairing)
         !    ! ..... reconstruct all densities ....
         !    if(follow_com) call adapt_com(Density)
         !    ! .... and recalculate constrained quantities
         !    call CalculateMoments(Density, .true.)
         !
         !    if(check_cranking() .and. .not. crank_smooth) then
         !      call update_spwf_properties_HF()
         !      if(PairingType.eq.2) call update_spwf_properties_CAN()
         !      call updateAM(Density, .true.)
         !    endif
         !endif

         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Construct new potentials ...
         ! ...but only if Iter > FreezeIter AND d2H < d2H_freeze
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         if ((iter .gt. freezeiter) .and. (d2H .lt. d2H_freeze)) then
            potentials_frozen = .false.
         end if

         if (.not. potentials_frozen) then
            ! calculate new values for the potentials from the densities
            potentials_out = calcPotentials(Density)

            if (scfscheme .eq. 0) then
               potentials_out = precondition_potentials(potentials, potentials_out)
            end if
            ! Save the information to memory, throwing out older information.
            ! This information is not used in any further part of the evolution
            ! if mixingscheme = 0.
            call save_potential_history(potentials, potentials_out)

            select case (mixingscheme)
            case (0)
               ! No mixing
               potentials = potentials_out
            case (1)
               ! Mixing with Anderson acceleration
               potentials_out = AndersonMixPotentials(Potential_iterates, &
               &                                      Potential_updates,  &
               &                                      mixstepsize, iter)
               potentials = potentials_out
            end select

         elseif (iter .eq. freezeiter) then
            ! Recalculate the Coulomb potential at the last iteration for
            ! comparison purposes with other codes.
            call solve_coulomb(Density, Potentials, sx_rho, sy_rho, sz_rho)
         end if
         !-----------------------------------------------------------------------
         ! Above: actual evolution of physical quantities
         ! Below: administration/bookkeeping
         !-----------------------------------------------------------------------
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Recalculate the energy with one of two options:
         ! - cheap calculation that omits the recalculation of some parts of the
         !   energy that are computationally intensive
         ! - expensive, complete calculation
         if ((mod(iter, PrintIter) .eq. 0) .or. (iter .eq. maxiter)) then
            iprint = 1; calc_expensive = .true.
         else
            iprint = 0
            calc_expensive = .false.
            if (check_cranking()) calc_expensive = .true. ! need the details of the spwf
            ! to calculate the cranking quantities
         end if

         call CalcEnergy(Density, Potentials, calc_expensive)
         ! Calculate the average pairing gap
         call calc_avg_gap()
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Check for convergence or a failed calculation
         ! TODO: what is this?
         if (ifail .ne. 0) then
            iomsg = 'FERMI'
            ConvergenceAchieved = .false.
            exit
         else
            call Converged(ConvergenceAchieved, iter)
         end if

         if (convergenceAchieved) then
            iprint = 1
            ! Recalculate the energy with all parts included at the end, don't
            ! skimp on the expensive parts
            call CalcEnergy(Density, Potentials, .true.)
         end if
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         !  update all spwf properties first to ensure correct printout of spwfs
#if(PASTA == 0)
         ! the memory and CPU time requirements of these routine scale very badly...
         print_all_spwf_properties = print_adv_spwf_properties .or. &
         &                           (iter .eq. maxiter) .or. &
         &                           convergenceachieved
         if (print_all_spwf_properties) then
            call update_spwf_properties_HF()
            if (PairingType .eq. 2) call update_spwf_properties_CAN()
         end if
#else
         print_all_spwf_properties = .false.
#endif
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Decide whether to do a full printout
         ! .... but do a summary printout anyway to enable for "complete" output
         !      when grepping on quantities included in the summary
         if (MPI_RANK .eq. 0) call printsummary(iter, potentials_frozen)
         if (iprint .eq. 1) then
            call full_printout(iter, convergenceachieved, print_all_spwf_properties)
         end if
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Exit the loop if convergence is achieved.
         if (ConvergenceAchieved) then
            call print_convergence_message(iter)
            iomsg = 'CONVERGED'
            exit
         end if
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Write a wavefunction file at each multiple of checkpointiter
         if (checkpointiter .ne. 0) then
         if (mod(iter, checkpointiter) .eq. 0) then
#if(DEBUG_LEVEL==1)
            ! Output densities and potentials to specific files at every checkpoint
            ! TODO: refactor this into a subroutine in the IO.f90 module!
            write (denfile_iter, '("iter=",i5.5,".den")') iter
            write (potfile_iter, '("iter=",i5.5,".pot")') iter
            if (MPI_RANK .eq. 0) then
               call write_densities(Density, denfile_iter)
               call write_potentialfile(potentials, potfile_iter)
            end if
#endif
            if (MPI_RANK .eq. 0) print 9, iter, outputfilename
            iomsg = 'CHECKPOINT'
            if (trim(to_upper(OutputFileName(len_trim(OutputFileName) - 3:))) .eq. 'HDF5') then
#if(USE_HDF5>0)
               call write_tantalus_hdf5(outputfilename) !new hdf5 format
#else
               call stp('HDF5 support was not enabled at compilation.')
#endif
            else
               call write_tantalus_wf(12, outputfilename) ! old style in .wf file
            end if
         end if
         end if
      end do
   end subroutine ReachForWaterAndFood
#endif

   subroutine printsummary(iter, potentials_frozen)
      !---------------------------------------------------------------------------
      ! Short printout after an iteration with sufficient information to follow
      ! somewhat the convergence of the calculation.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Input:
      !    iter              : iteration count
      !    potentials_frozen : whether or not the potentials were updated
      !---------------------------------------------------------------------------
      use functional
      use evolution
      use moments
      use pairing

      implicit none

      integer, intent(in)   :: iter
      logical, intent(in)   :: potentials_frozen
      type(Moment), pointer :: current, part
      real(KIND=dp)         :: dF(2), DN(2), dQ, dL, dev, val, devJ
      character(len=1)      :: t, spec

1     format(86('-'))
2     format(' Iteration = ', i4)
21    format(' Potentials frozen.')
3     format(' dt    = ', f10.4, 4x, '  mu   = ', f10.4, ' gradn = ', es12.3, ' D2H  = ', es12.3)
31    format(' dtg   = ', f10.4, 4x, '  mug  = ', f10.4, ' gradn = ', es12.3)
4     format(' E     = ', f20.10, 2x, '  DE   = ', e12.5)
41    format(' R     = ', f20.10, 2x, '  DR   = ', e12.5)
42    format(' R-E   = ', f20.10, 2x, 'D(R-E) = ', e12.5)
#if(PASTA == 1)
43    format(' Epasta= ', f20.10, 2x, 'DEpasta= ', e12.5)
#endif
5     format(' ', a1, 'Q', 2i1, a1, ' = ', f12.4, 3x, 'dQ = ', es8.1, 2x,         &
                          &          'L = ', f12.4, 2x, ' dL = ', es8.1, 2x, 'dev = ', es8.1)

6     format(' dmun  = ', es8.1, 4x, '  dmup= ', es8.1)
7     format(' dN    = ', es8.1, 4x, '  dZ  = ', es8.1)
8     format(' Jz    = ', f12.4, 2x, ' dJZ= ', es8.1,  &
                          &         ' Om = ', f12.4, 2x, ' dO = ', e8.1, 2x, 'dev = ', es8.1)

      part => FindMoment(0, 0, .false.)

      if (iter .eq. 1) print 1
      print 2, iter
      if (potentials_frozen) print 21
      print 3, dt, momentum, gradientnorm, d2h
      if (pairingscheme .eq. 1) then
         print 31, gradient_stepsize, gradient_mu, sqrt(sum(HFBGradnorm**2))
      end if
      print 4, totalE, (totalE - Ehistory(1))/abs(totalE)
      print 41, Routhian, (Routhian - Rhistory(1))/abs(Routhian)
      print 42, Routhian - totalE, &
      &  ((Routhian - Rhistory(1)) - (totalE - Ehistory(1)))/abs(totalE)
#if(PASTA == 1)
      print 43, calculate_epasta(totalE), (calculate_epasta(totalE) - calculate_epasta(Ehistory(1)))/abs(calculate_epasta(totalE))
#endif
      if (fixfermi) then
         dN = part%value - part%history
         print 7, dN
      else
         dN(1) = sum(rho_can(1:nwn)) - neutrons
         dN(2) = sum(rho_can(nwn + 1:nwt)) - protons
         if (any(dN(:)*dN(:) .gt. 1.d-14)) then
            print 7, dN
         end if
         dF = FermiEnergy - FermiHistory
         print 6, dF
      end if

      Current => Root
      do while (associated(Current%next))
         Current => Current%next

         if ((Current%constrainttype .ne. 0) .or. (Current%l .eq. 2)) then
            dL = Current%multiplier - Current%mult_hist
            if (Current%Impart) then
               t = 'I'
            else
               t = 'R'
            end if

            if (Current%constrainttype .ne. 0) then
               dev = Current%deviation
            else
               dev = 0
            end if

            if (current%l .eq. 2 .and. current%isoswitch .ne. 0) then
               dQ = sum(Current%value) - sum(Current%history)
               print 5, t, Current%l, Current%m, 't', sum(Current%value), dQ, &
               &              0.0d0, 0.0d0, 0.0d0
            end if

            select case (Current%isoswitch)
            case (0)
               dQ = sum(Current%value) - sum(Current%history)
               spec = 't'
               val = sum(Current%value)
            case (1, 2)
               dQ = Current%value(Current%isoswitch)&
                & - Current%history(Current%isoswitch)
               val = Current%value(Current%isoswitch)
               select case (Current%isoswitch)
               case (1)
                  spec = 'n'
               case (2)
                  spec = 'p'
               end select
            end select
            print 5, t, Current%l, Current%m, spec, val, &
              &         dQ, Current%multiplier, dL, dev
         end if
      end do

      if (.not. crank_smooth) then
         if (cranktype(3) .eq. 1) then
            devJ = TotalAngMom(3) - CrankValues(3)
         else
            devJ = 0.0d0
         end if
         print 8, totalangmom(3), totalangmom(3) - angmomold(3), &
         &        omega(3), omega(3) - omega_prev(3), devJ
      else
         if (cranktype(3) .eq. 1) then
            devJ = TotalAngMom_dens(3) - CrankValues(3)
         else
            devJ = 0.0d0
         end if
         print 8, totalangmom_dens(3), totalangmom_dens(3) - angmomold_dens(3), &
         &        omega(3), omega(3) - omega_prev(3), devJ
      end if

      print 1

   end subroutine printsummary

   subroutine update_spwf_properties_HF()
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Calculate all relevant properties of the spwf in the HFbasis.
      ! This is just a wrapper function that calls update_spwf_properties with
      ! the correct input depending on the type of calculation we are performing.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use evolution, only: diagsphamil
      use wavefunctions, only: update_spwf_properties, HFTransfo
      use wavefunctions, only: HF_J, HF_JTR, HF_JTI, HF_JJ, HF_J2
      use wavefunctions, only: HF_spin, HF_STR, HF_STI
      use wavefunctions, only: spwf_r2_HF, P_HF

      if (diagsphamil) then
         call update_spwf_properties('HF', HFTRANSFO, .false.,           &
                                 &  HF_J, HF_JTR, HF_JTI, HF_J2, HF_JJ, & ! J-like stuff
                                 &  HF_spin, HF_STR, HF_STI,            & ! spin-stuff
                                 &  spwf_r2_hf,                         & ! radii
                                 &  P_hf)                                 ! symmetry-stuff
      else
         call update_spwf_properties('HF', HFTRANSFO, .true.,    &
                                 &  HF_J, HF_JTR, HF_JTI, HF_J2, HF_JJ, & ! J-like stuff
                                 &  HF_spin, HF_STR, HF_STI,            & ! spin-stuff
                                 &  spwf_r2_hf,                         & ! radii
                                 &  P_hf)                                 ! symmetry-stuff
      end if

   end subroutine update_spwf_properties_HF

   subroutine update_spwf_properties_CAN()
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Calculate all relevant properties of the spwf in the HFbasis.
      ! This is just a wrapper function that calls update_spwf_properties with
      ! the correct input depending on the type of calculation we are performing.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use pairing, only: cantransfo
      use wavefunctions, only: update_spwf_properties, nwn
      use wavefunctions, only: CAN_J, CAN_JTR, CAN_JTI, CAN_JJ, CAN_J2
      use wavefunctions, only: CAN_spin, CAN_STR, can_STI
      use wavefunctions, only: spwf_r2_can, P_can

      call update_spwf_properties('CAN', CANTRANSFO, .false.,                &
                                &  CAN_J, CAN_JTR, CAN_JTI, CAN_J2, CAN_JJ, & ! J-like stuff
                                &  CAN_spin, CAN_STR, CAN_STI,              & ! spin-stuff
                                &  spwf_r2_can,                             & ! radii
                                &  P_can)                                     ! symmetry-stuff
   end subroutine update_spwf_properties_CAN

   subroutine full_printout(iter, converged, print_all_spwf_properties)
      !-----------------------------------------------------------------------------
      ! Perform a complete print out of the entire state of the code
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Input:
      !    iter                      : iteration count
      !    converged                 : logical, if .true. the calculation converged
      !    print_all_spwf_properties : logical, if .true. print ALL details on the
      !                                spwf
      ! Output:
      !    NONE
      !-----------------------------------------------------------------------------
      use pairing, only: printpairing
      use moments, only: printallmoments
      use densities, only: density
#if($FAM == 0)
      use densities, only: print_boxsize_check
#endif
      use momentsofinertia, only: printMomentsOfInertia
      use printing, only: printqps, print_adv_spwf_properties, printspwfs
      use cranking, only: printcranking
      use functional, only: printenergy, PairStabFactor

1     format(86('-'))
2     format(30x, 'Iteration = ', i5,/)
3     format(24x, 'FINAL Iteration = ', i5,/)

      integer, intent(in) :: iter
      logical, intent(in) :: print_all_spwf_properties, converged

      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Drive all routines to print their output to STDOUT
      if (MPI_RANK .eq. 0) then
         print 1
         if ((iter .eq. maxiter) .or. converged) then
            ! Add a clear indication this is the FINAL iteration
            print 3, iter
         else
            print 2, iter
         end if
#if(PASTA == 0 && DEBUG_LEVEL== 0)
         ! Pasta calculations typically involve TONS of spwfs
         ! .... but we might be interested in their properties when debugging!
         call printspwfs(print_all_spwf_properties, print_last=0)
         call printqps(print_all_spwf_properties)
#else
         call printspwfs(print_all_spwf_properties, print_last=100)
#endif
         call printallmoments
#if(PASTA == 0)
#if($FAM == 0)
         call print_boxsize_check(Density)
#endif
         call printmomentsofinertia
         call printcranking(Density)
#endif
         call printpairing(pairstabfactor)
         call printenergy()
      end if

   end subroutine full_printout

   subroutine print_convergence_message(iter)
      !-----------------------------------------------------------------------------
      ! Print a convergence message if the calculation converged.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Input:
      !    iter : iteration count
      ! Output:
      !    NONE
      !-----------------------------------------------------------------------------

      use geninfo

      integer, intent(in) :: iter

1     format('----------------------------------')
2     format('| Convergence criteria satisfied.|')
3     format('| Needed ', i4, ' iterations.', 8x, '|')
4     format('| dE    < ', es10.3, 12x, ' | ')
5     format('| dQ2   < ', es10.3, 12x, ' | ')
6     format('| d2H   < ', es10.3, 12x, ' | ')
61    format('| |spg| < ', es10.3, 12x, ' | ')
7     format('| dmu   < ', es10.3, 12x, ' | ')
71    format('| dJz   < ', es10.3, 12x, ' | ')
8     format('| Ending the iterative proces.   |')

      if (MPI_RANK .eq. 0) then
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
      end if

   end subroutine print_convergence_message

   subroutine initialize_all_timers(fam)
      !----------------------------------------------------------------------------
      ! Initialize all the timers that have been defined.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Input:
      !       fam (logical), optional :  add fam timers
      ! Output:
      !       NONE
      !----------------------------------------------------------------------------
      use timing

      logical, intent(in), optional :: fam

      call add_timer('Tantalus', T_tantalus)
      call add_timer('Wavefunction initialisation', T_wfini)
      call add_timer('Wavefunction output', T_wfoutput)
      call add_timer('Wavefunction reading', T_wfinput)
      call add_timer('HF-basis Derivatives', T_derivatives)
      call add_timer('Canonical basis Derivatives', T_derivatives_can)
      call add_timer('Spwf evolution', T_evolution)
      call add_timer('Orthonormalization', T_ortho)
      call add_timer('Construction of norm matrix', T_norm_ortho)
      call add_timer('Diagonalisation norm matrix', T_diag_ortho)
      call add_timer('Density calculations', T_densities)
      call add_timer('Density: pp', T_den_pp)
      call add_timer('Density: ph', T_den_ph)
      call add_timer('Density: derivatives', T_den_der)
      call add_timer('Potential calculations', T_potentials)
      call add_timer('Potential preconditioning', T_pot_precon)
      call add_timer('Energy calculations', T_energy)
      call add_timer('Pairing solver ', T_pairing)
      call add_timer('Sp. Hamiltonian ', T_sphamil)
      call add_timer('Coulomb solver', T_coulomb)
      call add_timer('Can. basis construction', T_den_can)
      call add_timer('Moments of inertia ', T_MOI)
      call add_timer('Centre-of-mass correction ', T_COM)
      call add_timer('COM one-body ', T_COM1)
      call add_timer('COM two-body ', T_COM2)
      call add_timer('Matrix elements summation', T_COM2_summation)
      call add_timer('Matrix elements of \nabla ', T_NablaMElements)
      call add_timer('Pairing gaps ', T_gaps)
      call add_timer('Multipole moments ', T_moments)
      call add_timer('Multipole moments cutoff', T_moment_cutoff)
      call add_timer('Feas. Proj. step ', T_feasible)
      call add_timer('Spwf angular momentum ', T_spwfangmom)
      call add_timer('Charge density folding', T_chargedensity)
      call add_timer('Collective MOIs', T_collective_moi)
      call add_timer('Microscopic pairing', T_microscopic_pairing)
      call add_timer('Subspace rotation          ', T_Hortho)
      call add_timer('Construction HF transfo', T_HFDiag)
      call add_timer('Basis transformation', T_Basistransfo)
      call add_timer('Subspace rotation', T_subspace_rotation)
      call add_timer('Spwf transformation', T_subrot_transfo)
      call add_timer('Matrix diagonalisation', T_subrot_diag)
      call add_timer('Calculation of h in subspace', T_calc_sph)
      call add_timer('Matrix elements of h', T_calc_sph_me)
      call add_timer('Update of h in subspace', T_update_sph)
#if( USE_MPI > 0)
      call add_timer('Layout transfer: 1D -> 2D', T_transfer_psi_1to2)
      call add_timer('Layout transfer: 2D -> 1D', T_transfer_psi_2to1)
      call add_timer('MPI_ALLREDUCE calls     ', T_allreduce)
#endif
      if (present(fam)) then
         call add_timer('FAM', T_fam)
         call add_timer('Perturbation densities', T_den_perturbed)
         call add_timer('Sym. pert. densities', T_den_perturbed_sym)
         call add_timer('Anti pert. densities', T_den_perturbed_asym)
         call add_timer('Matrix elements \delta h', T_spme_perturbed)
         call add_timer('Sym. \delta h', T_spme_perturbed_sym)
         call add_timer('Anti \delta h', T_spme_perturbed_asym)
         call add_timer('PP  per. densities', T_den_perturbed_pp)
         call add_timer('Matrix elements \Delta', T_spme_perturbed_pp)
      end if

   end subroutine initialize_all_timers

   subroutine cleanupthemess()
      !-----------------------------------------------------------------------------
      ! Driver routine calling all routines to cleanup memory in different modules.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Input:
      !       NONE
      ! Output:
      !       NONE
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
      call clean_moments
      call clean_coulomb
      call clean_evolution

   end subroutine cleanupthemess
end module Tantalus
