!===============================================================================
!     __  __  ___   ____ ____
!    |  \/  |/ _ \ / ___/ ___|__ _
!    | |\/| | | | | |  | |   / _` |
!    | |  | | |_| | |__| |__| (_| |
!    |_|  |_|\___/ \____\____\__,_|
!
!    Copyright (C) 2026 W. Ryssens and M. Bender
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU Affero General Public License as published
!    by the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU Affero General Public License for more details.
!
!    You should have received a copy of the GNU Affero General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!
!===============================================================================
module MOCCa

   implicit none (type, external)
   public

contains

#if( $FAM == 0)
   subroutine Run_MOCCa(file_number, input_file)
      !==============================================================================
      ! TODO: describe input/output of this routine
      !
      ! While I (W.R.) like to think about MOCCa as a standalone code, this
      ! 'main' routine is now written as a subroutine (with inputs!) to accomodate
      ! meta-codes that want to run MOCCa multiple times.
      !==============================================================================
      use compilation,   only: dp
      use geninfo,       only: NPROCS, MPI_RANK
      use IO,            only: ReadInput, PrintInput, write_advanced_output
      use IO,            only: readwavefunction, outputfilename, writewavefunction
      use fission_MOI,   only: calc_collective_inertia, print_collective_inertia
      use fission_MOI,   only: N_inertia
      use version,       only: print_header
      use derivatives,   only: inilag
      use timing,        only: start_timer, stop_timer, print_all_timers, T_MOCCa,&
           &                   initialize_all_timers
#if(USE_MPI > 0)
      use mpi_f08
      ! This include statement is not particularly elegant, but appending it with an
      ! 'only'-list seems to generate behaviour that is not consistent across compilers.
#endif


      !------------------------------------------------------------------------------
      ! Convergence signals
      ! Iteration counter
      integer           :: iteration
      ! Message for the output of the code, useful for the Brussels group.
      character(len=99) :: iomsg
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
      call start_timer(T_MOCCa)

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
      call run_mean_field(iteration, iomsg)
      !------------------------------------------------------------------------------
      ! Perform analysis on the final many-body state
      ! (i) calculate and print the collective moment of inertias
      if (N_inertia > 0) then
         call calc_collective_inertia
         if (MPI_RANK == 0) then
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
      ! Print all timing info
#if(USE_MPI > 0)
      ! guarantee that timing info is only at the very end
      call MPI_Barrier(MPI_COMM_WORLD, mpi_err)
#endif
      call stop_timer(T_MOCCa)
      call print_all_timers()
      !------------------------------------------------------------------------------
      ! end the processes across MPI ranks
#if(USE_MPI > 0)
      call mpi_finalize(mpi_err)
#endif

      ! end of one mean-field calculation..;
   end subroutine Run_MOCCa

   subroutine run_mean_field(iter, iomsg)
     !---------------------------------------------------------------------------
     !  We start the mean-field process.
     !
     !  Output:
     !   iter : integer, the final iteration performed in the SCF process
     !   iomsg: character, message about convergence
     !----------------------------------------------------------------------------
     use scfiteration, only : scfscheme

     integer, intent(out)           :: iter
     character(len=99), intent(out) :: iomsg

     ! 1. set-up phase, construction of the starting point
     call set_up_mean_field()
     ! 2. Iteration phase, start of iterations
     if(scfscheme /= 2) then
        call iterate_mean_field(iter, iomsg)
     else
        call iterate_spectrum(iter, iomsg)
     endif

   end subroutine run_mean_field

   subroutine set_up_mean_field()
      !
      ! Use statements with ONLY clauses
      ! TODO document
      !
      use pairing,       only: BogoFromFile, PairingType, pairingscheme, &
           &                   guessgaps, SolvePairing, rho_pairing,     &
           &                   kappa_pairing, rho_can, kappa_can,        &
           &                   calc_avg_gap
      use geninfo,       only: store_derivatives
      use densities,     only: Density, densit, construct_canonical_basis
      use moments,       only: CalculateMoments, adapt_COM
      use functional,    only: potentials_read,calcPotentials,           &
           &                   CalcEnergy, potentials
      use printing,      only: print_adv_spwf_properties
      use cranking,      only: updateAM
      use IO,            only: potentials_from_file
      use IO_wf,         only: readHFBinfofile
      use momentsofinertia, only: setBelyaevProcedure
      use wavefunctions, only: allocate_memory_derivatives, deriveHF


      integer :: iprint, scheme, ifail

#if(USE_MPI > 0)
      integer :: mpi_err
#endif
#if(DEBUG_LEVEL == 1)
      character(len=40) :: denfile_iter, potfile_iter
#endif

      ifail = 0

      ! Provide memory for the derivatives of the spwfs
      call allocate_memory_derivatives(PairingType)

      !---------------------------------------------------------------------------
      ! Initial calculations
      !---------------------------------------------------------------------------
      if ((Bogofromfile .and. readHFBinfofile) .and. pairingscheme == 1) then
         ! If using a gradient strategy and we want to continue from file.
         ! Only allowed of course if we have actually read a Bogoliubov transfo.
         scheme = -1
      else
         scheme = 0
      end if
      call SolvePairing(scheme, ifail)
      if (ifail /= 0) then
         ! Solve the pairing, with the current values of <h> and the pairing gaps.
         ! Note that this is ALWAYS a direct solve, i.e. we diagonalise the HFB
         ! Hamiltonian with a LAPACK call. We do this if the code did not receive
         ! explicit instructions to start from the Bogoliubov transformation on
         ! file.
         print *, 'WARNING! Pairing solver failed.'
      end if
      if (bogofromfile .and. guessgaps .and. pairingscheme == 1) then
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
      if (PairingType == 2) call update_spwf_properties_CAN()
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

   end subroutine set_up_mean_field

   subroutine iterate_mean_field(iter, iomsg)
     !-----------------------------------------------------------------------------
     ! This routine solves the self-consistent mean-field problem, i.e. it performs
     ! the whole internested updating of both
     !
     !   1. single-particle wavefunctions
     !   2. mean-field potentials
     !
     ! as described in
     !   W. Ryssens, M. Bender, M. and P.-H. Heenen,
     !   Eur. Phys. J. A, 55, 93. https://doi.org/10.1140/epja/i2019-12766-6
     !
     ! Output
     !   iter : integer, the final iteration performed in the SCF process
     !   iomsg: character, message about convergence
     !
     !---------------------------------------------------------------------------
      use compilation,   only: dp
      use geninfo,       only: MaxIter, d2H_freeze, FreezeIter, PrintIter, &
           &                   store_derivatives, to_upper, stp
      use wavefunctions, only: sphamil, HFTransfo, spenergies, deriveHF
      use densities,     only: Density, densit, construct_canonical_basis
      use densities,     only: sx_rho, sy_rho, sz_rho
      use functional,    only: calcPotentials, potentials, potentials_out,      &
           &                   precondition_potentials, save_potential_history, &
           &                   CalcEnergy, PairStabFactor, update_E_history
      use evolution,     only: Evolve_subspace, Calc_Sphamil, subspace_rotation,&
           &                   apply_subspace_rotation, d2H, gradientnorm, &
           &                   feasibleproject
      use IO,            only: OutputFileName, MPI_RANK, checkpointiter, N_inertia
      use IO_wf,         only: write_MOCCa_wf
#if(USE_HDF5>0)
      use IO_wf,         only: write_MOCCa_hdf5
#endif
      use moments,       only: CalculateMoments, ReadjustAllMoments, follow_com, &
           &                   checkconstraints, turnoffconstraints, adapt_com
      use coulombmod,    only: solve_coulomb
      use pairing,       only: SolvePairing, pairingscheme, calcgaps,           &
      &                        rho_pairing, kappa_pairing, rho_can, kappa_can,  &
      &                        PairingType, calc_avg_gap, FermiEnergy, FermiHistory
      use printing,      only: print_adv_spwf_properties
      use convergence,   only: Converged
      use scfiteration,  only: scfscheme, mixingscheme, mixstepsize, freezeiter, &
      &                        maxiter, AndersonMixPotentials,               &
      &                        Potential_iterates, Potential_updates

      use timing,        only: start_timer, stop_timer
      use cranking,      only: check_cranking, crank_smooth, cranktype, CrankValues, &
                               TotalAngMom, omega, angmomold, omega_prev, totalangmom_dens, &
                               angmomold_dens, updateAM, readjustCranking, check_cranking


     integer, intent(out)           :: iter
     character(len=99), intent(out) :: iomsg

     logical :: ConvergenceAchieved, calc_expensive, print_all_spwf_properties
     logical :: potentials_frozen
     ! Logical to see if any moments with feasible set projection are necessary
     logical :: projectpresent

     integer :: ifail, iprint

9    format(' Iter =', i5, '; writing checkpoint to file ', a20, '.')

     ConvergenceAchieved = .false.
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
            if (PairingType == 2) call update_spwf_properties_CAN()
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
         if ((iter > freezeiter) .and. (d2H < d2H_freeze)) then
            potentials_frozen = .false.
         end if

         if (.not. potentials_frozen) then
            ! calculate new values for the potentials from the densities
            potentials_out = calcPotentials(Density)

            if (scfscheme == 0) then
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

         elseif (iter == freezeiter) then
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
         if ((mod(iter, PrintIter) == 0) .or. (iter == maxiter)) then
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
         if (ifail /= 0) then
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
         print_all_spwf_properties = (print_adv_spwf_properties .or. &
         &                           (iter == maxiter)) .or. &
         &                           convergenceachieved
         if (print_all_spwf_properties) then
            call update_spwf_properties_HF()
            if (PairingType == 2) call update_spwf_properties_CAN()
         end if
#else
         print_all_spwf_properties = .false.
#endif
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Decide whether to do a full printout
         ! .... but do a summary printout anyway to enable for "complete" output
         !      when grepping on quantities included in the summary
         if (MPI_RANK == 0) call printsummary(iter, potentials_frozen)
         if (iprint == 1) then
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
         if (checkpointiter /= 0) then
            if (mod(iter, checkpointiter) == 0) then
#if(DEBUG_LEVEL==1)
               ! Output densities and potentials to specific files at every checkpoint
               ! TODO: refactor this into a subroutine in the IO.f90 module!
               write (denfile_iter, '("iter=",i5.5,".den")') iter
               write (potfile_iter, '("iter=",i5.5,".pot")') iter
               if (MPI_RANK == 0) then
                  call write_densities(Density, denfile_iter)
                  call write_potentialfile(potentials, potfile_iter)
               end if
#endif
               if (MPI_RANK == 0) print 9, iter, outputfilename
               iomsg = 'CHECKPOINT'
               if (trim(to_upper(OutputFileName(len_trim(OutputFileName) - 3:))) == 'HDF5') then
#if(USE_HDF5>0)
                  call write_MOCCa_hdf5(outputfilename) !new hdf5 format
#else
                  call stp('HDF5 support was not enabled at compilation.')
#endif
               else
                  call write_MOCCa_wf(12, outputfilename) ! old style in .wf file
               end if
            end if
         end if
      end do
   end subroutine iterate_mean_field

   subroutine iterate_spectrum(iter, iomsg)
     !-------------------------------------------------------------------------
     !
     ! Output:
     !   iter : integer, the final iteration count
     !   iomsg: string, message about (non)convergence
     !-------------------------------------------------------------------------

     use compilation,   only: dp
     use geninfo,       only: MaxIter, PrintIter, store_derivatives
     use wavefunctions, only: locked, deriveHF, sphamil, HFtransfo, spenergies
     use evolution,     only: Evolve_subspace, calc_sphamil, subspace_rotation
     use evolution,     only: apply_subspace_rotation
     use functional,    only: potentials
     use pairing,       only: pairingtype

     integer, intent(out)           :: iter
     character(len=99), intent(out) :: iomsg
     logical                        :: convergence_achieved
     logical                        :: print_all_spwf_properties

     do iter=1, maxiter
         ! Just keep doing heavy-ball steps
         call Evolve_subspace(potentials, iter)
         if (store_derivatives) call deriveHF() ! and update derivatives

         ! Full recalculation of h with the spwfs AFTER evolution
         sphamil = Calc_Sphamil(potentials, .true.)
         call apply_subspace_rotation(sphamil, HFTransfo, spenergies)
         if (store_derivatives) call deriveHF() ! and update derivatives

         call print_summary_spectrum(iter)
         ! Convergence is reached if ALL the spwfs are locked inside the
         !  evolution,i.e. if their dispersion is small enough
         convergence_achieved = ALL(locked)

         if (convergence_achieved .or. iter .eq. maxiter) then
            call update_spwf_properties_HF()
            if (PairingType == 2) call update_spwf_properties_CAN()
            call full_printout(iter, convergence_achieved, .true.)
         end if
         ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         ! Exit the loop if convergence is achieved.
         if (convergence_achieved) then
            call print_convergence_message(iter)
            iomsg = 'CONVERGED'
            exit
         end if
     enddo

   end subroutine iterate_spectrum

   subroutine printsummary(iter, potentials_frozen)
      !---------------------------------------------------------------------------
      ! Short printout after an iteration with sufficient information to follow
      ! somewhat the convergence of the calculation.
      ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      ! Input:
      !    iter              : iteration count
      !    potentials_frozen : whether or not the potentials were updated
      !---------------------------------------------------------------------------
      use compilation,   only: dp
      use geninfo,       only: neutrons, protons, MPI_RANK
      use wavefunctions, only: nwn, nwt
      use functional,    only: totalE, Ehistory, Routhian, Rhistory
#if(PASTA == 1)
      use functional,   only: calculate_epasta
#endif
      use evolution,    only: dt, momentum, gradientnorm, d2h
      use moments,      only: moment, findmoment, Root
      use pairing,      only: gradient_stepsize, gradient_mu, HFBgradnorm, rho_can, &
           &                  fixfermi, pairingscheme, FermiEnergy, FermiHistory
      use cranking,     only: crank_smooth, cranktype, TotalAngMom, CrankValues, &
           &                  angmomold, omega, omega_prev, TotalAngMom_dens,    &
           &                  angmomold_dens

      integer, intent(in)   :: iter
      logical, intent(in)   :: potentials_frozen
      type(Moment), pointer :: current, part
      real(KIND=dp)         :: dF(2), DN(2), dQ, dL, dev, val, devJ, dE_pasta
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

      if (iter == 1) print 1
      print 2, iter
      if (potentials_frozen) print 21
      print 3, dt, momentum, gradientnorm, d2h
      if (pairingscheme == 1) then
         print 31, gradient_stepsize, gradient_mu, sqrt(sum(HFBGradnorm**2))
      end if
      print 4, totalE, (totalE - Ehistory(1))/abs(totalE)
      print 41, Routhian, (Routhian - Rhistory(1))/abs(Routhian)
      print 42, Routhian - totalE, &
      &  ((Routhian - Rhistory(1)) - (totalE - Ehistory(1)))/abs(totalE)
#if(PASTA == 1)
      dE_pasta = (calculate_epasta(totalE) - calculate_epasta(Ehistory(1))) &
           &     /abs(calculate_epasta(totalE))
      print 43, calculate_epasta(totalE),dE_pasta
#endif
      if (fixfermi) then
         dN = part%value - part%history
         print 7, dN
      else
         dN(1) = sum(rho_can(1:nwn)) - neutrons
         dN(2) = sum(rho_can(nwn + 1:nwt)) - protons
         if (any(dN(:)*dN(:) > 1.d-14)) then
            print 7, dN
         end if
         dF = FermiEnergy - FermiHistory
         print 6, dF
      end if

      Current => Root
      do while (associated(Current%next))
         Current => Current%next

         if ((Current%constrainttype /= 0) .or. (Current%l == 2)) then
            dL = Current%multiplier - Current%mult_hist
            if (Current%Impart) then
               t = 'I'
            else
               t = 'R'
            end if

            if (Current%constrainttype /= 0) then
               dev = Current%deviation
            else
               dev = 0
            end if

            if (current%l == 2 .and. current%isoswitch /= 0) then
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
         if (cranktype(3) == 1) then
            devJ = TotalAngMom(3) - CrankValues(3)
         else
            devJ = 0.0d0
         end if
         print 8, totalangmom(3), totalangmom(3) - angmomold(3), &
         &        omega(3), omega(3) - omega_prev(3), devJ
      else
         if (cranktype(3) == 1) then
            devJ = TotalAngMom_dens(3) - CrankValues(3)
         else
            devJ = 0.0d0
         end if
         print 8, totalangmom_dens(3), totalangmom_dens(3) - angmomold_dens(3), &
         &        omega(3), omega(3) - omega_prev(3), devJ
      end if

      print 1

   end subroutine printsummary

   subroutine print_summary_spectrum(iter)
     !------------------------------------------------------------------
     ! Print a summary on the progress of the iterate_spectrum routine.
     !
     !------------------------------------------------------------------

     use evolution,     only : dt, momentum, gradientnorm, d2h
     use wavefunctions, only : dispersions

     integer, intent(in) :: iter

1    format(86('-'))
2    format(' Iteration = ', i4)
21   format(' Potentials frozen.')
3    format(' dt    = ', f10.4, 4x, '  mu   = ', f10.4, ' gradn = ', es12.3, ' D2H  = ', es12.3)
4    format(' Maximum dispersion =', es12.3)

     print 1
     print 2, iter
     print 3, dt, momentum, gradientnorm, d2h
     print 4,  maxval(dispersions)

   end subroutine print_summary_spectrum
#endif

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
      use geninfo, only: MPI_RANK, NPROCS, maxiter
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
      if (MPI_RANK == 0) then
         print 1
         if ((iter == maxiter) .or. converged) then
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
      use geninfo, only : energy_prec, moment_prec, disp_prec, gradient_prec, &
           &              fermi_prec, angmom_prec, MPI_RANK

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

      if (MPI_RANK == 0) then
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
end module MOCCa
