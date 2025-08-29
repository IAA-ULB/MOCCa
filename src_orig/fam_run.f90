program run_FAM

  use compilation
  use IO
  use Tantalus, only : print_header, initialize_all_timers, full_printout
  use Tantalus, only : update_spwf_properties_HF, update_spwf_properties_CAN
  use fam
  use fam_testing, only : run_FAM_tests, test_gmres, test_gmres_affine, test_linearity_T
  use gmres 
  use timing

  1 format(86('-'))
  11 format(10(' '), 20('='), ' omega =', f5.2, ' MeV ', 20('='), 10(' '))
  2 format('FAM iteration = ', i5) 
  3 format(' S_',i1,i1,' (', f5.2, ') = ', es10.3)

  implicit none
  integer :: iteration, nbprod
  logical :: is_converged, is_divergent
  real(kind=dp) :: omega_curr
  integer :: omega_num, omega_index
  real(kind=dp), allocatable :: omega_arr(:), S_arr(:), S_free_arr(:)
  integer, allocatable :: iter_arr(:)
  character(len=100) :: famfilename
  ! integer :: i, B, si,N

  complex(KIND=dp), allocatable :: dH_flat(:), dH_flat_next(:)
  real(KIND=dp) :: res

  ! integer :: ifail ! Future dev: required for HFB

  ! Print a nice header with all kinds of relevant info
  call print_header(.true.)

  !------------------------------------------------------------------------------
  ! starting all timers
  ! 
  ! -> This is necessary since subroutines below make use of the timers
  ! 
  call initialize_all_timers(.true.) ! optinal argument .true. starts FAM timers
  call start_timer(T_fam)


  !-----------------------------------------------------------------------------
  ! Read input from STDIN
  call ReadInput()

  !-----------------------------------------------------------------------------
  ! Initalize the matrices for performing derivatives on the mesh
  call inilag()

  !------------------------------------------------------------------------------
  ! Read all information from a .wf file
  call ReadWavefunction()
  !------------------------------------------------------------------------------
  ! Print all relevant input gleaned from STDIN and the wf file.
  call PrintInput()

  ! Provide memory for the derivatives of the spwfs
  call allocate_memory_derivatives(PairingType)

  ! Future dev: required for HFB
  ! ifail = 0
  ! call SolvePairing(pairingscheme, ifail)

  ! Derive all single-particle wavefunctions on the mesh
  if(store_derivatives) call deriveHF()

  !--------------------------------------------------------------------------------
  ! Step 0a: build explicitly the matrix of the single-particle hamiltonian and
  !          diagonalize it within the subspace spanned by the spwfs read from file
  Density     = densit(rho_can, kappa_pairing)
  Potentials  = calcPotentials(Density)
  sphamil     = Calc_Sphamil(potentials, .true.)

  ! ATTENTION: this explicit diagonalisation can break the apparent agreement
  !            between proton and neutron matices since the LAPACK diagonalisation
  !            might perform different rotations of the spwfs dependent on small
  !            numerical details.
  ! TODO: reenable once visual inspections are no longer necessary.

  !  call apply_subspace_rotation(sphamil, HFTransfo, spenergies)
  !  ! diagonalisation done; now recalculate other quantities
  !  if(store_derivatives) call deriveHF() ! and update derivatives
  !  Density     = densit(rho_can, kappa_pairing)
  !  Potentials  = calcPotentials(Density)
  !  sphamil     = Calc_Sphamil(potentials, .true.)

  ! Note: there is a silent assumption here that the HF-spectrum is sufficiently
  !       well-converged such that an explicit orthonormalisation will not change
  !       our mean-field state in any meaningful way. In the future, we might want
  !       to resolve the whole "pairing subproblem" again here and check that the
  !       structure does not vary too much.


  !----------------------------------------------------------------------------------
  ! Step 0b: calculate all relevant quantities on the meanfield level to enable a
  !          complete printout
  call update_spwf_properties_HF () !
  if(PairingType.eq.2) call update_spwf_properties_CAN()
  print_adv_spwf_properties = .true.
  call setBelyaevProcedure()
  call CalcEnergy(Density,Potentials,.true.) ! expensive parts included
  call calc_avg_gap()
  call full_printout(0,.false.,print_adv_spwf_properties)



  !---------------------------------------------------------------------------------
  ! construct the full HF densities rather than the merely the vector rho_can
  if (pairingtype .eq. 0) call iniHFdensities()

  ! call test_gmres_affine()
  ! stop

  !---------------------------------------------------------------------------------
  ! allocate the single-particle hamitonians 

  if(.not. allocated(dH_flat)) then
    allocate(dH_flat(nwt*nwt))
  endif

  if(.not. allocated(dH_flat_next)) then
    allocate(dH_flat_next(nwt*nwt))
  endif

  !---------------------------------------------------------------------------------
  ! solving FAM for a range of omega frequencies

  omega_num = int((omega_max - omega_min) / omega_step) + 1

  allocate(omega_arr(omega_num))
  allocate(S_arr(omega_num))
  allocate(S_free_arr(omega_num))
  allocate(iter_arr(omega_num))
  iter_arr = 0

  omega_curr = omega_min

  do omega_index=1, omega_num

    print 11, omega_curr

    !-------------------------------------------------------------------------------
    ! initialise FAM matrices end set perturbing external field
    !-------------------------------------------------------------------------------

    call inifam(omega_curr, Density, Potentials)

    S_free_arr(omega_index) = calc_strength()

    is_converged = .false.
    is_divergent = .false.

    if (fam_mixingscheme == 0) then

      !---------------------------------------------------------------------------------
      ! via GMRES on implicit matrix*vector procedure one_minus_T()
      !---------------------------------------------------------------------------------

      call alloc_gmres(one_minus_T, dH_free_flat, fam_maxiter, fam_maxiter, 1e-6_dp, size(dH_free_flat, 1), norm_dH, ScProd_dH)
      
      fam_verbose = 0

      ! initiliase the GMRES solver, using the free response as the initial guess x0
      call init_gmres(dH_free_flat)

      do iteration=1, gmres_itmax
        call iterate_gmres()

        if (gmres_res < gmres_tol) then 
          print 1
          print *, "Hooray! GMRES is converged! "
          iter_arr(omega_index) = iteration
          print 1
          exit
        endif

        if (iteration == gmres_itmax) then
          print 1
          print 1
          print *, "   Reached maximal number of iterations, ", fam_maxiter
          iter_arr(omega_index) = -gmres_itmax
        endif

      enddo

      print *, "One final FAM iteration based on GMRES solution:  "
      fam_verbose = 1
      call iterate_dHsp(x_gmres, dH_flat_next)
      print *, "Convergence check : || FAM(dH) - dH || / ||dH|| = ", norm_dH(dH_flat_next - x_gmres) / norm_dH(x_gmres)

    else if (fam_mixingscheme == 1) then

      !---------------------------------------------------------------------------------
      ! linear mixing while employing iterate_dHsp()
      !---------------------------------------------------------------------------------

      ! initialise the sp hamiltonians to the ones of the free response 
      dH_flat = dH_free_flat
      dH_flat_next = 0

      ! Start of the iterations 
      do iteration=1, fam_maxiter

        print 1
        print 2, iteration

        ! iterate the single-particle Hamiltonian by one complete FAM loop dH -> T(dH) + dH_free
        call iterate_dHsp(dH_flat, dH_flat_next)

        ! Run all kinds of unit tests; should be made optional as this includes a stop statement
        ! call run_FAM_tests(X,Y)

        ! simple linear mixing of sp hamiltonians dH[i+1] = a * dH[i+1] + (1-a) * dH[i]
        dH_flat_next = fam_lin_mix * dH_flat_next + (1.0_dp - fam_lin_mix) * dH_flat

        ! shift dH to prepare for the next iteration
        dH_flat = dH_flat_next

        !---------------------------------------------------------------------------------
        ! test convergenence

        ! Exit the loop if convergence is achieved.
        if (iteration > 1) then ! at least two iterations to be able to compare
         call test_convergence(is_converged, is_divergent)
          if(is_converged) then
            print 1
            print 1
            print *, "   Hooray! FAM is converged! "
            iter_arr(omega_index) = iteration
            exit
          endif
          if(is_divergent) then
            print 1
            print 1
            print *, "   FAM diverges, exiting"
            iter_arr(omega_index) = -iteration
            exit
          endif
        endif
        if (iteration == fam_maxiter) then
          print 1
          print 1
          print *, "   Reached maximal number of iterations, ", fam_maxiter
          iter_arr(omega_index) = -fam_maxiter
        endif
      enddo

      ! fixed-point check
      call iterate_dHsp(dH_flat, dH_flat_next)
      print *, "Convergence check : || FAM(dH) - dH || / ||dH|| = ", norm_dH(dH_flat_next - dH_flat) / norm_dH(dH_flat)

    endif


    !---------------------------------------------------------------------------------
    ! store the converged strength
    !---------------------------------------------------------------------------------
      
    omega_arr(omega_index) = omega_curr
    S_arr(omega_index) = calc_strength()

    print 1
    print *, "   number of iterations: ", iteration
    print *, "   converged strength:  "
    print *, "          l, m  = ", l, m
    print *, "          omega = ", omega_arr(omega_index)
    print *, "          S     = ", S_arr(omega_index) 
    print 1

    omega_curr = omega_curr + omega_step
    call dealloc_gmres()

  enddo

  if(fam_maxiter.eq.0) then
    write (famfilename, fmt='(a2,2i1,a10)') "S_", l, m, "_unper.fam"
  else 
    write (famfilename, fmt='(a2,2i1,a4)') "S_", l, m, ".fam"
  endif

  call write_fam_strength(omega_arr, S_arr, iter_arr, S_free_arr, l, m, famfilename)

  print *, "Reached the end successfully" 

  call stop_timer(T_fam)
  call print_all_timers()


  ! end of one FAM calculation;

end program run_FAM
