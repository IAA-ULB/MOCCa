program run_FAM

  use compilation
  use IO
  use version,  only : print_header
  use Tantalus, only : initialize_all_timers, full_printout
  use Tantalus, only : update_spwf_properties_HF, update_spwf_properties_CAN
  use fam
  use fam_testing, only : run_FAM_tests, test_gmres, test_gmres_affine, test_linearity_T, test_linearity_FAM_coulomb
  use gmres 
  use timing

  1 format(86('-'))
  11 format(/,24('='), ' omega = ', f5.2, ' MeV ', 24('='),/)
  2 format('FAM iteration = ', i5) 
  3 format(' S_',i1,i1,' (', f5.2, ') = ', es10.3)

  implicit none
  integer :: iter, num_iter, ifail
  logical :: is_converged, is_divergent
  real(kind=dp) :: omega_curr
  integer :: omega_num, omega_index
  real(kind=dp) :: strength_free

  complex(KIND=dp) :: S_complex_decomp(8) = 0
  real(KIND=dp) :: S_decomp(8) = 0

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
  call SolvePairing(pairingscheme, ifail)

  ! Derive all single-particle wavefunctions on the mesh
  if(store_derivatives) call deriveHF()

  !--------------------------------------------------------------------------------
  ! Step 0a: build explicitly the matrix of the single-particle hamiltonian and
  !          diagonalize it within the subspace spanned by the spwfs read from file
  Density     = densit(rho_can, kappa_pairing)
  call CalculateMoments(Density,.true.)           ! necessary here if constraints are included
  Potentials  = calcPotentials(Density)
  sphamil     = Calc_Sphamil(potentials, .true.)

  ! ATTENTION: this explicit diagonalisation can break the apparent agreement
  !            between proton and neutron matices since the LAPACK diagonalisation
  !            might perform different rotations of the spwfs dependent on small
  !            numerical details.
  dispersions = calculate_spwf_dispersions(potentials)
  call apply_subspace_rotation(sphamil, HFTransfo, spenergies)
  if(store_derivatives) call deriveHF() ! and update derivatives
  ! Explicitly recalculate dispersion to provide an idea of the quality of the mean-field state
  dispersions = calculate_spwf_dispersions(potentials)
  ! diagonalisation done; now recalculate other quantities
  call SolvePairing(pairingscheme, ifail)
  Density     = densit(rho_can, kappa_pairing)
  Potentials  = calcPotentials(Density)
  sphamil     = Calc_Sphamil(potentials, .true.)
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

  !---------------------------------------------------------------------------------
  ! Evaluate the energy weighted sum rule
  ewsr = calc_EWSR()


  !---------------------------------------------------------------------------------
  ! create the FAM output file
  call init_fam_file_new(famfile)

  if(xyfile .ne. '') then
    call init_xy_file(xyfile)
  endif

  if(DENFILE .ne. '') then
    call init_perturbed_denfile(DENFILE)
  endif



  !---------------------------------------------------------------------------------
  ! allocate the single-particle hamiltonians 
  if(.not. allocated(dH_flat)) then
    allocate(dH_flat(nwt*nwt))
  endif

  if(.not. allocated(dH_flat_next)) then
    allocate(dH_flat_next(nwt*nwt))
  endif

  !---------------------------------------------------------------------------------
  ! solving FAM for a range of omega frequencies

  omega_num = int((omega_max - omega_min) / omega_step) + 1

  omega_curr = omega_min

  do omega_index=1, omega_num

    print 11, omega_curr

    !-------------------------------------------------------------------------------
    ! initialise FAM matrices end set perturbing external field
    !-------------------------------------------------------------------------------

    call inifam(omega_curr, Density, Potentials)

    strength_free = calc_strength()

    is_converged = .false.
    is_divergent = .false.

    if (fam_mixingscheme == 0) then

      !---------------------------------------------------------------------------------
      ! via GMRES on implicit matrix*vector procedure one_minus_T()
      !---------------------------------------------------------------------------------

      call alloc_gmres(one_minus_T, dH_free_flat, fam_maxiter, fam_maxhist, fam_precision, norm_dH, ScProd_dH)
      
      fam_verbose = 0

      ! initiliase the GMRES solver, using the free response as the initial guess x0
      call init_gmres(dH_free_flat)

      do iter=1, gmres_itermax
        call iterate_gmres()

        if (gmres_res < gmres_precision) then 
          print 1
          print *, "Hooray! GMRES is converged! "
          num_iter = iter
          print 1
          exit
        endif

        if (iter == gmres_itermax) then
          print 1
          print 1
          print *, "   Reached maximal number of iterations, ", fam_maxiter
          num_iter = - gmres_itermax
        endif

      enddo

      call extract_x_gmres()

      print *, "One final FAM iteration based on GMRES solution:  "
      fam_verbose = 1
      call iterate_dHsp(x_gmres, dH_flat_next)
      print *, "Convergence check : || FAM(dH) - dH || / ||dH|| = ", norm_dH(dH_flat_next - x_gmres) / norm_dH(x_gmres)

      call dealloc_gmres()


    else if (fam_mixingscheme == 1) then

      !---------------------------------------------------------------------------------
      ! linear mixing while employing iterate_dHsp()
      !---------------------------------------------------------------------------------

      ! initialise the sp hamiltonians to the ones of the free response 
      dH_flat = dH_free_flat
      dH_flat_next = 0

      ! Start of the iterations 
      do iter=1, fam_maxiter

        print 1
        print 2, iter

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
        if (iter > 1) then ! at least two iterations to be able to compare
         call test_convergence(is_converged, is_divergent)
          if(is_converged) then
            print 1
            print 1
            print *, "   Hooray! FAM is converged! "
            num_iter = iter
            exit
          endif
          if(is_divergent) then
            print 1
            print 1
            print *, "   FAM diverges, exiting"
            num_iter = - iter
            exit
          endif
        endif
        if (iter == fam_maxiter) then
          print 1
          print 1
          print *, "   Reached maximal number of iterations, ", fam_maxiter
          num_iter = - fam_maxiter
        endif
      enddo

      ! fixed-point check
      call iterate_dHsp(dH_flat, dH_flat_next)
      print *, "Convergence check : || FAM(dH) - dH || / ||dH|| = ", norm_dH(dH_flat_next - dH_flat) / norm_dH(dH_flat)

    endif


    !---------------------------------------------------------------------------------
    ! store the converged strength
    !---------------------------------------------------------------------------------
      
    strength = calc_strength()

    print 1
    print *, "   number of iterations: ", num_iter
    print *, "   converged strength:  "
    print *, "          l, m  = ", l, m
    print *, "          omega = ", omega_curr
    print *, "          S     = ", strength 
    print 1

    call calc_strength_decomp(S_complex_decomp, S_decomp)

    ! call append_fam_file(omega_curr, strength, num_iter, strength_free, famfile)

    call append_fam_file_new(S_decomp, num_iter, famfile)

    if(xyfile .ne. '') then
      call append_xy_file(xyfile)
    endif


    if(DENFILE .ne. '') then
      call append_perturbed_denfile(dRs, dRa, DENFILE)
    endif

    omega_curr = omega_curr + omega_step

  enddo

  print *, "Reached the end successfully" 

  call stop_timer(T_fam)
  call print_all_timers()


  ! end of one FAM calculation;

end program run_FAM
