program run_FAM

  use compilation
  use IO
  use Tantalus, only : print_header, initialize_all_timers
  use fam
  use fam_testing, only : run_FAM_tests
  use gmres, only : do_gmres


  implicit none
  integer :: iteration
  logical :: is_converged, is_divergent
  real(kind=dp)  :: lin_mix_coeff=1.0d-2
  real(kind=dp) :: omega_curr
  integer :: omega_num, omega_index
  real(kind=dp), allocatable :: omega_arr(:), S_arr(:)
  integer, allocatable :: iter_arr(:)
  character(len=100) :: famfilename
  integer :: i, B, si,N

  ! integer :: ifail ! Future dev: required for HFB

  ! Print a nice header with all kinds of relevant info
  call print_header(.true.)

  !------------------------------------------------------------------------------
  ! starting all timers
  ! 
  ! -> This is necessary since subroutines below make use of the timers
  ! 
  call initialize_all_timers

  !-----------------------------------------------------------------------------
  ! Read input from STDIN
  ! 
  ! For FAMQRPA, the code should read in addition:
  ! 
  ! -  the type of perturbing operator/external field: E1, E2, M1, M2, ...
  !    and more complicated stuff when targetting beta-decay
  !    Important note: we will need to distinguish
  ! -  the frequency \omega_fam of the perturbing field
  ! -  the 'size' of the perturbation to perform the finite differencing
  ! -  a smearing parameter to avoid discontinuities at the poles of the 
  !    response function
  ! 
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
  ! Step 0: build explicitly the matrix of the single-particle hamiltonian and
  !         diagonalize it within the subspace spanned by the spwfs read from file
  Density     = densit(rho_can, kappa_pairing)
  Potentials  = calcPotentials(Density)
  sphamil     = Calc_Sphamil(potentials, .true.)

  ! ! Testing printout - left for now
  ! print *, 'BEFORE'
  ! si = 0
  ! do B=1,8
  !   N = HFBLocks(B)

  !   print *, 'BLOCK', B, N
  !   do i=si+1,si+N
  !     print ('(99f10.3)'), sphamil(i, si+1:si+N)
  !   enddo
  !   print *
  !   si = si + N
  ! enddo

  call apply_subspace_rotation(sphamil, HFTransfo, spenergies)
  ! diagonalisation done; now recalculate other quantities
  if(store_derivatives) call deriveHF() ! and update derivatives
  Density     = densit(rho_can, kappa_pairing)
  Potentials  = calcPotentials(Density)
  sphamil     = Calc_Sphamil(potentials, .true.)

  ! ! Testing printout - left for now
  ! print *, 'AFTER'
  ! si = 0
  ! do B=1,8
  !   N = HFBLocks(B)

  !   print *, 'BLOCK', B, N
  !   do i=si+1,si+N
  !     print ('(99f10.3)'), sphamil(i, si+1:si+N)
  !   enddo
  !   print *
  !   si = si + N
  ! enddo
  ! stop

  ! print ('(99f10.3)'), rho_can(:)

  !
  ! Note: there is a silent assumption here that the HF-spectrum is sufficiently
  !       well-converged such that an explicit orthonormalisation will not change
  !       our mean-field state in any meaningful way. In the future, we might want
  !       to resolve the whole "pairing subproblem" again here and check that the
  !       structure does not vary too much.


  !---------------------------------------------------------------------------------
  ! construct the full HF densities rather than the merely the vector rho_can
  if (pairingtype .eq. 0) call iniHFdensities()


  !---------------------------------------------------------------------------------
  ! solving FAM for a range of omega frequencies

  print *, omega_max, omega_min,  omega_step

  omega_num = int((omega_max - omega_min) / omega_step) + 1

  allocate(omega_arr(omega_num))
  allocate(S_arr(omega_num))
  allocate(iter_arr(omega_num))
  iter_arr = 0

  omega_curr = omega_min

  ! maxfamiter = 1000

  do omega_index=1, omega_num

    !-------------------------------------------------------------------------------
    ! initialise FAM matrices end set perturbing external field
    call inifam(omega_curr, Potentials)

    if( calc_strength() .ge. 0.1) then
      lin_mix_coeff=1.0d-3
    else
      if( calc_strength() .ge. 0.001) then
        lin_mix_coeff=1.0d-2
      else
        lin_mix_coeff=1.0d-1
      endif
    endif

    ! Run all kinds of unit tests; should be made optional as this includes a stop statement
    ! call run_FAM_tests(X,Y)

    is_converged = .false.
    is_divergent = .false.

    ! Start of the iterations 
    do iteration=1, maxfamiter

      print *, "FAM iteration : ", iteration

      ! build the perturbed hamiltonian using explicit linearisation of the field
      call build_dH_explicit(Density, dRs, dRa)

      ! calculate X and Y from the perturbed sp Hamil dH
      call calculate_XY(dH)
      
      ! Apply simple linear mixing of X and Y. 
      call mix_XY_linear(lin_mix_coeff)
      ! To be replaced with something more fancy in the future

      ! build the perturbed densities on the mesh dRs, dRa from X and Y
      call build_perturbed_densities(X, Y, dRs, dRa)

      call store_XY_hist()

      ! call print_all_fam_spmat()

      print *, " S(", omega_curr, ") = ",  calc_strength()

      ! Exit the loop if convergence is achieved.
      if (iteration > 1) then ! at least two iterations to be able to compare
       call test_convergence(is_converged, is_divergent)
        if(is_converged) then
          print *, "Hooray! FAM is converged! "
          ! omega_arr(omega_index) = omega_curr
          ! S_arr(omega_index) = calc_strength()
          iter_arr(omega_index) = iteration
          exit
        endif
        if(is_divergent) then
          print *, "FAM diverges, exiting"
          ! omega_arr(omega_index) = omega_curr
          ! S_arr(omega_index) = calc_strength()
          iter_arr(omega_index) = -iteration
          exit
        endif
      endif
      if (iteration == maxfamiter) then
        print *, "Reached maximal number of iterations, ", maxfamiter
        iter_arr(omega_index) = -maxfamiter
      endif
    enddo

    omega_arr(omega_index) = omega_curr
    S_arr(omega_index) = calc_strength()


    print *, " S(", omega_arr(omega_index), ") = ", S_arr(omega_index)


    omega_curr = omega_curr + omega_step

  enddo

  if(maxfamiter.eq.0) then
    write (famfilename, fmt='(a2,2i1,a10)') "S_", l, m, "_unper.fam"
  else 
    write (famfilename, fmt='(a2,2i1,a4)') "S_", l, m, ".fam"
  endif

  print *, omega_arr
  print *, S_arr
  print *, iter_arr
  call write_fam_strength(omega_arr, S_arr, iter_arr, l, m, famfilename)

  print *, "Reached the end successfully" 


  ! end of one FAM calculation;

end program run_FAM
