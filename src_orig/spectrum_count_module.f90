!===============================================================================
! Stochastic count of the eigenvalues of a mesh single-particle hamiltonian
! below a cutoff, by Chebyshev trace estimation. Sizes the fixed-block LOBPCG
! subspace of spectrum runs; see development item 007.
!
! The count is a property of the DISCRETIZED hamiltonian: the cutoff
! comparison lives on the mesh, as trim_spectrum applies it. The application
! of h is supplied as a procedure argument, so this module is self-contained:
! the spectrum driver passes apply_sphamil, the unit test a standalone
! kinetic operator. MPI is the caller's business: ranks holding disjoint
! probes may sum their partial moments before combining them here.
!
! Method: n_hat = (1/S) sum_s v_s^T p(h) v_s with random +-1 probes v_s and
! p the Chebyshev series of the indicator of [lower, cutoff]. The expansion
! lives on the FULL spectral interval [lower, upper] - upper must be the
! TRUE spectral maximum of h, not the cutoff: eigenvalues above the cutoff
! sit inside the recurrence's stability region only if upper covers them.
!
! Damping: the series is Jackson-damped. The Jackson kernel is positive,
! so the damped step lies in [0,1] everywhere - no Gibbs overshoot - and
! the estimate deviates from the exact count only by the leakage across
! the step, confined to the kernel width pi/(degree+1) in theta. On the
! mesh that width is far below the level spacing at the cutoff, so the
! bias is negligible against the caller's overallocation; it is NOT a
! rigorous one-sided (under-) approximation - do not pretend otherwise.
!
! Caveats:
!  - n_hat is an expectation value; its variance scales as 1/S. The
!    caller's overallocation (one integer) covers it.
!  - The damping bias is O(rho/degree) around the step (rho = level
!    density); the degree must stay well above the level density in the
!    window.
!===============================================================================
module spectrum_count_module

  use compilation, only : dp
  use geninfo,     only : stp

  implicit none

  private
  public :: chebyshev_count

  integer, parameter    :: i64 = selected_int_kind(18)
  integer(i64), private :: count_seed = 0

  abstract interface
     !-------------------------------------------------------------------------
     ! Apply h to m mesh spinors, each of the form (mv,4); common symmetry
     ! properties per block and isospin are the caller's business.
     !-------------------------------------------------------------------------
     function apply_h_interface(mv, x) result(hx)
       import dp
       integer, intent(in)  :: mv
       real(dp), intent(in) :: x(mv,4)
       real(dp)             :: hx(mv,4)
     end function apply_h_interface
  end interface

contains

  !-----------------------------------------------------------------------------
  ! Counts the eigenvalues of h below cutoff, from probe trace estimation
  ! with a Jackson-damped step polynomial.
  !-----------------------------------------------------------------------------
  ! Input:
  !   apply_h - applies h to a single mesh spinor (see apply_h_interface)
  !   mv      - mesh-spinor dimension of the probes
  !   lower   - lower spectral bound of h (below its lowest eigenvalue)
  !   upper   - upper spectral bound: the TRUE maximum of h, not the cutoff
  !   cutoff  - the step of the indicator; eigenvalues above it count zero
  !   probes  - number of random +-1 probe vectors
  !   degree  - Chebyshev degree of the damped step polynomial
  ! Output:
  !   n_hat   - the stochastic count estimate; ceiling it for the
  !             subspace size and add the caller's overallocation
  !-----------------------------------------------------------------------------
  function chebyshev_count(apply_h, mv, lower, upper, cutoff, probes,          &
       &                                                   degree) result(n_hat)

    procedure(apply_h_interface)      :: apply_h
    integer, intent(in)               :: mv, probes, degree
    real(dp), intent(in)              :: lower, upper, cutoff
    real(dp)                          :: n_hat

    real(dp)       :: v(mv,4), u0(mv,4), u1(mv,4), u2(mv,4)
    real(dp)       :: t(0:degree), b(0:degree)
    real(dp)       :: scale, shift, theta_c, pi
    integer(i64)   :: stream(mv*4)
    integer        :: s, j

    pi    = acos(-1.0_dp)
    scale = 2.0_dp / (upper - lower)
    shift = -(upper + lower) / (upper - lower)
    ! Mapped step; theta_c = arccos(step), with h = lower at theta = pi.
    theta_c = acos(cutoff*scale + shift)
    if (theta_c <= 0.0_dp .or. theta_c >= pi)                            &
         &  call stp('chebyshev_count: cutoff outside the spectral interval.', &
         &          'chebyshev_count')

    call step_coefficients(theta_c, degree, b)

    ! One recurrence pass per probe, accumulating the Chebyshev moments
    ! t_j = (1/S) sum_s v_s^T T_j(h_scaled) v_s against the FIXED probe v.
    t = 0.0_dp
    do s=1,probes
       ! Rademacher probe, seeded per probe from a fixed base: the s-th probe
       ! is independent of the probe count, so runs are reproducible and
       ! golden_diff-friendly.
       count_seed = ishft(1_i64, 20) + s
       call next_probe_entries(mv, stream)
       v  = reshape(stream, (/ mv, 4 /))
       u0 = v

       ! Three-term recurrence u_0 = v, u_1 = h_scaled u_0,
       ! u_{j+1} = 2 h_scaled u_j - u_{j-1}: one application of h per degree.
       u1 = scale*apply_h(mv, u0) + shift*u0
       t(0) = t(0) + sum(v*v)
       t(1) = t(1) + sum(v*u1)
       do j=1,degree-1
          u2 = 2.0_dp*(scale*apply_h(mv, u1) + shift*u1) - u0
          t(j+1) = t(j+1) + sum(v*u2)
          u0 = u1 ; u1 = u2
       enddo
    enddo
    t = t / real(probes, dp)

    n_hat = 0.5_dp*b(0)*t(0)
    do j=1,degree
       n_hat = n_hat + b(j)*t(j)
    enddo

    print '(a,i0,a,f10.3)', ' Chebyshev count: degree ', degree,            &
         &                  ', estimate ', n_hat

  end function chebyshev_count

  !-----------------------------------------------------------------------------
  ! Chebyshev coefficients b_j of the Jackson-damped indicator of
  ! theta in [theta_c, pi] (h <= cutoff). The Jackson kernel is positive,
  ! so the damped polynomial lies in [0,1] everywhere: no Gibbs overshoot.
  ! Its error is leakage across the step, confined to the kernel width
  ! ~ pi/(degree+1) in theta; eigenvalue spacing on the mesh is orders of
  ! magnitude coarser, so the count estimate deviates from the exact count
  ! by far less than the caller's overallocation.
  !-----------------------------------------------------------------------------
  subroutine step_coefficients(theta_c, degree, b)
    real(dp), intent(in)  :: theta_c
    integer, intent(in)   :: degree
    real(dp), intent(out) :: b(0:degree)

    real(dp) :: pi, cot, jackson
    integer  :: j

    pi  = acos(-1.0_dp)
    cot = 1.0_dp / tan(pi/real(degree + 1, dp))

    ! Exact Chebyshev coefficients of the step, times the Jackson kernel
    ! factors (Weisse, Wellein, Alvermann, Fehske 2006): g_0 = 1, g_degree = 0.
    b(0) = 2.0_dp*(1.0_dp - theta_c/pi)
    do j=1,degree
       jackson = (real(degree - j + 1, dp)*cos(pi*real(j, dp)              &
       &                  /real(degree + 1, dp))                           &
       &                  + sin(pi*real(j, dp)/real(degree + 1, dp))*cot)  &
       &                  /real(degree + 1, dp)
       b(j) = -2.0_dp*sin(real(j, dp)*theta_c)/(pi*real(j, dp)) * jackson
    enddo

  end subroutine step_coefficients

  !-----------------------------------------------------------------------------
  ! Draws mv*4 independent Rademacher entries from the module-level xorshift64
  ! stream.
  !-----------------------------------------------------------------------------
  subroutine next_probe_entries(mv, stream)
    integer, intent(in)         :: mv
    integer(i64), intent(out)   :: stream(mv*4)

    integer :: i

    do i=1,mv*4
       count_seed = ieor(count_seed, ishft(count_seed, 13))
       count_seed = ieor(count_seed, ishft(count_seed, -7))
       count_seed = ieor(count_seed, ishft(count_seed, 17))
       stream(i) = merge(1_i64, -1_i64, btest(count_seed, 0))
    enddo

  end subroutine next_probe_entries

end module spectrum_count_module
