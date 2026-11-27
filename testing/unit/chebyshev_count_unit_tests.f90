!===============================================================================
! Unit test of spectrum_count_module: count the eigenvalues of a small mesh
! hamiltonian by Chebyshev trace estimation and compare against the exact
! count from a dense diagonalization.
!
! Test hamiltonian: the 3-point periodic Laplacian stencil acting on the 4
! spinor components, plus a per-component constant shift, so the spinor
! blocks decouple and the exact count follows from a dense diagonalization
! of the full spinor hamiltonian. The spectral bounds handed to the counter
! come from Gershgorin of the stencil only, never from the dense solve.
!
! Tiers (argument):
!   8  - smoke test (seconds, suitable for CI)
!   16 - reference tier: level density near the cutoff comparable to a real
!        run; stresses probe variance and the Chebyshev approximation gap.
!        Too slow for CI; run locally before committing.
!===============================================================================
program chebyshev_count_unit_tests

  use compilation, only : dp
  use geninfo,               only : stp
  use spectrum_count_module, only : chebyshev_count

  implicit none

  external :: DSYEV

  integer               :: nx, mv, dim, i, j, k, c, idim, n_exact, lwork, info
  real(dp), allocatable :: mat(:,:), work(:), eigenvalues(:)
  real(dp)              :: dx, cutoff, lower, upper, n_hat, mu, nh
  real(dp)              :: comp_shift(4) = (/ 0.0_dp, 0.7_dp, 1.3_dp, 2.1_dp /)
  character(len=16)     :: arg

  nx = 8
  if (command_argument_count() >= 1) then
     call get_command_argument(1, arg)
     read(arg, *) nx
  endif
  if (nx /= 8 .and. nx /= 16) then
     print *, 'Usage: chebyshev_count_unit_tests [8|16]'
     stop 1
  endif

  mv  = nx**3
  dim = 4*mv
  dx  = 0.8_dp
  cutoff = 0.4_dp * (12.0_dp/(dx*dx)) + 1.0_dp

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Dense reference: 3-point periodic Laplacian stencil acting on the 4 spinor
  ! components, plus a per-component constant shift (the components decouple).
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  allocate(mat(dim,dim))
  mat = 0.0_dp
  do c=1,4
  do k=1,nx
  do j=1,nx
  do i=1,nx
     idim = gid(i, j, k, c)
     mat(idim,idim) = 6.0_dp/(dx*dx) + comp_shift(c)
     mat(idim, gid(mod(i,nx)+1, j, k, c)) = -1.0_dp/(dx*dx)
     mat(idim, gid(mod(i+nx-2, nx)+1, j, k, c)) = -1.0_dp/(dx*dx)
     mat(idim, gid(i, mod(j,nx)+1, k, c)) = -1.0_dp/(dx*dx)
     mat(idim, gid(i, mod(j+nx-2, nx)+1, k, c)) = -1.0_dp/(dx*dx)
     mat(idim, gid(i, j, mod(k,nx)+1, c)) = -1.0_dp/(dx*dx)
     mat(idim, gid(i, j, mod(k+nx-2, nx)+1, c)) = -1.0_dp/(dx*dx)
  enddo
  enddo
  enddo
  enddo

  ! Exact count from a dense diagonalization.
  allocate(eigenvalues(dim), work(1))
  call DSYEV('N', 'U', dim, mat, dim, eigenvalues, work, -1, info)
  lwork = int(work(1)); deallocate(work); allocate(work(lwork))
  call DSYEV('N', 'U', dim, mat, dim, eigenvalues, work, lwork, info)
  n_exact = count(eigenvalues <= cutoff)

  ! The bounds the counter sees must come from cheap information only:
  ! Gershgorin of the stencil, never from the dense diagonalization.
  lower = minval(comp_shift) - 0.1_dp
  upper = 12.0_dp/(dx*dx) + maxval(comp_shift)

  n_hat = chebyshev_count(apply_test_h, mv, lower, upper, cutoff, 32, 1024)

  print '(a,i0,a,i0,2(a,f8.3))', ' Mesh ', nx, '^3: exact count ', n_exact,    &
       &       ', Chebyshev estimate ', n_hat, ' at cutoff ', cutoff

  ! The estimate must not overshoot the exact count by more than one state,
  ! and must not miss more than a small fraction of the window.
  if (n_hat > real(n_exact, dp) + 0.5_dp) then
     print *, 'FAIL: Chebyshev count overshoots the exact count'
     stop 1
  endif
  if (n_hat < real(n_exact, dp) - 0.05_dp*real(n_exact, dp) - 2.0_dp) then
     print *, 'FAIL: Chebyshev count misses too many states of the window'
     stop 1
  endif

  ! Reproducibility: a second call with identical bounds must return the
  ! identical estimate.
  if (abs(chebyshev_count(apply_test_h, mv, lower, upper, cutoff, 32, 1024)      &
       &                   - n_hat) > 1.0e-12_dp) then
     print *, 'FAIL: Chebyshev count is not reproducible'
     stop 1
  endif

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Trivial-h scenario, checked before the dense tier: h = mu * identity on a
  ! small spinor. Every eigenvalue equals mu, so the counter's estimate is
  ! mu * p(mu_scaled) * spinor-size exactly - no diagonalization needed.
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mu = 2.0_dp
  nh = chebyshev_count(identity_h, 8, 0.0_dp, 10.0_dp, 3.0_dp, 8, 64)
  ! mu below the cutoff: full count; mu above it: zero.
  if (abs(nh - 32.0_dp) > 1.0_dp .or.                                 &
       &      abs(chebyshev_count(identity_h, 8, 0.0_dp, 10.0_dp, 1.0_dp,    &
       &                         8, 64)) > 1.0_dp) then
     print *, 'FAIL: trivial-h scenario'
     stop 1
  endif
  print *, 'trivial-h scenario PASS'

  print *, 'PASS'
  stop

contains

  !---------------------------------------------------------------------------
  ! The test hamiltonian: 3-point periodic Laplacian stencil acting on the
  ! 4 spinor components, plus a per-component constant shift.
  !---------------------------------------------------------------------------
  function apply_test_h(mv, x) result(hx)
    use compilation, only : dp
    integer, intent(in)     :: mv
    real(dp), intent(in)    :: x(mv,4)
    real(dp)                :: hx(mv,4)

    integer :: i, j, k, c, l

    do c=1,4
    do k=1,nx
    do j=1,nx
    do i=1,nx
       l = (k-1)*nx*nx + (j-1)*nx + i
       hx(l,c) = ( 6.0_dp*x(l,c)                                             &
            &                       - x(lidx(mod(i,nx)+1, j, k),c)           &
            &                       - x(lidx(mod(i+nx-2, nx)+1, j, k),c)     &
            &                       - x(lidx(i, mod(j,nx)+1, k),c)           &
            &                       - x(lidx(i, mod(j+nx-2, nx)+1, k),c)     &
            &                       - x(lidx(i, j, mod(k,nx)+1),c)           &
            &                       - x(lidx(i, j, mod(k+nx-2, nx)+1),c)     &
            &              ) / (dx*dx) + comp_shift(c)*x(l,c)
    enddo
    enddo
    enddo
    enddo
  end function apply_test_h

  ! Mesh index of site (i,j,k) within one spinor component, 1-based.
  function lidx(i, j, k)
    integer, intent(in) :: i, j, k
    integer             :: lidx
    lidx = (k-1)*nx*nx + (j-1)*nx + i
  end function lidx

  ! Global index of spinor component c at mesh site (i,j,k) in the dense
  ! reference matrix (component-major), 1-based.
  function gid(i, j, k, c)
    integer, intent(in) :: i, j, k, c
    integer             :: gid
    gid = ((c-1)*nx + k-1)*nx*nx + (j-1)*nx + i
  end function gid

  ! The trivial-h scenario's hamiltonian: h = mu * identity.
  function identity_h(mvv, x) result(hx)
    use compilation, only : dp
    integer, intent(in)     :: mvv
    real(dp), intent(in)    :: x(mvv,4)
    real(dp)                :: hx(mvv,4)
    hx = 2.0_dp * x
  end function identity_h

end program chebyshev_count_unit_tests
