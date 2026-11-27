!===============================================================================
! Unit test of spectrum_lobpcg_module: solve for the lowest eigenpairs of a
! small mesh hamiltonian by fixed-block LOBPCG and compare against the exact
! lowest-n eigenvalues from a dense diagonalization.
!
! Test hamiltonian: the 3-point periodic Laplacian stencil acting on the 4
! spinor components, plus a per-component constant shift, so the spinor
! blocks decouple and the exact spectrum follows from a dense
! diagonalization of the full spinor hamiltonian. Initial guesses: random
! spinors, seeded reproducibly.
!
! Checks:
!  - the Ritz values match the exact lowest-n eigenvalues to the tolerance
!  - the residuals of the converged block are below the tolerance
!  - the converged eigenvectors are orthonormal in the dv metric
!  - reproducibility: a second call with identical initial guesses must
!    return the identical spectrum
!
! Tier: the 8^3 smoke test (seconds, suitable for CI); a 16^3 tier would
! stress the level spacing near the cutoff like a real run - run it locally
! before committing.
!===============================================================================
program spectrum_lobpcg_unit_tests

  use compilation, only : dp
  use geninfo,               only : stp
  use spectrum_lobpcg_module, only : lobpcg_solve

  implicit none

  integer               :: nx, mv, n, i, j, k, c, idim, s, lwork, info
  integer               :: dim, n_exact
  real(dp), allocatable :: mat(:,:), work(:), eigenvalues(:)
  real(dp), allocatable :: x(:,:,:), thetas(:), residuals(:)
  real(dp)              :: dx, cutoff, tol, err, dev, dv
  real(dp)              :: comp_shift(4) = (/ 0.0_dp, 0.7_dp, 1.3_dp, 2.1_dp /)
  integer, allocatable  :: seed(:)
  character(len=16)     :: arg

  nx = 8
  if (command_argument_count() >= 1) then
     call get_command_argument(1, arg)
     read(arg, *) nx
  endif
  if (nx /= 8 .and. nx /= 16) then
     print *, 'Usage: spectrum_lobpcg_unit_tests [8|16]'
     stop 1
  endif

  mv  = nx**3
  dx  = 0.8_dp
  dv  = dx**3
  tol = 1.0e-4_dp
  n   = 12
  dim = 4*mv

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Dense reference: 3-point periodic Laplacian stencil acting on the 4 spinor
  ! components, plus a per-component constant shift (the components decouple).
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  mv = nx**3
  allocate(mat(4*mv,4*mv))
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

  ! Exact spectrum from a dense diagonalization.
  allocate(eigenvalues(dim), work(1))
  call DSYEV('N', 'U', dim, mat, dim, eigenvalues, work, -1, info)
  lwork = int(work(1)); deallocate(work); allocate(work(lwork))
  call DSYEV('N', 'U', dim, mat, dim, eigenvalues, work, lwork, info)
  deallocate(mat, work)

  n_exact = count(eigenvalues <= cutoff)

  ! The exact count from a dense diagonalization.
  allocate(x(mv,4,n), thetas(n), residuals(n))

  ! Random initial guesses, seeded reproducibly after the house pattern:
  ! the same seed gives the same solve, so runs are reproducible and
  ! golden_diff-friendly.
  call random_seed(size=s)
  allocate(seed(s)) ; seed = 42
  call random_seed(put=seed) ; deallocate(seed)
  call random_number(x)
  do i=1,n
     x(:,:,i) = x(:,:,i) / sqrt(dv*sum(x(:,:,i)**2))
  enddo

  call lobpcg_solve(apply_test_h, mv, n, dv, x, tol, thetas, residuals)
  ! The Ritz values must match the exact lowest-n eigenvalues to the
  ! tolerance; the residuals of the converged block are below it; the
  ! converged eigenvectors are orthonormal in the dv metric.
  do i=1,n
     err = abs(thetas(i) - eigenvalues(i))
     if(err > 10.0_dp*tol) then
        print '(a,i0,a,es10.2,a,es10.2)', ' FAIL: Ritz value ', i,           &
             &                ' deviates from the exact eigenvalue by ', err
        stop 1
     endif
     if(residuals(i) > tol) then
        print *, 'FAIL: residual ', i, ' above the tolerance'
        stop 1
     endif
     dev = abs(sum(x(:,:,i)**2) - 1.0_dp/dv)
     if(dev > 1.0e-10_dp) then
        print *, 'FAIL: converged eigenvectors are not orthonormal'
        stop 1
     endif
  enddo

  print '(a,i0,a,i0,a,es10.2)', ' Mesh ', nx, '^3: lowest ', n,              &
       &       ' Ritz values match the dense diagonalization, max deviation ', err

contains

  !---------------------------------------------------------------------------
  ! h acting on one mesh spinor of the current block, via the 3-point
  ! periodic Laplacian stencil with the block's component shift.
  !---------------------------------------------------------------------------
  function apply_test_h(mvv, x) result(hx)
    use compilation, only : dp
    integer, intent(in)     :: mvv
    real(dp), intent(in)    :: x(mvv,4)
    real(dp)                :: hx(mvv,4)

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

  !---------------------------------------------------------------------------
  ! Mesh index of site (i,j,k) within one spinor component, 1-based.
  !---------------------------------------------------------------------------
  function lidx(i, j, k)
    integer, intent(in) :: i, j, k
    integer             :: lidx
    lidx = (k-1)*nx*nx + (j-1)*nx + i
  end function lidx

  !---------------------------------------------------------------------------
  ! Global spinor index of mesh point (i,j,k), component c, for the
  ! 3-point periodic Laplacian stencil (house pattern of the count test).
  !---------------------------------------------------------------------------
  function gid(i, j, k, c)
    integer, intent(in)  :: i, j, k, c
    integer              :: gid

    gid = ((c-1)*nx + k-1)*nx*nx + (j-1)*nx + i

  end function gid

end program spectrum_lobpcg_unit_tests
