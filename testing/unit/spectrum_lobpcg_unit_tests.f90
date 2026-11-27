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
!  - preconditioning: the same solve with a kinetic CG preconditioner
!    M^{-1} = (b_k + T)^{-1} (per state, b_k = |theta_k| exact, T the
!    kinetic stencil) must converge in fewer iterations
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
  integer               :: dim
  integer               :: iters_plain, iters_precond
  real(dp), allocatable :: mat(:,:), work(:), eigenvalues(:)
  real(dp), allocatable :: x(:,:,:), thetas(:), residuals(:)
  real(dp), allocatable :: x0(:,:,:), thetas2(:), residuals2(:)
  real(dp)              :: dx, tol, err, dev, dv
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

  allocate(x(mv,4,n), thetas(n), residuals(n))
  allocate(x0(mv,4,n), thetas2(n), residuals2(n))

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

  ! Unpreconditioned solve.
  x0 = x
  call lobpcg_solve(apply_test_h, mv, n, dv, x, tol, thetas, residuals,      &
       &            iterations = iters_plain)

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

  ! Preconditioned solve: same initial guesses (x0 saved before the first
  ! call), M^{-1} = (b_k + T)^{-1} per state, b_k = |theta_k| exact. Must
  ! converge to the same spectrum, in strictly fewer iterations.
  x0 = x
  call lobpcg_solve(apply_test_h, mv, n, dv, x0, tol, thetas2, residuals2,   &
       &            apply_m = apply_test_m, iterations = iters_precond)
  do i=1,n
     err = abs(thetas2(i) - eigenvalues(i))
     if(err > 10.0_dp*tol) then
        print '(a,i0,a,es10.2,a,es10.2)', ' FAIL: preconditioned Ritz ', i,  &
             &                ' deviates from the exact eigenvalue by ', err
        stop 1
     endif
     if(residuals2(i) > tol) then
        print *, 'FAIL: preconditioned residual ', i, ' above the tolerance'
        stop 1
     endif
     dev = abs(sum(x0(:,:,i)**2) - 1.0_dp/dv)
     if(dev > 1.0e-10_dp) then
        print *, 'FAIL: preconditioned eigenvectors are not orthonormal'
        stop 1
     endif
  enddo
  if(iters_precond >= iters_plain) then
     print '(a,i0,a,i0)', ' FAIL: preconditioned solve took ', iters_precond,&
          & ' iterations, not fewer than the plain solve''s ', iters_plain
     stop 1
  endif
  print '(a,i0,a,i0,a)', ' Preconditioned solve: ', iters_precond,           &
       &   ' iterations vs ', iters_plain, ' without (same spectrum, tol ok)'

contains

  !---------------------------------------------------------------------------
  ! h acting on one mesh spinor of the current block: the 3-point periodic
  ! Laplacian stencil with the block's component shift.
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
  ! M^{-1} acting on one residual column: per state k, a conjugate-gradient
  ! solve of (b_k I + T) p = r, b_k = |theta_k| (the exact eigenvalue of that
  ! state, known here), T the shift-free kinetic stencil of apply_test_h.
  ! The operator is positive definite (T >= 0, b_k > 0); the solve damps the
  ! high-mesh-frequency content of the residual.
  !---------------------------------------------------------------------------
  subroutine apply_test_m(mvv, r, k, p)
    use compilation, only : dp
    integer, intent(in)     :: mvv, k
    real(dp), intent(in)    :: r(mvv,4)
    real(dp), intent(out)   :: p(mvv,4)

    real(dp) :: res(mvv), dir(mvv), tdir(mvv)
    real(dp) :: initnorm, newnorm, oldnorm, acg, bcg, b_k
    integer  :: c, iter

    b_k = abs(eigenvalues(k))
    p = 0.0_dp
    do c=1,4
       res      = r(:,c)
       initnorm = sum(res**2)*dv
       if(initnorm < 1.0e-30_dp) cycle
       newnorm  = initnorm
       dir      = res
       do iter=1,300
          call apply_test_T(mvv, dir, tdir)
          tdir = b_k*dir + tdir
          acg  = newnorm/(sum(dir*tdir)*dv)
          p(:,c)   = p(:,c)   + acg*dir
          res      = res      - acg*tdir
          oldnorm  = newnorm
          newnorm  = sum(res**2)*dv
          bcg      = newnorm/oldnorm
          dir      = res + bcg*dir
          if(newnorm < 1.0e-10_dp*initnorm) exit
       enddo
    enddo

  end subroutine apply_test_m

  !---------------------------------------------------------------------------
  ! T = the shift-free kinetic stencil: (6I - neighbours)/dx^2.
  !---------------------------------------------------------------------------
  subroutine apply_test_T(mvv, v, tv)
    use compilation, only : dp
    integer, intent(in)     :: mvv
    real(dp), intent(in)    :: v(mvv)
    real(dp), intent(out)   :: tv(mvv)

    integer :: i, j, k, l

    do k=1,nx
    do j=1,nx
    do i=1,nx
       l = lidx(i, j, k)
       tv(l) = ( 6.0_dp*v(l)                                                &
            &                     - v(lidx(mod(i,nx)+1, j, k))              &
            &                     - v(lidx(mod(i+nx-2, nx)+1, j, k))        &
            &                     - v(lidx(i, mod(j,nx)+1, k))              &
            &                     - v(lidx(i, mod(j+nx-2, nx)+1, k))        &
            &                     - v(lidx(i, j, mod(k,nx)+1))              &
            &                     - v(lidx(i, j, mod(k+nx-2, nx)+1))        &
            &              ) / (dx*dx)
    enddo
    enddo
    enddo
  end subroutine apply_test_T

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
