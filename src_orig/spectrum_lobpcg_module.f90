!===============================================================================
! Fixed-block LOBPCG for the lowest eigenvalues of a mesh single-particle
! hamiltonian, per symmetry block of spectrum runs; see development item 007.
!
! The subspace has dimension exactly n (the block's spwf count, from the
! Chebyshev count) with the basis {X, R, P}: X the current eigenvector
! approximations, R = hX - X Theta the residuals, P the previous step
! directions. One dense projected eigensolve per iteration (DGEMM + DSYEV);
! the preconditioner M enters through an optional procedure argument applying
! M^{-1} to one residual column (column index passed along, so the caller can
! use per-state parameters). Without it the residual block carries raw R and
! convergence is locally linear (Knyazev 2001); with it, the block carries
! M^{-1} R and the rate improves accordingly. Convergence when all
! residuals ||h x - theta x|| are below the tolerance, uniformly over the
! block.
!
! The basis is orthonormalized (modified Gram-Schmidt over all its columns,
! in the dv-weighted mesh metric) BEFORE h acts on it, so the Ritz rotation
! keeps the new X dv-orthonormal and the residual identity
! ||hX||^2_dv - theta^2 is exact.
!
! A residual already below the tolerance contributes a nearly degenerate R
! column that MGS would renormalize from round-off noise; such a column is
! replaced by a fresh random spinor to keep the basis full rank.
!
! The application of h is supplied as a procedure argument, so this module
! is self-contained: the spectrum driver passes apply_sphamil with the
! block's reflection signs and the frozen potentials, a unit test a
! standalone operator. h preserves the block's symmetry sector by
! construction, so the projected matrix is symmetric and every basis
! vector stays in-sector - no symmetrization of basis states.
!
! Serial: spectrum runs are not MPI-ready (item 007); the caller keeps
! the serial-only guard.
!===============================================================================
module spectrum_lobpcg_module

  use compilation, only : dp
  use geninfo,     only : stp, dv

  implicit none

  private
  public :: lobpcg_solve

  integer, parameter :: lobpcg_max_iter = 200   ! Iteration ceiling of the solve

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

  abstract interface
     !-------------------------------------------------------------------------
     ! Apply M^{-1} to one mesh spinor (the residual of state k); the column
     ! index lets the caller use per-state parameters.
     !-------------------------------------------------------------------------
     subroutine apply_precond_interface(mv, r, k, mr)
       import dp
       integer, intent(in)  :: mv, k
       real(dp), intent(in) :: r(mv,4)
       real(dp), intent(out):: mr(mv,4)
     end subroutine apply_precond_interface
  end interface

contains

  !-----------------------------------------------------------------------------
  ! Solves for the n lowest eigenpairs of h in the block by fixed-block LOBPCG.
  ! The mesh integration weight dv enters only in the projected matrix and
  ! the residual norms, so thetas and residuals carry energy units and the
  ! residual norm matches the locking criterion ||r||^2 < locking_treshold.
  !-----------------------------------------------------------------------------
  ! Input:
  !   apply_h   - applies h to one mesh spinor (see apply_h_interface)
  !   apply_m   - optional: applies M^{-1} to one residual column (see
  !               apply_precond_interface); absent means M = identity
  !   mv        - mesh-spinor dimension
  !   n         - subspace dimension: the block's spwf count
  !   dv        - mesh integration weight
  !   x         - initial guesses, (mv,4,n)
  !   tolerance - uniform residual criterion on sqrt(||r||^2)
  ! Output:
  !   x         - the lowest-n eigenvectors, orthonormal, ascending theta
  !   thetas    - the Ritz eigenvalues
  !   residuals - the residual norms sqrt(dv ||h x||^2 - theta^2)
  !   iterations - optional: the iteration count at convergence
  ! Errors out (stp) at the iteration ceiling: the caller's overallocation
  ! covers the count, so non-convergence means the solve failed - no silent
  ! truncation of the spectrum.
  !-----------------------------------------------------------------------------
  subroutine lobpcg_solve(apply_h, mv, n, dv, x, tolerance, thetas, residuals, &
       &                 apply_m, tolerances, iterations)

    procedure(apply_h_interface)      :: apply_h
    procedure(apply_precond_interface), optional :: apply_m
    integer, intent(in)          :: mv, n
    real(dp), intent(in)         :: dv
    real(dp), intent(inout)      :: x(mv,4,n)
    real(dp), intent(in)         :: tolerance
    real(dp), intent(out)        :: thetas(n), residuals(n)
    ! Per-state tolerances: states the caller does not need converged (the
    ! overallocated block tail above the cutoff) get a loose one, so their
    ! slow tail convergence cannot stall the uniform criterion.
    real(dp), intent(in), optional :: tolerances(n)
    integer, intent(out), optional :: iterations

    real(dp), allocatable :: y(:,:,:), hy(:,:,:), g(:,:), c(:,:), work(:)
    real(dp), allocatable :: hx_store(:,:,:), p_store(:,:,:), residual(:)
    real(dp), allocatable :: ynew(:,:)
    real(dp), allocatable :: w(:)
    integer               :: m, k, it, lwork, info

    external :: DSYEV

    if(mv <= 0) call stp('lobpcg_solve: mesh not initialized.',                  &
         &               'lobpcg_solve')
    if(n  <= 0) call stp('lobpcg_solve: empty block.', 'lobpcg_solve')

    ! Allocatable locals persist across calls (one call per symmetry block);
    ! reallocate for this block's mv and n.
    if(allocated(y)) deallocate(y, hy, hx_store, p_store, residual)
    allocate(y(mv,4,3*n), hy(mv,4,3*n), hx_store(mv,4,n),                    &
         &   p_store(mv,4,n), residual(n), ynew(mv,4))
    y = 0.0_dp

    do it=1, lobpcg_max_iter

       m = 2*n ; if(it > 1) m = 3*n

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ! The basis: {X, hX} cold, {X, R, P} afterwards. With a preconditioner
       ! the residual block carries M^{-1} R. A residual already below the
       ! tolerance is replaced by a fresh random spinor to keep the basis full
       ! rank. Orthonormalized before h acts on it.
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       if(it == 1) then
          y(:,:,1:n) = x
          do k=1,n
             y(:,:,n+k) = apply_h(mv, y(:,:,k))
          enddo
       else
          y(:,:,1:n)     = x
          y(:,:,n+1:2*n) = hx_store - state_times_theta(mv, n, x, thetas)
          y(:,:,m-n+1:3*n) = p_store
          do k=1,n
             if(residuals(k) < tolerance) then
                ! The state meets the criterion: its residual column is
                ! round-off noise and its conjugation direction is stale,
                ! both destabilizing the projected solve. Fresh random
                ! columns keep the basis full rank for the others.
                y(:,:,n+k)   = random_state(mv, k)
                y(:,:,m-n+k) = random_state(mv, 100+k)
             endif
          enddo
          if(present(apply_m)) then
             do k=1,n
                call apply_m(mv, y(:,:,n+k), k, ynew)
                y(:,:,n+k) = ynew
             enddo
          endif
       endif
       call orthonormalize_block(mv, m, dv, y)

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ! The action of h on the whole basis. The X block is free: x is a
       ! rotation of the previous orthonormal basis and hx_store of the
       ! previous hY by the SAME coefficient matrix, so h x = hx_store;
       ! only the R and P columns (m - n of them) need applications.
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       if(it == 1) then
          do k=1,m
             hy(:,:,k) = apply_h(mv, y(:,:,k))
          enddo
       else
          hy(:,:,1:n) = hx_store
          do k=n+1,m
             hy(:,:,k) = apply_h(mv, y(:,:,k))
          enddo
       endif

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ! The projected matrix G = Y^T hY and its lowest 2n Ritz pairs: the
       ! first n update X, the next n seed the next P.
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       if(allocated(g)) deallocate(g, c, w)
       allocate(g(m,m), c(m,n), w(m))
       call DGEMM('t','n', m, m, 4*mv, dv, y, 4*mv, hy, 4*mv,                   &
            &      0.0_dp, g, m)
       lwork = -1 ; allocate(work(1))
       call DSYEV('V', 'U', m, g, m, w, work, lwork, info)
       lwork = int(work(1)) ; deallocate(work) ; allocate(work(lwork))
       call DSYEV('V', 'U', m, g, m, w, work, lwork, info)
       deallocate(work)
       if(info /= 0) call stp('lobpcg_solve: projected eigensolve '//           &
            &                 'failed.', 'lobpcg_solve')

       c = g
       call rotate_basis(mv, m, n, y, c, x)
       thetas = w(1:n)
       ! h on the new X comes from the rotation of hY, so the residual
       ! identity ||hX||^2 - theta^2 is exact; it also seeds the next R.
       call rotate_basis(mv, m, n, hy, c, hx_store)
       call rotate_basis(mv, m, n, y, c(:,n+1:2*n), p_store)

       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       ! Residuals of the new X, and the convergence measure.
       ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       do k=1,n
          residuals(k) = sqrt(max(0.0_dp, dv*sum(hx_store(:,:,k)**2)            &
           &                                      - thetas(k)**2))
       enddo

       print '(a,i0,a,es10.2)', ' LOBPCG iteration ', it,                       &
            &                  ', max residual ', maxval(residuals)
       if(present(tolerances)) then
          if(all(residuals < tolerances)) then
             if(present(iterations)) iterations = it
             return
          endif
       elseif(all(residuals < tolerance)) then
          if(present(iterations)) iterations = it
          return
       endif

    enddo

    call stp('lobpcg_solve: the solve did not converge within the '//           &
         &   'iteration ceiling; the spwf set does not span the cutoff.',       &
         &   'lobpcg_solve')

  end subroutine lobpcg_solve

  !-----------------------------------------------------------------------------
  ! Rotates the basis y by the coefficient columns c: x = y c.
  !-----------------------------------------------------------------------------
  subroutine rotate_basis(mv, m, n, y, c, x)
    integer, intent(in)  :: mv, m, n
    real(dp), intent(in)  :: y(mv,4,m), c(m,n)
    real(dp), intent(out) :: x(mv,4,n)

    call DGEMM('n','n', 4*mv, n, m, 1.0_dp, y, 4*mv, c, m, 0.0_dp, x, 4*mv)

  end subroutine rotate_basis

  !-----------------------------------------------------------------------------
  ! X Theta as mesh-spinor columns.
  !-----------------------------------------------------------------------------
  function state_times_theta(mv, n, x, thetas) result(xt)
    integer, intent(in)  :: mv, n
    real(dp), intent(in)  :: x(mv,4,n), thetas(n)
    real(dp)              :: xt(mv,4,n)

    integer :: k

    do k=1,n
       xt(:,:,k) = thetas(k)*x(:,:,k)
    enddo

  end function state_times_theta

  !-----------------------------------------------------------------------------
  ! A fresh random mesh spinor, seeded per column after the house pattern:
  ! reproducible, golden_diff-safe.
  !-----------------------------------------------------------------------------
  function random_state(mv, k) result(v)
    integer, intent(in)  :: mv, k
    real(dp)             :: v(mv,4)

    integer               :: s
    integer, allocatable  :: seed(:)
    real(dp)              :: u(mv,4)

    call random_seed(size=s)
    allocate(seed(s)) ; seed = 961 + k
    call random_seed(put=seed) ; deallocate(seed)
    call random_number(u)
    v = u / sqrt(sum(u*u))

  end function random_state

  !-----------------------------------------------------------------------------
  ! Modified Gram-Schmidt orthonormalization of m mesh-spinor columns in the
  ! dv-weighted mesh metric: y_i^T dv y_j = delta_ij.
  !-----------------------------------------------------------------------------
  subroutine orthonormalize_block(mv, m, dv, y)
    integer, intent(in)     :: mv, m
    real(dp), intent(in)    :: dv
    real(dp), intent(inout) :: y(mv,4,m)

    real(dp) :: hn
    integer  :: k, j

    do k=1,m
       do j=1,k-1
          y(:,:,k) = y(:,:,k) - dv*sum(y(:,:,j) * y(:,:,k)) * y(:,:,j)
       enddo
       hn = sqrt(dv*sum(y(:,:,k)**2))
       if(hn < tiny(1.0_dp))                                                  &
            &   call stp('orthonormalize_block: degenerate basis column.',    &
            &                'orthonormalize_block')
       y(:,:,k) = y(:,:,k) / hn
    enddo

  end subroutine orthonormalize_block

end module spectrum_lobpcg_module
