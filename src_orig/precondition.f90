!===============================================================================
!  #######   ##   #    # #####   ##   #      #    #  ####
!     #     #  #  ##   #   #    #  #  #      #    # #
!     #    #    # # #  #   #   #    # #      #    #  ####
!     #    ###### #  # #   #   ###### #      #    #      #
!     #    #    # #   ##   #   #    # #      #    # #    #
!     #    #    # #    #   #   #    # ######  ####   ####
!
!  Copyright W. Ryssens & M. Bender
!
!===============================================================================
!
! Module that governs all of the possible preconditioning that can be applied.
!
!
! Two different possibilities are currently presented. 
! 
!===============================================================================
! (1) P.G. Style preconditioning for the spwfs
!===============================================================================
!     Constructing the approximate inverse matrices of the second derivatives in
!     every Cartesian direction with an added constant.
!     I.e. that
!      A_x{x/y/z} (C * Delta_{x/y/z} - eps)^{-1}
!     Then we suppose that 
!      (Delta + D)^{-1} \approx A_{x} + A_{y} + A_{z}
!     for some specific values of C,D and eps.
!  
!===============================================================================
! (2) W.R.-style preconditioning for the potentials
!===============================================================================
! 
!     Invert the following operator on the mesh
!         Pf(r) = b*f(r) + a*Delta[f(r)]
!
!     Not by constructing its inverse, but by repeatedly applying that 
!     operator in a conjugate gradient scheme.
!===============================================================================

module preconditioning

    use derivatives

    implicit none
    
contains

 !------------------------------------------------------------------------------
 !   (1)
 !------------------------------------------------------------------------------
 subroutine InvertDerivatives(eps,C,invX,invY,invZ)
    !---------------------------------------------------------------------------
    ! Construct the inverse matrices of the second Lagrange derivative matrices.
    ! With some parameters ( epsilon and c) to make the thing look like a true
    ! inverted Hamiltonian.
    !
    ! Calculates and returns matrices invX, invY and invZ that are
    !    
    !      (C * Delta_{x/y/z} - eps)^{-1}
    !---------------------------------------------------------------------------
    ! OUTPUT
    real(KIND=dp),intent(out) :: invX(nx,nx,2),invY(ny,ny,2),invZ(nz,nz,2)
    ! INPUT 
    real(KIND=dp),intent(in)  :: eps, C 
    
    integer :: pivotx(nx)
    integer :: pivoty(ny)
    integer :: pivotz(nz)    
    integer :: ierror, pm,i
    real(KIND=dp) :: work(nz)

    !---------------------------------------------------------------------------
    ! Invert the shifted Laplacians
    do pm=1,2
        invX(:,:,pm) = C*laplaX(:,:,pm)
        do i=1,nx
            invX(i,i,pm) = invX(i,i,pm) - eps
        enddo
        call dgetrf (nx, nx, invX(:,:,pm), nx,pivotx, ierror)
        call dgetri (nx, invX(:,:,pm), nx, pivotx, work, nx, ierror) 
    enddo

    do pm=1,2
        invY(:,:,pm) = C*laplaY(:,:,pm)
        do i=1,ny
            invY(i,i,pm) = invY(i,i,pm) - eps
        enddo
        call dgetrf (ny, ny, invY(:,:,pm), ny,pivoty, ierror)
        call dgetri (ny, invY(:,:,pm), nx, pivoty, work, ny, ierror) 
    enddo

    do pm=1,2
        invZ(:,:,pm) = C*laplaZ(:,:,pm)
        do i=1,nz
            invZ(i,i,pm) = invZ(i,i,pm) - eps
        enddo
        call dgetrf (nz, nz, invZ(:,:,pm), nz,pivotz, ierror)
        call dgetri (nz, invZ(:,:,pm), nz, pivotz, work, nz, ierror) 
    enddo
    !---------------------------------------------------------------------------
    if(ierror.ne.0) then
        print *, 'Error in inverting the shifted laplacians.'
        stop
    endif
 end subroutine InvertDerivatives
 !------------------------------------------------------------------------------
 !   (2)
 !------------------------------------------------------------------------------

 function PreconditionPotential(pot,a,b,sx,sy,sz) result(invpot)
    !---------------------------------------------------------------------------
    ! Precondition a potential with the matrix
    !
    !   ( B  + A * \Delta)^{-1}
    !
    ! Calculated through the repeated application of its inverse in a 
    ! conjugate gradient algorithm. 
    !---------------------------------------------------------------------------
    use Derivatives
    
    real(KIND=dp), intent(in) :: a, b
    real(KIND=dp), intent(in) :: pot(nx*ny*nz,2)
    integer, intent(in) :: sx,sy,sz

    real(KIND=dp) :: residual(nx*ny*nz), update(nx*ny*nz), aCG
    real(KIND=dp) :: invpot(nx*ny*nz,2), direction(nx*ny*nz), bCG
    real(KIND=dp) :: newresnorm, oldresnorm
    
    integer:: it, iter
    
    !---------------------------------------------------------------------------
    invpot  = 0.0
    do it=1,2 
        Residual         = pot(:,it) 
        Direction        = Residual
        newresnorm       = sum(direction**2)*dv
        !-----------------------------------------------------------------------
        do iter=1,300
          update   = preconoperator(direction,a,b,sx,sy,sz)

          aCG    = NewResNorm/(sum(Direction*update)*dv)
          
          invpot(:,it)   = invpot(:,it)  + aCG * Direction
          residual       = residual      - aCG * update
          
          oldresnorm = newresnorm
          newresnorm = sum(residual**2)*dv
          
          BCG        = NewResNorm/OldResNorm
          Direction  = Residual + bCG * Direction
          if(newresnorm.lt.1d-8) exit
        enddo
        !-----------------------------------------------------------------------
    enddo
    !---------------------------------------------------------------------------
!    print *, 'Inv. Pot., iter = ', iter, newresnorm, sum(invpot(:,1))*dv,  &
!    &                     sum(invpot(:,2))*dv
    
 end function Preconditionpotential


 function preconoperator(f,a,b,sx,sy,sz) result(Pf)
    !---------------------------------------------------------------------------
    ! Implement the preconditioning operator
    !
    ! Pf(r) = B(r)*f(r) + a(r) * Delta(f(r)) 
    !---------------------------------------------------------------------------
    use Derivatives
    
    real(KIND=dp), intent(in) :: f(nx*ny*nz)
    real(KIND=dp)             :: Pf(nx*ny*nz), df(nx*ny*nz)
    real(KIND=dp),intent(in)  :: a, b
    integer, intent(in)       :: sx,sy,sz

    call Derive_lap(f,sx,sy,sz,df) 
    Pf = b*f + a*df
    
  end function preconoperator

end module preconditioning
