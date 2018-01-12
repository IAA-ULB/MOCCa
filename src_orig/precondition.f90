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
!===============================================================================

module preconditioning

    use derivatives

    implicit none
    
contains

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
    integer :: ierror, pm,i, k
    real(KIND=dp), allocatable:: toinvert(:,:)
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

end module preconditioning
