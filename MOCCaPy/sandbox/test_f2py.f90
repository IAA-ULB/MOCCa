subroutine fill(a, a0)
  ! fill 1D array with a0
  ! Arguments
    real*8, intent(inout):: a(:)
    real*8, intent(in)   :: a0
  ! variables
    integer :: i,n
    n = size(a)
    do i=1,n
        a(i) = a0
    end do
end subroutine fill