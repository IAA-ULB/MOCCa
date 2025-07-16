program circle_test
  use class_Circle
  implicit none

  type(Circle) :: c     ! Declare a variable of type Circle.
  type(Circle), pointer :: pc
  integer :: status

  c = Circle(1.5)       ! Use the implicit constructor, radius = 1.5.
  call circle_print(c)  ! Call a class subroutine

  allocate(pc, stat=status)
  write (*,*) 'allocate   status :',status
  pc%radius = 5.
  call circle_print(pc)
  deallocate(pc, stat=status)
  write (*,*) 'deallocate status :',status

end program circle_test
