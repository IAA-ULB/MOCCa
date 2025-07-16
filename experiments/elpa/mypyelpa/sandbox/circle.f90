module class_Circle
  use iso_c_binding

  implicit none
  private
  public :: Circle, circle_area, circle_print

  real :: pi = 3.1415926535897931d0 ! Class-wide private constant

  type Circle
     real :: radius 
  end type Circle
contains
  function circle_area(this) result(area)
    type(Circle), intent(in) :: this
    real :: area
    area = pi * this%radius**2
  end function circle_area

  subroutine circle_print(this)

    type(Circle), intent(in) :: this
    real :: area

    area = circle_area(this)  ! Call the circle_area function
    print *, 'Circle: r = ', this%radius, ' area = ', area
  end subroutine circle_print

  function circle_allocate_c(error) result(ptr) bind(C, name="alloc_")

    integer(kind=c_int)        :: error
    integer                    :: error2
    type(c_ptr)                :: ptr
    type(circle), pointer      :: obj

    ! obj => elpa_impl_allocate(error)
    allocate(obj, stat=error2)
    if (error2 .ne. 0) then
      write(*, *) "alloc(): could not allocate object"
    endif
    error = error2
    ptr = c_loc(obj)
  end function

  subroutine circle_deallocate_c(handle, error) bind(C, name="free_")

    type(c_ptr), value         :: handle
    type(circle), pointer :: self
    integer(kind=c_int)        :: error

    call c_f_pointer(handle, self)
    ! call self%destroy(error)
    deallocate(self)
  end subroutine
  
subroutine set_radius(handle,radius) bind(C, name="set_radius_")

    type(c_ptr), value    :: handle
    type(circle), pointer :: self
    real, intent(in)      :: radius

    call c_f_pointer(handle, self)
    ! call self%destroy(error)
    self%radius = radius
    write(*,*) self%radius
  
  end subroutine set_radius

  subroutine print(handle) bind(C, name="print_")
  
    type(c_ptr), value    :: handle
    type(circle), pointer :: self

    call c_f_pointer(handle, self)
    call circle_print(self)

  end subroutine print

  subroutine test
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

  end subroutine test
  
!   subroutine test2
!     type(c_ptr) :: pc
!     write(*,*) 'test2'
!     pc => allocate_circle()
!     pc%radius = 6
!     call circle_print(pc)
!     call deallocate_circle(pc)
!   end subroutine test2
  
  
end module class_Circle
