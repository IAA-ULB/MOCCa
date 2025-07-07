program test_abs
  integer :: i = 1
  integer :: j = 2
  real :: x = -1.e0
  
  complex :: z = (-1.e0,0.e0)
  i = sqrt(float(abs(i-j)))
  x = abs(x)
  x = abs(z)
end program test_abs