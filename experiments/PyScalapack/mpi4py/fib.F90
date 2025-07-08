! fib.f90
subroutine fib(a, n)
  use iso_c_binding
   integer(c_int), intent(in) :: n
   integer(c_int), intent(out) :: a(n)
   do i = 1, n
      if (i .eq. 1) then
         a(i) = 0.0d0
      elseif (i .eq. 2) then
         a(i) = 1.0d0
      else
         a(i) = a(i - 1) + a(i - 2)
      end if
   end do
end