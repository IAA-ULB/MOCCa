program tantalus_multiple

  use Tantalus

  implicit none
	
  integer*8 :: a=1, b =3
  
  call Run_Tantalus(a, 'data.one.in') 
  call Run_Tantalus(b, 'data.two.in') 

end program 
