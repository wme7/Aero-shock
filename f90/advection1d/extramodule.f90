module extramodule
implicit none
 !
 ! here goes any global variables
 !
contains

 function linspace(a,b,n)
 integer :: n,i
 real :: a,b,dx
 real :: linspace(n)

 dx = (b-a)/(n-1)
 do i = 0,n-1
  linspace(i+1) = a + i*dx
 end do 

 end function linspace

end module extramodule
