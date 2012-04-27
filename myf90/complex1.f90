program complex1
  implicit none

  !variables & constants
  complex, parameter :: i = (0,1) ! sqrt(-1)
  complex :: x, y

  x = (1,1); y = (1,-1)
  write(*,*) i*x*y

end program complex1
