program swapmain2

  use swapmod ! use statements must come first

  implicit none

  real :: a, b

  ! read in two values
  read(*,*) a, b
 
  call swap(a,b)
  write(*,*) a, b

end program swapmain2
