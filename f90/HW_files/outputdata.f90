program outputdata
  implicit none
  real, dimension(100) :: x, y
  integer :: i

  ! setup x and y with some data 
  do i=1,100
     x(i) = i*0.1
     y(i) = sin(x(i))*(1-cos(x(i)/3.0))
  end do

  ! output data to file
  open(1, file='data1.dat', status='new')
  do i = 1,100
     write(1,*) x(i), y(i)
  end do
  close(1)

end program outputdata
