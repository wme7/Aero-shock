program fileio
  implicit none
  integer :: i

  open(20, file='cubes.dat')
  open(30, file='cubes_copy.dat')
  do i=1,100
     write(20,1) i, i**2, i**3
     write(30,1) i, i**2, i**3
  end do
  close(20)
  close(30)
1 format(i4,i6,i8)

end program fileio
