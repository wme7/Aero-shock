program format1
  implicit none
  integer :: i

  do i=1,20
     write(*,2) i, i**2, i**3
  end do
2 format(i4,i6,i8)

end program format1
