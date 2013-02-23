program arraytest
implicit none
integer :: i,j,k
integer, parameter :: idim = 10
real, dimension(idim) :: x, y, z

x = (/ 1,2,3,4,5,6,7,8,9,0 /)
y = (/ 0,9,8,7,6,5,4,3,2,1 /)

z = x + y

print*, 'result =', x

end program
