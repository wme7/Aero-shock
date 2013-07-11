program matrix_operations
! Example code on how to deal with array operations in Fortran 95.
! coded by Manuel Diaz, NTU, 2013.07.11
! manuel.ade'at'gmail.com

!PACKAGES
use dispmodule
use mathmodule, only: ones, zeros, inverse

!Preamble
implicit none

integer, parameter :: n=3, m=3
real, dimension(n,m) :: A,B,C 	!Matrix Arrays
real, dimension(1,m) :: x,z 	!Row vectors
real, dimension(n,1) :: y 	!Column vectors
real, dimension(:,:), allocatable :: w,xx,zz
! allocate variables in memory
allocate (w(n,1))
allocate (zz(n,1))
allocate (xx(n,1))

! Build arrays
A = transpose( reshape( (/1,2,3,4,5,6,7,8,9/),shape(A)))
B = transpose( reshape( (/2,0,0,0,2,0,0,0,2/),shape(B)))
call disp('A = ',A)
call disp('B = ',B)

! Build vector arrays
x = reshape( (/1,2,3/),shape(x)) 
y = reshape( (/2,2,2/),shape(y))
w = reshape(x,(/size(x,2),size(x,1)/) )
zz= reshape(x,(/size(x,1),size(x,2)/) )
xx= transpose(x)
call disp('x = ',x)
call disp('y = ',y)
call disp('w = ',w)
call disp('zz= ',zz)
call disp('xx= ',xx)

! Vector Addition
call disp(' A + B = ',A+B)
call disp('x+T(y) = ',x+transpose(y))

! Dot product
call disp(' A * B = ',A*B)
call disp('x*T(y) = ',x*transpose(y))
call disp(' A*A*B = ',A*A*B)

! matrix multiplication and inner product
call disp('matmul(A,B) = ',matmul(A,B))
call disp('matmul(x,y) = ',matmul(x,y))

! Create zeros array
call disp('ones(3) = ',ones(3,3))

! Create ones array
call disp('zeros(4)= ',zeros(4,4))

! Invert array
call disp('invB = ',inverse(B))

! Print array using dispmodule.f90
!call disp('A = ',C)
!call disp('c = ',z)

end program matrix_operations
