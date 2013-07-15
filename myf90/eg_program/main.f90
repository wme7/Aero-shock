program matrix_operations
! Example code on how to deal with array operations in Fortran 95. I'm also creating some functions so that I can imitate Matlab stile of programming. ;D
! coded by Manuel Diaz, NTU, 2013.07.11
! manuel.ade'at'gmail.com

!PACKAGES
use dispmodule
use mathmodule!, only: ones, zeros, inverse, transposev

!Preamble
implicit none

integer, parameter :: n=3, m=3
real, dimension(n,m) :: A,B,C 	!Matrix Arrays
real, dimension(1,m) :: x,z,ty 	!Row vectors
real, dimension(n,1) :: y,w,tz 	!Column vectors
real, dimension(:,:), allocatable :: xx,zz !unknown dimension

! allocate variables in memory
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
call disp('A*B*T(A)=',A*B*transpose(A))
call disp('x*T(y) = ',x*transpose(y))
call disp(' A*A*B = ',A*A*B)
C = sum(B*A*transpose(B))
call disp('   C   = ',C)

! matrix multiplication and inner product
call disp('matmul(tA,B)= ',matmul(transpose(A),B))
call disp('matmul(A,B) = ',matmul(A,B))
call disp('matmul(x,y) = ',matmul(x,y))
call disp('matmul(A,x) = ',matmul(A,y))
call disp('matmul(ty,B)= ',matmul(transpose(y),B))

! Sum of elements in arrays
call disp('sum(A,1) = ',sum(A,dim=1))
z = reshape(sum(A,dim=1),shape(z))
call disp('sum(A,1) = ',z)
!call disp('sum(A,1) = ',reshape(sum(A,dim=1),shape(z)) <-can't
y = reshape(sum(A,dim=1),shape(y))
call transposev(y,ty)
call disp('tsum(A,1)= ',ty)
call disp('sum(A,2) = ',sum(A,dim=2))

! Create zeros array
call disp('ones(3,3) = ',ones(3,3))
call disp('ones(1,3) = ',ones(1,3))

! Create ones array
call disp('zeros(4,4)= ',zeros(4,4))
call disp('zeros(4,1)= ',zeros(4,1))

! Invert array
call disp('inv(B^2) = ',inverse(matmul(B,B)))

! Print array using dispmodule.f90
!call disp('C = ',C)
!call disp('z = ',z)

end program matrix_operations
