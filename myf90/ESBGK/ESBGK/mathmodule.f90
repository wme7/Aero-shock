module mathmodule
implicit none
    !
    ! Global declarations would appear here it any
    !
contains
! routines provided by this module:

function ones(n,m)
implicit none
integer :: n,m
real, dimension(n,m) :: ones
! Using fortran 95:
ones = 1.0
end function ones

function ones1d(n)
implicit none
integer :: n
real, dimension(n) :: ones1d
! Using fortran 95:
ones1d = 1.0
end function ones1d

function zeros(n,m)
implicit none
integer :: n,m
real, dimension(n,m) :: zeros
! Using fortran 95:
zeros = 0.0
end function zeros

function linspace(a,b,n)
implicit none
integer :: n,i
real :: a,b,dx
real linspace(n)
! using fortran 95
dx = (b-a)/(n-1)
linspace(1) = a
do i = 1,n
    linspace(i) = a + (i-1)*dx
end do
end function linspace

function kron(a,b)
implicit none
real :: a(:),b(:)
integer :: i,j,Na,Nb
real :: kron(size(a)*size(b))
Na = ubound(a,1)
Nb = ubound(b,1)
do i=1,Na
kron((i-1)*Nb+1:i*Nb) = a(i)*b
end do
end function kron

subroutine transposev(a,ta)
! transpose and vector array
! becuase fortran can olny transpose rank 2 arrays
implicit none
integer :: n,m,i,j
real, dimension(:,:), intent(in)  :: a
real, dimension(:,:), intent(out) :: ta
n = size(a,dim=1)
m = size(a,dim=2)
do i=1,n
  do j =1,m
  ta(j,i) = a(i,j)
  end do
end do
end subroutine transposev

function inverse(a)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed
! during the calculation
!===========================================================
implicit none
! External variables
real, dimension(:,:) :: a

! Internal variables
integer i, j, k, n, m
real, allocatable, dimension(:,:) :: inverse, L, U, c
real, allocatable, dimension(:) :: b, d, x
real coeff

! Allocate variables
n = size(a,dim=1)
m = size(a,dim=2)
allocate(inverse(n,m))
allocate(L(n,m))
allocate(U(n,m))
allocate(c(n,m))
allocate(b(n))
allocate(d(n))
allocate(x(n))

! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
L=0.0
U=0.0
b=0.0
! step 1: forward elimination
do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
end do

! Step 2: prepare L and U matrices
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
do i=1,n
  L(i,i) = 1.0
end do
! U matrix is the upper triangular part of A
do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
end do

! Step 3: compute columns of the inverse matrix C
do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
end do
inverse = c
end function inverse

end module mathmodule
