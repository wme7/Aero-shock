Program Matrix_Inversion
implicit none
! this routine is designed for matrix inversion
! matrix inversion = LU decomps. + backsubstitution 
integer n,m, nmax
parameter (n=10,m=10,nmax=500)
real(8) a(n),b(n),c(n),r(n,m),u(n,m),tu(n,m),t(n,m)
!	This program inverses matrix [A]
!	[A] is a tridiagonal matrix with [a,b,c]
!	[A]^{-1} is stored in u(i,j)
!	r(i,j) is an identity matrix [I]
integer i,j,k
real bet,gam(nmax)
open (unit = 10, file = 'u.tec', status = 'unknown')
open (unit = 20, file = 't.tec', status = 'unknown')
open (unit = 30, file = 'tu.tec', status = 'unknown')
do i = 2,n
	a(i) = 0
end do
do i = 1,n
	b(i) = 2
end do
do	i = 1,n-1
	c(i) = 0
end do
do i = 1,n
	do j = 1,m
		r(i,j) = 0
	end do
	r(i,i) = 1
end do
do i = 1,n
	do j = 1,m
		t(i,j) = 0
	end do
	t(i,i) = 2
end do

do  k = 1,m
	bet=b(1)
	u(1,k) = r(1,k)/bet
	do  j = 2,n
		gam(j) = c(j-1)/bet
		bet=b(j)-a(j)*gam(j)
		u(j,k) = (r(j,k)-a(j)*u(j-1,k))/bet
	end do 
	do  j = n-1,1,-1
		u(j,k) = u(j,k) - gam(j+1)*u(j+1,k)
	end do 
end do
!
do i = 1,n
  do j = 1,n
    tu(i,j) = 0.
  end do
end do

do i = 1,n
  do j = 1,n
    do k = 1,n
	  tu(i,j)=tu(i,j) + t(i,k)*u(k,j)
    end do
  end do
end do
do i = 1,n
do j = 1,m
   write (10,*) i,j,u(i,j)
   write (20,*) i,j,t(i,j)
   write (30,*) i,j,tu(i,j)
end do 
end do 
end program 