program LUD
parameter(np = 3,n = 3,nmax = 100)
dimension a(np,np), b(np), indx(n), vv(nmax)

open (unit=8,file='a.tec',status='unknown')
open (unit=9,file='x.tec',status='unknown')

!do i = 1, np
!	do j = 1,np
!		a(i,j) = 2.0
!	end do
!end do

    A(1,1) = 2
	A(1,2) = 0
	A(1,3) = 0
	A(2,1) = 2
	A(2,2) = 2
	A(2,3) = 0
	A(3,1) = 2
	A(3,2) = 2
	A(3,3) = 2

!do i = 1,np
!	b(i) =1
!end do

    b(1)=3
	b(2)=2
	b(3)=1

do i = 1,n
	do j = 1,n
		write(8,*) a(i,j),indx,d
	end do
end do

call LUDCMP (a,n,np,indx,d)
call LUBKSB (a,n,np,indx,b)

do i = 1,n
	write(9,*)b(i)
end do
end program

subroutine LUDCMP (a,n,np,indx,d)
parameter(nmax = 100, tiny = 1.0E-20)
dimension a(np,np), indx(n), vv(nmax)
d = 1
do 12 I = 1, n
	aamax = 0
	do 11 j = 1, n
		if (abs(a(i,j)) .gt. aamax) aamax = abs(a(i,j))
	11 continue
	if (aamax .eq. 0.) pause 'singular matrix.' 
	vv(i) = 1./aamax
12 continue

do 19 j = 1, n
	do 14 i = 1, j-1
		sum = a(i,j)
		do 13 k = 1,i-1
			sum = sum - a(i,k) * a(k,j)
		13 continue
		a(i,j) = sum
	14 continue
	aamax = 0

	do 16 i = j, n
		sum = a(i,j)
		do 15 k = 1,j-1
			sum = sum - a(i,k) * a(k,j)
		15 continue
		a(i,j) = sum
		dum = vv(i)*abs(sum)
		if (dum.ge.aamax) then
			imax = i
			aamax = dum
		end if
	16 continue

	if (j.ne.imax) then
		do 17 k = 1,n
			dum = a(imax,k)
			a(imax,k) = a(j,k)
			a(j,k) = dum
		17 continue
		d = -d
		vv(imax) = vv(j)
	end if

	indx(j) = imax
	if (a(j,j).eq.0.) a(j,j) = tiny
	if (j.ne.n) then
		dum = 1./a(j,j)
		do 18 i = j+1,n
			a(i,j) = a(i,j) * dum
		18 continue
	end if
19 continue
return
end

subroutine LUBKSB (a,n,np,indx,b)
dimension a(np,np), indx(n), b(n)
ii = 0
do 12 i = 1,n
	ll = indx(i)
	sum = b(ll)
	b(ll) = b(i)
	if (ii.ne.0) then
		do 11 j = ii, i - 1
			sum = sum - a(i,j) * b(j)
		11 continue
	else if (sum .ne. 0) then
	ii = i
	end if
	b(i) = sum
12 continue

do 14 i = n, 1, -1
	sum = b(i)
	if (i.lt.n)	then
		do 13 j = i + 1,n
			sum = sum - a(i,j) * b(j)
		13 continue
	end if
	b(i) = sum/a(i,i)
14 continue
return
end