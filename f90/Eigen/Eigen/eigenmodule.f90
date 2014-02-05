module eigenmodule
implicit none
!
! Define Global Data here
!
contains
    SUBROUTINE tqli(d,e,z)
    USE nrtype; USE nrutil, ONLY : assert_eq, nrerror
    USE nr, ONLY : pythag
    IMPLICIT NONE
    REAL(SP), DIMENSION(:), INTENT(INOUT) :: d,e
    REAL(SP), DIMENSION(:,:), OPTIONAL, INTENT(INOUT) :: z
    INTEGER(I4B) :: i,iter,l,m,n,ndum
	REAL(SP) :: b,c,dd,f,g,p,r,s
	REAL(SP), DIMENSION(size(e)) :: ff
	n=assert_eq(size(d),size(e),'tqli: n')
	if (present(z)) ndum=assert_eq(n,size(z,1),size(z,2),'tqli: ndum')
	e(:)=eoshift(e(:),1)
	do l=1,n
		iter=0
		iterate: do
			do m=l,n-1
				dd=abs(d(m))+abs(d(m+1))
				if (abs(e(m))+dd == dd) exit
			end do
			if (m == l) exit iterate
			if (iter == 30) call nrerror('too many iterations in tqli')
			iter=iter+1
			g=(d(l+1)-d(l))/(2.0_sp*e(l))
			r=pythag(g,1.0_sp)
			g=d(m)-d(l)+e(l)/(g+sign(r,g))
			s=1.0
			c=1.0
			p=0.0
			do i=m-1,l,-1
				f=s*e(i)
				b=c*e(i)
				r=pythag(f,g)
				e(i+1)=r
				if (r == 0.0) then
					d(i+1)=d(i+1)-p
					e(m)=0.0
					cycle iterate
				end if
				s=f/r
				c=g/r
				g=d(i+1)-p
				r=(d(i)-g)*s+2.0_sp*c*b
				p=s*r
				d(i+1)=g+p
				g=c*r-b
				if (present(z)) then
					ff(1:n)=z(1:n,i+1)
					z(1:n,i+1)=s*z(1:n,i)+c*ff(1:n)
					z(1:n,i)=c*z(1:n,i)-s*ff(1:n)
				end if
			end do
			d(l)=d(l)-p
			e(l)=g
			e(m)=0.0
		end do iterate
	end do
	END SUBROUTINE tqli

	SUBROUTINE tred2(a,d,e,novectors)
	USE nrtype; USE nrutil, ONLY : assert_eq,outerprod
	IMPLICIT NONE
	REAL(SP), DIMENSION(:,:), INTENT(INOUT) :: a
	REAL(SP), DIMENSION(:), INTENT(OUT) :: d,e
	LOGICAL(LGT), OPTIONAL, INTENT(IN) :: novectors
	INTEGER(I4B) :: i,j,l,n
	REAL(SP) :: f,g,h,hh,scale
	REAL(SP), DIMENSION(size(a,1)) :: gg
	LOGICAL(LGT) :: yesvec
	n=assert_eq(size(a,1),size(a,2),size(d),size(e),'tred2')
	if (present(novectors)) then
		yesvec=.not. novectors
	else
		yesvec=.true.
	end if
	do i=n,2,-1
		l=i-1
		h=0.0
		if (l > 1) then
			scale=sum(abs(a(i,1:l)))
			if (scale == 0.0) then
				e(i)=a(i,l)
			else
				a(i,1:l)=a(i,1:l)/scale
				h=sum(a(i,1:l)**2)
				f=a(i,l)
				g=-sign(sqrt(h),f)
				e(i)=scale*g
				h=h-f*g
				a(i,l)=f-g
				if (yesvec) a(1:l,i)=a(i,1:l)/h
				do j=1,l
					e(j)=(dot_product(a(j,1:j),a(i,1:j)) &
					+dot_product(a(j+1:l,j),a(i,j+1:l)))/h
				end do
				f=dot_product(e(1:l),a(i,1:l))
				hh=f/(h+h)
				e(1:l)=e(1:l)-hh*a(i,1:l)
				do j=1,l
					a(j,1:j)=a(j,1:j)-a(i,j)*e(1:j)-e(j)*a(i,1:j)
				end do
			end if
		else
			e(i)=a(i,l)
		end if
		d(i)=h
	end do
	if (yesvec) d(1)=0.0
	e(1)=0.0
	do i=1,n
		if (yesvec) then
			l=i-1
			if (d(i) /= 0.0) then
				gg(1:l)=matmul(a(i,1:l),a(1:l,1:l))
				a(1:l,1:l)=a(1:l,1:l)-outerprod(a(1:l,i),gg(1:l))
			end if
			d(i)=a(i,i)
			a(i,i)=1.0
			a(i,1:l)=0.0
			a(1:l,i)=0.0
		else
			d(i)=a(i,i)
		end if
	end do
	END SUBROUTINE tred2

	SUBROUTINE sort2(arr,slave)
	USE nrtype; USE nrutil, ONLY : assert_eq
	USE nr, ONLY : indexx
	IMPLICIT NONE
	REAL(SP), DIMENSION(:), INTENT(INOUT) :: arr,slave
	INTEGER(I4B) :: ndum
	INTEGER(I4B), DIMENSION(size(arr)) :: idx
	ndum=assert_eq(size(arr),size(slave),'sort2')
	call indexx(arr,idx)
	arr=arr(idx)
	slave=slave(idx)
	END SUBROUTINE sort2

end module eigenmodule
