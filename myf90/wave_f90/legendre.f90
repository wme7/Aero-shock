      SUBROUTINE legendre(pl,x,ks,kn)
      INTEGER,intent(in):: kn,ks
      double precision, intent(in)::x(ks)
      double precision, intent(out)::pl(ks,kn)
      INTEGER i,j
!      evaluating the Legendre polynomials at
!      the quadrature points
        do i=1,ks
         pl(i,1)=1.d00
         pl(i,2)=x(i)
         do j=2,kn-1
          pl(i,j+1)=((2.d00*j-1.d00)*x(i)*pl(i,j)&
                     -(j-1.d00)*pl(i,j-1))/j
         enddo
        enddo
!      write(6,*) pl
      return
      end subroutine legendre