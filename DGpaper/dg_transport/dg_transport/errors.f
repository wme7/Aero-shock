      subroutine errors
c     *****************
c
      include'com'
c
      double precision pl(KK+1,KK+1)
      double precision exact(KK+1)
      double precision apprx(KK+1)
      double precision errrr(KK+1)
c
c      --- computing the Gaussian points 
       kgauss=k+4
       call gauleg(-1.d00,1.d00,x,w,kgauss)  
c
c      evaluating the Legendre polynomials at
c      the quadrature points
        do i=1,kgauss
         pl(i,1)=1.d00
         pl(i,2)=x(i)
         do j=2,k
          pl(i,j+1)=((2.d00*j-1.d00)*x(i)*pl(i,j)
     *                -(j-1.d00)*pl(i,j-1))/j
         enddo
        enddo
c
c      --- computation of the L-infinity and L1-errors
c
c      --- initialization    
      errl1(ne)=0.d00
      errli(ne)=0.d00
c
c      --- loop on the elements
      do  i=1, nx
c
c        the midpoint of the element
        xic=(2.d00*i-1.d00)*dx/2.d00
c 
c        evaluating the errors  at the quadrature points
        do m=1,kgauss
          exact(m)=func(xic+x(m)*dx/2.d00,t)
          apprx(m)=0.d00
          do j=1,k+1
            apprx(m)=apprx(m)+u(i,j,0)*pl(m,j)
          enddo
          errrr(m)=abs(exact(m)-apprx(m))
        enddo
c
c        --- the errors
c        L1:
        aux=0.d00
        do m=1,kgauss
          aux=aux+errrr(m)*w(m)
        enddo
c        errl1(ne)=errl1(ne)+aux
        errl1(ne)=errl1(ne)+aux**2
c
c       L-infinity:
        do  m=1,kgauss
          if (errrr(m).gt.errli(ne)) errli(ne)=errrr(m)
        enddo
c
      enddo
c
      errl1(ne)=sqrt(ddx(ne)*errl1(ne))
c
      end


