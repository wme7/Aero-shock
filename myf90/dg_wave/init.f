      subroutine init
c     ***************
c                
      include 'com'
c
      double precision pl(KK+1,KK+1)
      double precision ffunc(KK+1),phi(NINT,KK+1)
c  
c     computing c
c     -----------
      if (nproblem.eq.1.or.nproblem.eq.2) then
       n1=int((a-rl)/dx)
       do i=1,n1
        c(i)=1.d00
        c(nx+1-i)=1.d00
       enddo
       do i=n1+1,nx-n1
        c(i)=cc
       enddo
      endif
c
c     computing the cij-s for the interior
c     numerical fluxes
c     ------------------------------------
      do i=1,nx-1
       cplus=c(i+1)
       cmins=c(i)
c
       cave=(cplus+cmins)/2.d00
       ciave=(1.d00/cplus+1.d00/cmins)/2.d00
       cjmp=cplus-cmins
c
       c11(i)=1.d00/(2.d00*cave)
       c12(i)=cjmp/(4.d00*cave)
       c22(i)=1.d00/(2.d00*ciave)
      enddo
c
c      computing the alpha coefficients of the
c      Runge-Kutta methods of order kt
c      ---------------------------------------
      alpha(1,0)=1.d00
      do m=2,kt
       alpha(m,0)=1.d00
       do j=1,m-2
        alpha(m,j)=alpha(m-1,j-1)/dble(j)
        alpha(m,0)=alpha(m,0)-alpha(m,j)
       enddo
       alpha(m,m-1)=alpha(m-1,m-2)/dble(m)
       alpha(m,0)=alpha(m,0)-alpha(m,m-1)
      enddo
c
c      nproblem=1,2
c      Computing the initial data
c      -------------------------
      if (nproblem.eq.1.or.nproblem.eq.2) then
c
c      computing the Guassian points
       call gauleg(-1.d00,1.d00,x,w,k+1)  
c
c      evaluating the Legendre polynomials at
c      the quadrature points
        do i=1,k+1
         pl(i,1)=1.d00
         pl(i,2)=x(i)
         do j=2,k+1
          pl(i,j+1)=((2.d00*j-1.d00)*x(i)*pl(i,j)
     *                -(j-1.d00)*pl(i,j-1))/j
         enddo
        enddo
c
c      computing the moments of the initial data in `func'
c      and storing them in `phi'
        do i=1,nx
         xic=-a+(2.d00*i-1.d00)*dx/2.d00
c
c        evaluating the function `func' at the quadrature points
          do m=1,k+1
            ffunc(m)=func(xic+x(m)*dx/2.d00)
          enddo
c
          do j=0,k
           phi(i,j+1)=0.d00
           do m=1,k+1
            phi(i,j+1)= phi(i,j+1)+ffunc(m)*pl(m,j+1)*w(m)
           enddo
            phi(i,j+1)= phi(i,j+1)*(2.d00*j+1.d00)/2.d00
          end do
        end do
c
c      setting the initial degrees of freedom
        do i=1,nx
          do j=1,k+1
             u(i,j,0)= phi(i,j)
             v(i,j,0)=-phi(i,j)
          end do
        end do
      endif 
c
      return
      end
