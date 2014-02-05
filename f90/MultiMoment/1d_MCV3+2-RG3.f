C     A Fortran program to test three 3-point multi-moment constrained 
C     collocation schemes, i.e. MCV3, MCV3_UPCC and MCV3_CPCC schemes
C     Feng Xiao July 12, 2012

      PROGRAM ONE_D
      PARAMETER (nx=100,ns=3)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 u(0:nx+1,ns),v(0:nx+1,ns),f_x(0:nx+1,ns),fm_x(0:nx+1,ns)
      real*8 ue(ns),xi(ns),dx
      real*8 x(0:nx+1,ns)
      real*4 ua0(nx*ns),ua(nx*ns),xa(nx*ns)
      real*8 f1(0:nx+1,ns),f2(0:nx+1,ns),f3(0:nx+1,ns),f4(0:nx+1,ns)
      real*8 u0(0:nx+1,ns)
      real*4 uout(nx),xout(nx)
      integer ischeme
C
      pi=4.0*atan(1.0d0)

      WRITE(*,*) 'Choose a scheme MCV3(1),UPCC(2),CPCC(3)'
      read(*,*) ischeme
      do k=1,ns
      if(ischeme.eq.1.or.ischeme.eq.2) then
        xi(k)=2.*dble(k-1)/dble(ns-1)-1.0
        else
        xi(k)=-dcos((dble(k)-0.5)/dble(ns)*Pi)
       endif
      enddo
      write(*,*) xi
      dx=1.0/dble(nx)
      dt=0.2*dx
      nmax=nint(1.0/dt)

      do i=0,nx
         do k=1,ns
          x(i,k)=dble(i)*dx+0.5*xi(k)*dx
         enddo
      enddo


      do i=1,nx
         do k=1,ns
            v(i,k)=-1.
            u(i,k)=dexp(-400.0*(x(i,k)-0.5)**2)
cx            u(i,k)=0.
cx            if(x(i,k).ge.0.4. and. x(i,k).le.0.6)  u(i,k)=1.0
         enddo
      enddo

      call bdc(u,nx,ns)
      call bdc(v,nx,ns)



      do nn=1,1000

      do i=0,nx+1
         do k=1,ns
            u0(i,k)=u(i,k)
         enddo
      enddo

      call  mflux_x(u,v,fm_x,xi,dx,nx,ns,ischeme)

      do i=0,nx+1
         do k=1,ns
            f1(i,k)=fm_x(i,k)
         enddo
      enddo

      do i=1,nx
         do k=1,ns
            u(i,k)= u0(i,k)-dt*fm_x(i,k) 
         enddo
      enddo
      call bdc(u,nx,ns)

      call  mflux_x(u,v,fm_x,xi,dx,nx,ns,ischeme)

      do i=0,nx+1
         do k=1,ns
            f2(i,k)=fm_x(i,k)
         enddo
      enddo

      do i=1,nx
         do k=1,ns
            u(i,k)= u0(i,k)-dt*(f1(i,k)+f2(i,k))/4.0 
         enddo
      enddo
      call bdc(u,nx,ns)

      call  mflux_x(u,v,fm_x,xi,dx,nx,ns,ischeme)

      do i=0,nx+1
         do k=1,ns
            f3(i,k)=fm_x(i,k)
         enddo
      enddo

      do i=1,nx
         do k=1,ns
            u(i,k)= u0(i,k)
     6              -dt/6.0*(f1(i,k)+f2(i,k)+4.*f3(i,k)) 
         enddo
      enddo
      call bdc(u,nx,ns)

cx ********** output results *******************************
      do i=1,nx
      xout(i)=dble(i)*dx-0.5*dx
      uout(i)=(1./3.*u(i,1)+4./3.*u(i,2)+1./3.*u(i,3))
cx      uout(i)=(4./9.*u(i,1)+10./9.*u(i,2)+4./9.*u(i,3))
      enddo
	open(21,file='mcv3.dat')
        do i=1,nx
        write(21,1012) 2.*xout(i)-1., uout(i)/2. 
        enddo
        close(21)

 1012	format(f10.8,1x,F20.16)
cx  ********************************************************

      write(*,*) 'N=',nn  


        enddo

cx      CALL XCLOSE(1)
C
C===========================
C

      STOP
      END


      subroutine mflux_x(u,v,fm_x,xi,dx,nx,ns,ischeme)
C     MCV3 & MCV3_UPCC 
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 u(0:nx+1,ns),v(0:nx+1,ns),fm_x(0:nx+1,ns)
      real*8 f(0:nx+1,ns),fb(0:nx+1),fb_x(0:nx+1)
      real*8 xi(ns),dx
      integer nx,ns
 
      do i=0,nx+1
       do k=1,ns
       f(i,k)=v(i,k)*u(i,k)
       enddo
      enddo
 
      call Riemann(u,v,fb,fb_x,xi,dx,nx,ns)

      do i=1,nx
      if(ischeme.eq.1) then
cx     MCV3
       fm_x(i,1)=2.0/dx*fb_x(i)
       fm_x(i,2)=2.0/dx*(3.0*(fb(i+1)-fb(i))-fb_x(i)-fb_x(i+1))/4.0
       fm_x(i,3)=2.0/dx*fb_x(i+1)

        else  if(ischeme.eq.2) then
cx     MCV3_UPCC
       fm_x(i,1)=2.0/dx*(2.0*(f(i,1)+f(i,2))-0.5*(7.0*fb(i)+fb(i+1)))
       fm_x(i,2)=1.0/dx*(f(i,3)-f(i,1))
       fm_x(i,3)=2.0/dx*(-2.0*(f(i,2)+f(i,3))+0.5*(7.0*fb(i+1)+fb(i)))

       else
cx     MCV3_CPCC
       s3=dsqrt(3.0d0)
       fm_x(i,1)=2.0/dx*(s3*(3./4.*f(i,1)+5./6.*f(i,2)-1./12.*f(i,3))
     1 +3./4.*(1.5-s3)*fb(i+1)-3./4.*(1.5+s3)*fb(i))
       fm_x(i,2)=2.0/dx*(f(i,3)-f(i,1))*s3/3.
       fm_x(i,3)=2.0/dx*(s3*(1./12.*f(i,1)-5./6.*f(i,2)-3./4.*f(i,3))
     1 +3./4.*(1.5+s3)*fb(i+1)-3./4.*(1.5-s3)*fb(i))

       endif
      enddo

      call bdc(fm_x,nx,ns)

      return
      end

      subroutine Riemann(u,v,fb,fb_x,xi,dx,nx,ns)
C     
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 u(0:nx+1,ns),v(0:nx+1,ns)
      real*8 f(0:nx+1,ns),fb(0:nx+1),fb_x(0:nx+1)
      real*8 f_l(0:nx+1),f_r(0:nx+1),fx_l(0:nx+1),fx_r(0:nx+1)
      real*8 xi(ns),fe(ns),dx
      integer nx,ns
 
      do i=1,nx
       do k=1,ns
       fe(k)=v(i,k)*u(i,k)
       enddo
       f_l(i)=f_primary(-1.d0,xi,ns,fe)
       f_r(i)=f_primary( 1.d0,xi,ns,fe)
       fx_l(i)=fx_primary(-1.d0,xi,ns,fe)
       fx_r(i)=fx_primary( 1.d0,xi,ns,fe)
      enddo

      call bdc_f(f_l,nx)
      call bdc_f(f_r,nx)
      call bdc_f(fx_l,nx)
      call bdc_f(fx_r,nx)

      do i=1,nx
       if(v(i,1).ge.0.0) then
       fb(i)=f_r(i-1)
       fb_x(i)=fx_r(i-1)
       else
       fb(i)=f_l(i)
       fb_x(i)=fx_l(i)
       endif
      enddo
      call bdc_f(fb,nx)
      call bdc_f(fb_x,nx)

      return
      end

      subroutine bdc(u,nx,ns)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 u(0:nx+1,ns)
         do k=1,ns
            u(0,k) = u(nx,k)
            u(nx+1,k) = u(1,k)
         enddo

      return
      end

      subroutine bdc_f(f,nx)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 f(0:nx+1)
       f(0) = f(nx)
       f(nx+1) = f(1)
       return
       end
 
      real*8 function f_guass_c(x,f1,f2,f3)
      IMPLICIT REAL*8(A-H,O-Z)     
      cc=dsqrt(3.0d0)/2.
      f_guass_c=2./3.*x*(x-cc)*f1-4./3.*(x+cc)*(x-cc)*f2
     1          +2./3.*(x+cc)*x*f3
      return
      end

      real*8 function fx_guass_c(x,f1,f2,f3)
      IMPLICIT REAL*8(A-H,O-Z)     
      cc=dsqrt(3.0d0)/3.
      fx_guass_c=(4./3.*x-cc)*f1-8./3.*x*f2+(4./3.*x+cc)*f3
      return
      end

      real*8 function f_primary(x,xi,ns,f)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 xi(ns),f(ns),p(ns)

      do j=1,ns 
      p(j)=1.0
      k=j
      do j1=1,ns 
      if (j1.ne.k) then
      p(j)=p(j)*(x-xi(j1))/(xi(k)-xi(j1))
      end if
      end do
      end do

      f_primary=0.0
      do j=1,ns 
       f_primary=f_primary+p(j)*f(j)
      end do
      return
      end

      real*8 function fx_primary(x,xi,ns,f)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 xi(ns),f(ns),p(ns),q(ns)

      do j=1,ns 
      p(j)=0.0
      q(j)=1.0
      k=j
      do j1=1,ns 
      if (j1.ne.k) then
      p(j)=p(j)+(x-xi(j1))
      q(j)=q(j)/(xi(k)-xi(j1))
      end if
      end do
      end do

      fx_primary=0.0
      do j=1,ns 
      fx_primary=fx_primary+p(j)/q(j)*f(j)
      end do

      fx_primary=(2.*x-xi(2)-xi(3))/(xi(1)-xi(2))/(xi(1)-xi(3))*f(1)
     1          +(2.*x-xi(1)-xi(3))/(xi(2)-xi(1))/(xi(2)-xi(3))*f(2)
     2          +(2.*x-xi(1)-xi(2))/(xi(3)-xi(1))/(xi(3)-xi(2))*f(3)
      return
      end

