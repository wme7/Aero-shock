      subroutine init 
!     *************** 
!                 
      use mod_wave 
! 
      implicit none 
      double precision,external::func 
      double precision::xic 
      double precision,dimension(ks):: ffunc 
      double precision,dimension(nx,k+1)::phi 
      double precision::cplus,cmins,cave,ciave,cjmp 
      integer::i,j,n1,m 
 
!   
!     computing c 
!     ----------- 
!      if (nproblem.eq.1.or.nproblem.eq.2) then 
       n1=int((a-rl)/dx) 
!       allocate(c(nx)) 
       do i=1,n1 
        c(i)=1.d00 
        c(nx+1-i)=1.d00 
       enddo 
       do i=n1+1,nx-n1 
        c(i)=cc 
       enddo 
!     ------------------------------------ 
      do i=1,nx-1 
       cplus=c(i+1) 
       cmins=c(i) 
! 
       cave=(cplus+cmins)/2.d00 
       ciave=(1.d00/cplus+1.d00/cmins)/2.d00 
       cjmp=cplus-cmins 
! 
       c11(i)=1.d00/(2.d00*cave) 
       c12(i)=cjmp/(4.d00*cave) 
       c22(i)=1.d00/(2.d00*ciave) 
      enddo 
! 
!      computing the alpha coefficients of the 
!      Runge-Kutta methods of order kt 
!      --------------------------------------- 
 
!      allocate(alpha(kt,0:kt-1)) 
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
         do j=1,kt 
          f(j)=0.d0 
         enddo 
 
 
! 
!      nproblem=1,2 
!      Computing the initial data 
!      ------------------------- 
!      if (nproblem.eq.1.or.nproblem.eq.2) then 
! 
!      computing the Guassian points 
!       write(6,*)"Gaussian" 
!       ks=2*k+1 
       call gauleg(-1.d00,1.d00,x,w,ks)   
       call gauleg(-1.d00,1.d00,xt,wt,kt)    
! 
       call legendre(pl,x,ks,k+1) 
       call legendre(plt,xt,kt,kt)        
	   if (ne.eq. 1)  call coefb(kt) 
!       allocate(u(nx,k+1,0:kt)) 
!       allocate(v(nx,k+1,0:kt)) 
!       allocate(fluxu(0:nx)) 
!       allocate(fluxv(0:nx)) 
!      computing the moments of the initial data in `func' 
!      and storing them in `phi' 
        do i=1,nx 
         xic=-a+(2.d00*i-1.d00)*dx/2.d00 
! 
 
 
!        evaluating the function `func' at the quadrature points 
          do m=1,ks 
            ffunc(m)=func(xic+x(m)*dx/2.d00-t0) 
          enddo 
! 
          do j=0,k 
           phi(i,j+1)=0.d00 
           do m=1,ks 
            phi(i,j+1)= phi(i,j+1)+ffunc(m)*pl(m,j+1)*w(m) 
           enddo 
            phi(i,j+1)= phi(i,j+1)*(2.d00*j+1.d00)/2.d00 
          end do
        end do
!

!      setting the initial degrees of freedom
        do i=1,nx
          do j=1,k+1
!       write(6,*)"u",i,j
             u(i,j,0)= phi(i,j)
             v(i,j,0)=-phi(i,j)
          end do
        end do
!      endif 
!
      return
      end subroutine init
