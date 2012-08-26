      subroutine rk
!     *************
!
      use mod_wave
      implicit none
      integer::ntt,j,l,i,nnt,flag
      double precision::tn1,tn,tn2
            double precision, external:: func

!
!      --- time iterations                       
!      ntt=0
      tn1=t0
      nf=1
!      --- controlling the time iterations

   10 if (dabs(tn1-t).le.1.d-6) return
!
!      --- computing nt, dt, tn, and tn1=tn+dt.
!      tn=t0+ntt*dt
         flag=-1

!        if (dabs(tn1-tf(nf)) .le. 1.d-7) then
!         tn=tn1
!         nf=nf+1
!         dt=dt0
!        else
         dt=dt0
         tn=tn1
         tn1=tn+dt
           if (tn1 .ge. tf(nf)) then
             dt=tf(nf)-tn
             tn1=tf(nf)
             nf=nf+1
           if (ne .eq. 4) flag=1
           endif
!        endif


!      if (ntt.eq.nt) then
!         dt=t-ntt*dt-t0
!         tn1=t
!      else
!         tn1=t0+(ntt+1)*dt
!      endif
!     
!      ntt=ntt+1           
!
!	   	  if (f(1) .gt. .85) write(6,*) tn,f

      do j=1,kt-1
       f(kt-j+1)=f(kt-j)
      enddo
       f(1)=func(-a-tn) 
!	   		  if (f(1) .gt. .95) write(6,*) tn,f

       call interpolation(tn)

!      --- Runge-Kutta inner steps
      do j=0,kt-1
!
!       compute u(.,.,j+1) from u(.,.,j)
!       compute v(.,.,j+1) from v(.,.,j)
!
       call euler(tn,j)
      enddo
!
       if (kt.eq.1) then
        do l=0,k
         do i=1,nx
          u(i,l+1,0)=u(i,l+1,1)     
          v(i,l+1,0)=v(i,l+1,1)
         enddo
        enddo
!
       else
        do l=0,k
         do i=1,nx
          u(i,l+1,0)=alpha(kt,0)*u(i,l+1,0)     
          v(i,l+1,0)=alpha(kt,0)*v(i,l+1,0)     
         enddo
        enddo
! 
       do j=1,kt-2
        do l=0,k
         do i=1,nx
          u(i,l+1,0)= u(i,l+1,0)+alpha(kt,j)*u(i,l+1,j)     
          v(i,l+1,0)= v(i,l+1,0)+alpha(kt,j)*v(i,l+1,j)     
         enddo
        enddo
       enddo
!
       do l=0,k
        do i=1,nx
         u(i,l+1,0)=u(i,l+1,0)+alpha(kt,kt-1)*u(i,l+1,kt)     
         v(i,l+1,0)=v(i,l+1,0)+alpha(kt,kt-1)*v(i,l+1,kt)     
        enddo
       enddo
      endif

!      --- end of the Runge-Kutta inner steps
!      if ((mod(ntt,100).eq.0) .and. (ne .eq. 4)) then
!        write(12,51) 
!          write(12,*) "set(axes,'fontsize', 20)" 
!        write(12,52) 
!          write(12,*) "set(axes,'fontsize', 20)" 
!      end if
!      if (nproblem .le. 2) then
       if (flag .gt. 0) then
       call impvar
!      if ((mod(ntt,20).eq.0) .and. (ne .eq. 4)) call impvar
!      if ((mod(ntt,100).eq.0) .and. (ne .eq. 4)) then
        write(12,51) 
        write(12,21) tn1
        write(12,31)
        write(12,41) nf 
        write(12,52)
        write(12,22) tn1
        write(12,32)   
        write(12,42) nf
!      end if
!       else
!       nnt=int(t/4)
!       if (abs(t-4.d0*nnt) .lt. 4.d0) then
!          if ((mod(ntt,20).eq.0) .and. (ne .eq. 4)) call impvar
!       endif 
       endif

      goto 10
!  
   51  format('figure(1)')
   52  format('figure(2)')
   21  format("title('U at T=",f7.3,"');")
   22  format("title('V at T=",f7.3,"');")
   31  format('axis([-2 2 -0.5 1.5]);')
   32  format('axis([-2 2 -1.5 0.5]);')
   41  format("saveas(gcf, 'u",I4.4,".eps');")
   42  format("saveas(gcf, 'v",I4.4,".eps');")

      end subroutine rk
