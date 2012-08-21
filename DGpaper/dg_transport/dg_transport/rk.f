      subroutine rk
c     *************
c
      include 'com'
c
c      --- time iterations                       
      ntt=0
      tn1=0.d00
c
c      --- controlling the time iterations
   10 if (dabs(tn1-t).le.1.d-6) then
        call impvar
        return
      endif
c
c      --- computing nt, dt, tn, and tn1=tn+dt.
      tn=ntt*dt
      if (ntt.eq.nt) then
         dt=t-ntt*dt
         tn1=t
      else
         tn1=(ntt+1)*dt
      endif     
      ntt=ntt+1                
c
c      --- Runge-Kutta inner steps
      do j=0,kt-1
c
c       compute u(.,.,j+1) from u(.,.,j)
c
       call euler(tn,j)
c
      enddo
c
       if (kt.eq.1) then
        do l=0,k
         do i=1,nx
          u(i,l+1,0)=u(i,l+1,1)     
         enddo
        enddo
c
       else
        do l=0,k
         do i=1,nx
          u(i,l+1,0)=alpha(kt,0)*u(i,l+1,0)
         enddo
        enddo
c 
       do j=1,kt-2
        do l=0,k
         do i=1,nx
          u(i,l+1,0)= u(i,l+1,0)+alpha(kt,j)*u(i,l+1,j)
         enddo
        enddo
       enddo
c
       do l=0,k
        do i=1,nx
         u(i,l+1,0)=u(i,l+1,0)+alpha(kt,kt-1)*u(i,l+1,kt)
        enddo
       enddo
      endif
c
c      --- end of the Runge-Kutta inner steps
      goto 10
c             
      end
