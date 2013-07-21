c WENO central scheme for 1d scalar conservation laws
c coded by Qiu Jianxian  May, 2001         

      dimension rab(3),a(3,3),b(3,3),rf(3),up(2800,2:4)
     &,wf(2800,3,4),w(3),s(3),wab(2800,3),cf(3,3)
     &,pp1(3),pp2(3),rc(3),u(2800),rcn(3),rcp(3),wcp(2800,3)
     &,ff(2800,4),c(3,3),v(2800),gi(2800,4),wcn(2800,3),ffd(2800)

      character*8 char_time,yc
      open(3,file='1d_single.time')
	 open(101,file='table32.err')
      call TIME(char_time)
      write(*,*) 'time' ,char_time
      write(3,*) 'time: ', char_time  
            
      cfl=0.4
      tf=10.0
      eps=1.e-8
      istop=0
      md=4
      n=10 
      dx=2.0/n

c  三点Gauss求积公式的系数
      a1=5.0/18.0
      a2=8.0/18.0
c     a3=a1
	tao1=0.5-sqrt(15.0)/10.0
      tao2=0.5
      tao3=0.5+sqrt(15.0)/10.0

c runge-kutta's coefficient
      a21=0.5
      a32=0.5
      a43=1.0
      b1=1.0/6.0
      b2=1.0/3.0
      b3=1.0/3.0
      b4=1.0/6.0
      c2=0.5
      c3=0.5
      c4=1.0

c coefficient of p(x) at [xj+1/2,xj+1]
      a(1,1)=1.0/16.0
      a(2,1)=-0.25
      a(3,1)=11.0/16.0
      a(1,2)=-1.0/16.0
      a(2,2)=0.5
      a(3,2)=1.0/16.0
      a(1,3)=5.0/16.0
      a(2,3)=0.25
      a(3,3)=-1.0/16.0

c coefficient of p(x) at [xj,xj+1/2]
      b(1,1)=-1.0/16.0
      b(2,1)=0.25
      b(3,1)=5.0/16.0
      b(1,2)=1.0/16.0
      b(2,2)=0.5
      b(3,2)=-1.0/16.0
      b(1,3)=11.0/16.0
      b(2,3)=-0.25
      b(3,3)=1.0/16.0

      c(1,1)=-1.0/24.0
      c(2,1)=1.0/12.0
      c(3,1)=23.0/24.0
      c(1,2)=-1.0/24.0
      c(2,2)=13.0/12.0
      c(3,2)=-1.0/24.0
      c(1,3)=23.0/24.0
      c(2,3)=1.0/12.0
      c(3,3)=-1.0/24.0

c coefficient of p'(x) for f at point xj
      cf(1,1)=0.5
      cf(2,1)=-2.0
      cf(3,1)=1.5
      cf(1,2)=-0.5
      cf(2,2)=0.0
      cf(3,2)=0.5
      cf(1,3)=-1.5
      cf(2,3)=2.0
      cf(3,3)=-0.5

c computation of linear power 
      rf(1)=4.0/24.0
      rf(2)=8.0/12.0
      rf(3)=4.0/24.0
      rc(1)=-9.0/80.0
      rc(2)=49.0/40.0
      rc(3)=-9.0/80.0
      rab(1)=3.0/16.0
      rab(2)=5.0/8.0
      rab(3)=rab(1)
	sigman=0.0
	sigmap=0.0
	do i=1,3
	rcn(i)=(3.0*abs(rc(i))-rc(i))/2.0
	rcp(i)=(3.0*abs(rc(i))+rc(i))/2.0
	sigmap=sigmap+rcp(i)
	sigman=sigman+rcn(i)
	enddo
	do i=1,3
	rcp(i)=rcp(i)/sigmap
	rcn(i)=rcn(i)/sigman
      enddo
      tc=0.0

c loop of different mesh sizes

      do 8 kkkk=1,8
      write(*,*) kkkk
		
c initial condition of cell averages

      do i=md+1,n+md	 
      x=dx*(i-md+0.5)-1.0
      u(i)=(u0(x)-u0(x-dx))/dx
      enddo

c time evoluation loop

      do nn=1,1000000											
      xxx1=1.0
	do i=md+1,n+md+1
c	xxx1=max(xxx1,abs(u(i)))
	enddo
      dt=cfl*dx/xxx1
      if(tc+2.*dt.ge.tf) then
       dt=(tf-tc)/2.0
       istop=1
      endif

c io=0 for first stagger step, io=1 for secondary staggered step    			

      do io=0,1 

c boundary condition for the cell average --- periodic condition 

      do i=1,md	
      u(i)=u(n+i)
      u(n+md+i)=u(md+i)
      enddo

c compt nonlinear smooth indicators for u     

      do i=md-1,n+md+2 
      s(1)=13.0/12.0*(u(i-2)-2.0*u(i-1)+u(i))**2
     &	 +0.25*(u(i-2)-4.0*u(i-1)+3.0*u(i))**2
      s(2)=13.0/12.0*(u(i-1)-2.0*u(i)+u(1+i))**2
     &	 +0.25*(u(i-1)-u(i+1))**2
      s(3)=13.0/12.0*(u(i)-2.0*u(i+1)+u(i+2))**2
     &	 +0.25*(u(i+2)-4.0*u(i+1)+3.0*u(i))**2

      wsum=0.0
      do j=1,3
      w(j)=rab(j)/(eps+s(j))**2
      wsum=wsum+w(j)
      enddo

      do j=1,3
      wab(i,j)=w(j)/wsum
      enddo

      wsum=0.0
      do j=1,3
      w(j)=rcp(j)/(eps+s(j))**2
      wsum=wsum+w(j)
      enddo

      do j=1,3
      wcp(i,j)=w(j)/wsum
      enddo
      wsum=0.0
      do j=1,3
      w(j)=rcn(j)/(eps+s(j))**2
      wsum=wsum+w(j)
      enddo

      do j=1,3
      wcn(i,j)=w(j)/wsum
      enddo
      enddo 

c compute point values

      do i=md-1,n+md+1

      do l=1,3
	  pp1(l)=0.0
      do j=1,3
      pp1(l)=pp1(l)+c(j,l)*u(i-4+j+l)
      enddo
      enddo

      ttt=0.0
      do j=1,3
      ttt=ttt+pp1(j)*(wcp(i,j)*sigmap-wcn(i,j)*sigman)
      enddo

c v(i)=point value of u at xi at tn
	  v(i)=ttt   
      ff(i,1)=f(v(i))
      enddo

c boundary condition for v and f

      do i=1,md
      v(i)=v(n+i)
      v(n+md+i)=v(md+i)
      ff(i,1)=f(v(i))
      ff(n+md+i,1)=f(v(n+md+i))
      enddo

c 4th order Runge-Kutta

	  do l_rk=1,4  

c compt nonlinear weight for f     

      do i=md-1,n+md+2 
      s(1)=13.0/12.0*(ff(i-2,l_rk)-2.0*ff(i-1,l_rk)+ff(i,l_rk))**2
     &+0.25*(ff(i-2,l_rk)-4.0*ff(i-1,l_rk)+3.0*ff(i,l_rk))**2
      s(2)=13.0/12.0*(ff(i-1,l_rk)-2.0*ff(i,l_rk)+ff(1+i,l_rk))**2
     &	 +0.25*(ff(i-1,l_rk)-ff(i+1,l_rk))**2
      s(3)=13.0/12.0*(ff(i,l_rk)-2.0*ff(i+1,l_rk)+ff(i+2,l_rk))**2
     &+0.25*(ff(i+2,l_rk)-4.0*ff(i+1,l_rk)+3.0*ff(i,l_rk))**2
c	s(1)=(ff(i-2,l_rk)-2.0*ff(i-1,l_rk)+ff(i,l_rk))**2
c      s(2)=(ff(i-1,l_rk)-2.0*ff(i,l_rk)+ff(1+i,l_rk))**2
c      s(3)=(ff(i,l_rk)-2.0*ff(i+1,l_rk)+ff(i+2,l_rk))**2
       wsum=0.0
       do j=1,3
       w(j)=rf(j)/(eps+s(j))**2
       wsum=wsum+w(j)
       enddo

       do j=1,3
       wf(i,j,l_rk)=w(j)/wsum
       enddo

       enddo 

c compute the derivative f_x at j+1/2

       do i=md-1,n+md+1 

       do l=1,3
       w(l)=0.0
       do k=1,3
       w(l)=w(l)+cf(k,l)*ff(i-4+l+k,l_rk)
       enddo 
       enddo 

       xxx1=0.0
       do l=1,3
       xxx1=xxx1+w(l)*wf(i,l,l_rk)
       enddo 
       ffd(i)=xxx1/dx

       enddo 

c compute the right-hand-side gi

       do i=md,n+md+1 
       gi(i,l_rk)=-ffd(i)
       if(l_rk.eq.1.or.l_rk.eq.2) then
	    up(i,l_rk+1)=v(i)+0.5*dt*gi(i,l_rk)
	    ff(i,l_rk+1)=f(up(i,l_rk+1))
       elseif(l_rk.eq.3) then 
        up(i,l_rk+1)=v(i)+dt*gi(i,l_rk)
        ff(i,l_rk+1)=f(up(i,l_rk+1))
	   endif
       enddo 

       if(l_rk.lt.4) then
c set boundary condition
	    do i=1,md	
         up(i,l_rk+1)=up(n+i,l_rk+1)
         ff(i,l_rk+1)=f(up(i,l_rk+1))
         up(n+md+i,l_rk+1)=up(md+i,l_rk+1)
         ff(n+md+i,l_rk+1)=f(up(n+md+i,l_rk+1))
        enddo
       endif
  
      enddo 
c end of RK

c comp of value u at gauss point 

      do i=md,n+md+1 
      up(i,2)=v(i)+dt*(bb1(tao1)*gi(i,1)+bb2(tao1)*gi(i,2)
     &	+bb3(tao1)*gi(i,3)+bb4(tao1)*gi(i,4))
      up(i,3)=v(i)+dt*(bb1(tao2)*gi(i,1)+bb2(tao2)*gi(i,2)
     &	+bb3(tao2)*gi(i,3)+bb4(tao2)*gi(i,4))
      up(i,4)=v(i)+dt*(bb1(tao3)*gi(i,1)+bb2(tao3)*gi(i,2)
     &	+bb3(tao3)*gi(i,3)+bb4(tao3)*gi(i,4))
      enddo  

      do 501 i=md+1-io,n+md-io

      do l=1,3
      pp1(l)=0.0
      pp2(l)=0.0
      do j=1,3
      pp1(l)=pp1(l)+a(j,l)*u(i-4+j+l)
      pp2(l)=pp2(l)+b(j,l)*u(i-3+j+l)
      enddo
      enddo

      ttt=0.0
      do j=1,3
      ttt=ttt+pp2(j)*wab(i+1,j)+pp1(j)*wab(i,j)
      enddo

      v(i)=ttt-dt/dx*((f(up(i+1,2))-f(up(i,2)))*a1+(f(up(i+1,3))
     &-f(up(i,3)))*a2+(f(up(i+1,4))-f(up(i,4)))*a1)
 501  continue

      do 502 i=md+1-io,n+md-io
 502  u(i+io)=v(i)

      tc=tc+dt
	  enddo   
c end of io

      if(istop.eq.1) goto 1000
      enddo  
c  enddo of nnn

 1000 do i=1,md	 
      u(i)=u(n+i)
      u(n+md+i)=u(md+i)
      enddo
      
      do i=md,n+md
      do l=1,3
	  pp1(l)=0.0
      do j=1,3
      pp1(l)=pp1(l)+c(j,l)*u(i-4+j+l)
      enddo
      enddo

      ttt=0.0
      do j=1,3
      ttt=ttt+pp1(j)*(wcp(i,j)*sigmap-wcn(i,j)*sigman)
      enddo

c v(i)=point value of u at xi at tn
	  v(i)=ttt   
      enddo
      yc='b'// char(kkkk+48) //'.plt'
      open(1,file=yc)		
      do kk0=md,md+n
      write(1,*) dx*(kk0-md)-1.0,v(kk0)
      enddo
      close(1)
	 error1=0.0
      error2=0.0
      do kk0=md+1,md+n 
      xxx=dx*(kk0-md)-1.0
      error1=error1+abs(fu1(xxx)-v(kk0))
      error2=max(error2,abs(fu1(xxx)-v(kk0)))
      enddo
      error1=error1*dx
      if(kkkk.eq.1) write(101,103) n,error1,error2
      write(*,*) error1,error2
      if(kkkk.gt.1) then
       rr1=log(er1/error1)/log(2.0)
       rr2=log(er2/error2)/log(2.0)
       write(101,102) n,error1,rr1,error2, rr2
       write(*,*) n,rr1,rr2
      endif
      er1=error1
      er2=error2
      tc=0.0
      istop=0
      n=n*2
	  dx=2.0/n
  8	  continue	
      close(101)
      call TIME(char_time)
      write(3,*)  'time: ', char_time 
      write(*,*) 'time' ,char_time
102    format(i6,1x,2('&',1x, es12.5e2,1x,'&' 1x,f8.5 ,1x)
     &,'\\',1x,'\hline')
103    format(i6,1x,2('&',1x,es12.5E2,1x,'&',1x),'\\',1x,'\hline')

      stop
      end


       function u0(x)
c primitive function of initial condition

      pi=4.0*atan(1.0)
c      u0=-cos(pi*x)/pi 
      u0=3.0/8.0*x-sin(2.0*pi*x)/4.0/pi
     &+sin(4.0*pi*x)/32.0/pi
      return
      end

         function fu1(x)
c initial condition
      pi=4.0*atan(1.0)
c      fu1=sin(pi*x)
      fu1=(sin(pi*x))**4   
      return
      end

         function f(x)
      f=x
      return
      end

	   function bb1(x)
 	bb1=0.666666666666666666666666667*x-1.5
	bb1=bb1*x*x+x
	return
	end
	  	function bb2(x)
        bb2=(-0.666666666666666666666666667*x+1.0)*x*x
	return
	end
	   function bb3(x)
        bb3=(-0.666666666666666666666666667*x+1.0)*x*x
	return
	end
	 function bb4(x)
        bb4=(0.6666666666666666666666666667*x-0.5)*x*x
	return
	end

