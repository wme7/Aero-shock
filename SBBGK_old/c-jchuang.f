      program test
      implicit real*8(a-h,o-z)
c
      Real(8), save, Allocatable, Dimension(:,:) :: f,fx
      Real(8), save, Allocatable, Dimension(:) :: c,v
      Real(8), save, Allocatable, Dimension(:) :: x,r,u,et,ei,p,t
      Real(8), save, Allocatable, Dimension(:) :: fn,fxn,sp,sm,fp,fm
c
      open(11,file='out.dat',status='unknown')
c
      nv = 201
      nx = 101
c
      tau =  .5
      ala =  3
      eps =  0.1
      bet =  0.15
      cfl     = 0.90d0
      outtime = 0.1
c
      nxp1 = nx + 1
      nxm1 = nx - 1
      mx   = nxp1
      ik   = 0.5*mx
c
      pi = atan2(1.,1.)*4.d0
      itwo    = 2
c
      Allocate(f(nv,mx),fx(nv,mx))
      Allocate(c(nv),v(nv))
      Allocate(x(mx),r(mx),u(mx),p(mx),t(mx),ei(mx),et(mx))
      Allocate(fn(mx),fxn(mx),sp(mx),sm(mx),fp(mx),fm(mx))
c
c     Newton-Cotes 4th order quadrature formula
c
      if(mod(nv,4).ne.1) then
         write(*,*)'nv error nv=',nv
      endif
C
      v1       = -20.
      v2       =  20.
      dv       = (v2-v1)/(nv-1.)
      do 100 i = 1,nv
      v(i)     = v1 + dv*(i-1.)
  100 continue
c
      do 110 i = 2,nv-1
      c(i)     = 64./45.*dv
      if(mod(i,4).EQ.1) c(i) = 28./45.*dv
      if(mod(i,4).EQ.3) c(i) = 24./45.*dv
  110 continue
      c(1)     = 14./45.*dv
      c(nv)    = c(1)
c
      do k = 1,nv
      write(11,97)k,v(k),c(k)
      enddo
   97 format(i5,2e15.5)
c
      dl       = 0.445
      ul       = 0.637
      tl       = 13.21
      dr       = 0.5
      ur       = 0.
      tr       = 1.9
c
      dx = 1.d0/dfloat(nx-1)
c
      x(1)   = -0.5d0*dx
      do 120 j =2,nxp1
        x(j) = x(j-1) + dx
  120 continue
c
      do 200 j = 1,nxp1
         if (x(j) .le. 0.5 ) then
            r(j) = dl
            u(j) = ul
            t(j) = tl
         else
            r(j) = dr
            u(j) = ur
            t(j) = tr 
         end if
         do 210 k = 1,nv
         pp       = -(v(k)-u(j))*(v(k)-u(j))/t(j)
         f(k,j)   = r(j)/sqrt(pi*t(j))*exp(pp)
         fx(k,j)  = 0.
  210    continue
  200 continue
c
      iter = 1
      time = 0
      istop = 0
 1000 continue
c
      dt    = dx/v(nv)
      dtcfl = dt*cfl
      time  = time + 0.5*dtcfl
      dtdx  = dtcfl/dx
      if(time .gt. outtime)then
      dtcfl = outtime - (time-0.5*dtcfl)
      time  = outtime
      dt    = dtcfl/cfl
      dtdx  = dtcfl/dx
      istop = 1
      endif
c
c     projection
c
      do 500 j = 1,nxp1
      do 500 k = 1,nv
         pp       = -(v(k)-u(j))*(v(k)-u(j))/t(j)
         f(k,j)   = r(j)/sqrt(pi*t(j))*exp(pp)
  500 continue
c
      eps1= 1.d-60
      ii  = (iter/2)*2
      iod = iter - ii
      m0  = nx + iod
      m1  = m0 - 1
      j1  =  1 - iod
c
      do 300 k = 1,nv
c
      anu  = v(k)*dtdx
      tau0 = abs(anu)
      tau1 = tau0*(1.-tau)+tau
      ap   = 1.d0 + anu
      am   = 1.d0 - anu
c
      do 310 j = 1,m0
      sp(j) = f(k,j) - ap*fx(k,j)
      sm(j) = f(k,j) + am*fx(k,j)
      fm(j) = f(k,j  ) - (1. - 2.*anu - tau1)*fx(k,j  )
  310 continue
      do 315 j = 1,m1
      fp(j) = f(k,j+1) - (1. + 2.*anu - tau1)*fx(k,j+1)
  315 continue
c
      do 320 j = 1,m1
      fn(j+j1) = 0.5d0*(ap*sm(j  ) + am*sp(j+1))
      fa       = 0.5d0*(   sp(j+1) -    sm(j  ))
      fc1      = (fp(j)  - f(k,j))/(1.+tau1)
      fc2      = (f(k,j) - fm(j) )/(1.+tau1)
      fc       = 0.5d0*(fc1 + fc2)
      afc1     = abs(fc1)**ala
      afc2     = abs(fc2)**ala
      fw       = (afc1*fc2 + afc2*fc1)/(afc1 + afc2 + eps1)
      fxn(j+j1)= fa + 2.d0*eps*(fc-fa) + bet*(fw-fc)
  320 continue
C
      if(iod.eq.0)then
      fn(1)    = f(k,1)
      fn(nxp1) = f(k,nx)
      fxn(1)   = fx(k,1)
      fxn(nxp1)= fx(k,nx)
      endif
      do 330 j = 1,m1
      f(k,j) = fn(j)
      fx(k,j) = fxn(j)
  330 continue
c
  300 continue
c
      do 400 j=1,nxp1
      sr     = 0.
      su     = 0.
      se     = 0.
      DO 410 k  = 1,nv
      sr     = sr   + c(k)*f(k,j)
      su     = su   + c(k)*f(k,j)*v(k)
      se     = se   + c(k)*f(k,j)*(0.5*v(k)*v(k))
  410 continue
      r(j)  = sr
      u(j)  = su/sr
      et(j) = 0.5*se
      ei(j) = se/sr-0.5*u(j)*u(j)
      t(j)  = 4.*ei(j)
      p(j)  = r(j)*t(j)
  400 continue
c
      write(*,98)iter,dt,time,x(ik),r(ik),u(ik),t(ik)
c
      if(istop .eq. 1) go to 2000
      iter = iter + 1
      go to 1000
 2000 continue
c       
      open (unit=12,file='b1.tec',status='unknown')
      write(12,96)
   96 format('VARIABLES="X","D","U","T","P"')
      write(12,95)
   95 format('zone')
      sft = 0.
      if(mod(iter,2).eq.1) sft=0.5*dx
      do j =1,nxp1
        write(12,99) x(j)+sft,r(j),u(j),t(j),r(j)*t(j)
      enddo
   99 format(1x,5e15.5)
   98 format(i5,6e15.5)
c
      stop
      end





           
