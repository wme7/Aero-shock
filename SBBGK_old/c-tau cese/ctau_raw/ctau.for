      implicit real*8(a-h,o-z)
      parameter (nxd=1000)
      dimension x(nxd), c(nxd), v(nxd)
      dimension z(nxd), u(nxd), T(nxd), et(nxd), r(nxd), P(nxd)
      dimension fij(nxd,nxd),fxij(nxd,nxd),fij_eq(nxd,nxd)
c      
      outtime = 0.1
      k = 100
      cn = 0.7
      nv = 201  
      ala = 4
      it = 1
c
      v1  = -20
      v2  = 20
      dv  = (v2-v1)/(nv-1.)
      do i = 1, nv
      v(i) = v1 + dv * (i-1.)
      end do      
      do i = 2, (nv-1)
      c(i) = 64./45. * dv
      if(mod(i,4).eq.1) c(i) = 28./45. * dv
      if(mod(i,4).eq.3) c(i) = 24./45. * dv
      end do
      c(1)    = 14./45. * dv
      c(nv)   = c(1)
c
      zl = 0.22534d0
      ul = 0.d0
      Tl = 4.3839d0
      zr = 0.12046d0
      ur = 0.d0
      Tr = 8.97254d0  
c      
      dx = 1.d0/dfloat(k)
c
      hdx = dx/2.d0
      k1 = 2*k + 1
      k2 = k1 + 1
      k3 = k1 - 2
      pi = 3.1415926535897932d0
c
      open (unit=8,file='ctau.tec',status='unknown')
      write(8,*) 'variables = "x","r","T","P","z"'
c      
      do 100 j = 1,k2
      x(j) = dfloat(j-1)*hdx
         if (x(j).le.0.5)then
         z(j) = zl
         u(j) = ul
         T(j) = Tl
         else
         z(j) = zr
         u(j) = ur
         T(j) = Tr
         end if
      do k = 1, nv
        pp = (v(k)-u(j))**2 / T(j)
        fij(k,j) = 1/((exp(pp)/z(j))+it)
        fxij(k,j) = 0
      end do
c
100   continue
c
      iter  =   1
      time  =   0
      istop =   0
c
1000  continue
      r_time = 0.001
      dt = cn*dx/v(nv)
      time = time + dt
            if (time.gt.outtime) then
            dtcfl = outtime - (time - dtcfl)
            time = outtime
            dt = dtcfl/cn
            istop = 1
            end if
c         
      do j = 1, k2
      do k = 1, nv
         pp = (v(k)-u(j))**2 / T(j)
         fij_eq(k,j) = 1/((exp(pp)/z(j))+it)
      end do
      end do  
c
      do 500 k = 1, nv
      anu = v(k)*dt/dx
      h = abs(anu)
      a1 = (1.d0 - anu)/2.d0
      a2 = (1.d0 + anu)/2.d0
      a3 = (1.d0 - anu**2)/2.d0
      a4 = 1.d0 - h + 2.d0*anu
      a5 = 1.d0 - h - 2.d0*anu
      a6 = 2.d0*(1.d0 + h)
c      
      do 300 i = 1,2
      j1 = 3 - i + (i/2)*2
      j2 = j1 + 2*(k-1)
      do 200 j = j1,j2,2
      fij(k,j) = a1*fij(k,j+1) + a2*fij(k,j-1) 
     +           + a3*(fxij(k,j-1) - fxij(k,j+1))
     +           + dt*(fij_eq(k,j)-fij(k,j))/(2*r_time)
      fp = fij(k,j+1) - a4*fxij(k,j+1)
      fm = fij(k,j-1) + a5*fxij(k,j-1)
      fmbar = (fij(k,j)-fm)/(0.5*a6)
      fpbar = (fp-fij(k,j))/(0.5*a6)
      wm = abs(fpbar)**ala/(abs(fmbar)**ala+abs(fpbar)**ala+1e-30)
      wp = abs(fmbar)**ala/(abs(fmbar)**ala+abs(fpbar)**ala+1e-30)
      fxij(k,j) = (wm*fmbar)+(wp*fpbar)
200   continue
      if (j1.eq.3) goto 250
      goto 300
250   continue
300   continue
500   continue
c
      do j = 1,k3
      sr = 0
      su = 0
      se = 0
      do k = 1, nv
      sr = sr + c(k) * fij(k,j)
      su = su + c(k) * fij(k,j) * v(k)
      se = se + c(k) * fij(k,j) * 0.5d0 * v(k) * v(k)
      end do
      r(j) = sr
      u(j) = su/sr
      et(j)= se 
      end do
      if (it.eq.0) goto 550
c
      do j = 1, k3
      za = 0.0001
      zb = 0.99
      do while (abs(za-zb) .gt. 0.00001)
      ga12 = 0
      gb12 = 0
      ga32 = 0
      gb32 = 0
        do l = 1, 50
        if (IT.eq.1) then 
        ga12 = ga12 + (za**l)*(-1)**(l-1)/(l**0.5)
        gb12 = gb12 + (zb**l)*(-1)**(l-1)/(l**0.5)
        ga32 = ga32 + (za**l)*(-1)**(l-1)/(l**1.5)
        gb32 = gb32 + (zb**l)*(-1)**(l-1)/(l**1.5)
        else
        ga12 = ga12 + (za**l)/(l**0.5)
        gb12 = gb12 + (zb**l)/(l**0.5)
        ga32 = ga32 + (za**l)/(l**1.5)
        gb32 = gb32 + (zb**l)/(l**1.5)
        end if
        end do
      psia = 2*et(j) - ga32*(r(j)/ga12)**3/(2*pi) - r(j)*u(j)**2
      psib = 2*et(j) - gb32*(r(j)/gb12)**3/(2*pi) - r(j)*u(j)**2
        zc = (za + zb)/2
        gc12 = 0
        gc32 = 0
        gc52 = 0
        do l = 1, 50
        if (IT.eq.1) then
        gc12 = gc12 + (zc**l)*(-1)**(l-1)/(l**0.5)
        gc32 = gc32 + (zc**l)*(-1)**(l-1)/(l**1.5)
        gc52 = gc52 + (zc**l)*(-1)**(l-1)/(l**2.5)
        else
        gc12 = gc12 + (zc**l)/(l**0.5)
        gc32 = gc32 + (zc**l)/(l**1.5)
        gc52 = gc52 + (zc**l)/(l**2.5)
        end if
        end do
      psic = 2*et(j) - gc32*(r(j)/gc12)**3/(2*pi) - r(j)*u(j)**2
        if ((psia*psic) .lt. 0) then
        zb = zc
        else
        za = zc
        end if
      end do
      z(j) = zc
      T(j) = r(j)**2 / (pi*gc12**2 )  
      P(j) = et(j) - 0.5 * r(j) * u(j)**2 
      end do
      go to 600 
c
550   continue
      do j = 1,k3
      T(j) = 4*et(j)/r(j) - 2*u(j)*u(j) 
      z(j) = r(j) / sqrt(pi * t(j))    
      P(j) = et(j) - 0.5 * r(j) * u(j)**2  
      end do 
600   continue    
c
      if (istop.eq.1) go to 2000
      iter = iter + 1
      write(*,*) time, r(k/2)
      go to 1000
c      
2000  continue
      do 400 j = 1,k3,2
      x(j) = dfloat(j-1)*hdx
      write (8,*) x(j),r(j),T(j),P(j),z(j)
400   continue
      stop
      end
