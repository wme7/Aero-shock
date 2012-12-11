      implicit real*8(a-h,o-z)
      parameter (nxd=1000)
      dimension x(nxd), c(nxd), v(nxd)
      dimension u(nxd), T(nxd), et(nxd), ei(nxd), r(nxd), P(nxd)
      dimension fij(nxd,nxd),fxij(nxd,nxd),fij_eq(nxd,nxd)
c	Authored by : Bagus Putra Muljadi (d97543016@ntu.edu.tw)
c 	One-Dimensional shock tube for arbitrary statistics 
c	Based on Semiclassical Boltzmann-BGK
c	Insert it = 0: Maxwell-Boltzmann, 1: Fermi-Dirac, -1: Bose-Einstein
c	CE/SE w-alpha scheme 
      outtime = 0.15
      kgrid = 100
      cn = 0.9
	nv = 201
      ala = 3
	en = 0
	r_time = 0.0005
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
      dl = 10.
      ul = 0.
      Tl = 1.667
      dr = 1.
      ur = 0.
      Tr = 1.333 
c      
      dx = 1.d0/dfloat(kgrid)
c
      hdx = dx/2.d0
      k1 = 2*kgrid + 1
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
         r(j) = dl
         u(j) = ul
         T(j) = Tl
         else
         r(j) = dr
         u(j) = ur
         T(j) = Tr
         end if
      do k = 1, nv
        pp = -(v(k)-u(j))*(v(k)-u(j))/t(j)
        fij(k,j) = r(j)/sqrt(pi*t(j))*exp(pp)
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
      
      dt = cn*dx/v(nv)
      time = time + dt
            if (time.gt.outtime) then
            dtcfl = outtime - (time - dtcfl)
            time = outtime
            dt = dtcfl/cn
            istop = 1
            end if
c         
      do k = 1, nv
	do j = 1,k2 
         pp = -(v(k)-u(j))*(v(k)-u(j))/t(j)
         fij_eq(k,j) = r(j)/sqrt(pi*t(j))*exp(pp)

		if (en.eq.0) then 
		fij(k,j) = fij_eq(k,j)
		end if

      end do
	end do
	
c
      do 500 k = 1, nv
      anu = v(k)*dt/dx
      h = 0.5*abs(anu)
      a1 = (1.d0 - anu)/2.d0
      a2 = (1.d0 + anu)/2.d0
      a3 = (1.d0 - anu**2)/2.d0
      a4 = 1.d0 - h + 2.d0*anu
      a5 = 1.d0 - h - 2.d0*anu
      a6 = 2.d0*(1.d0 + h)
c
c	CE/SE w-alpha
c      
      do 300 i = 1,2
      j1 = 3 - i + (i/2)*2
      j2 = j1 + 2*(kgrid-1)
      do 200 j = j1,j2,2
	if (en.eq.0) then
      fij(k,j) = a1*fij(k,j+1) + a2*fij(k,j-1) 
     +    + a3*(fxij(k,j-1) - fxij(k,j+1))
c     +    + dt*(fij_eq(k,j)-fij(k,j))/(2*r_time)
	else
	fij(k,j) = a1*fij(k,j+1) + a2*fij(k,j-1) 
     +    + a3*(fxij(k,j-1) - fxij(k,j+1))
     +    + dt*(fij_eq(k,j)-fij(k,j))/(2*r_time) 
      end if 
      fp = fij(k,j+1) - a4*fxij(k,j+1)
      fm = fij(k,j-1) + a5*fxij(k,j-1)
      fmbar = (fij(k,j)-fm)/(0.5*a6)
      fpbar = (fp-fij(k,j))/(0.5*a6)
      wm = abs(fpbar)**ala/(abs(fmbar)**ala+abs(fpbar)**ala+1e-20)
      wp = abs(fmbar)**ala/(abs(fmbar)**ala+abs(fpbar)**ala+1e-20)
      fxij(k,j) = (wm*fmbar)+(wp*fpbar)
200   continue
      if (j1.eq.3) goto 250
	fij(k,k2) = fij(k,k3)
	fxij(k,k2) = fxij(k,k3)
      goto 300
250   continue
	fij(k,1) = fij(k,3)
	fxij(k,1) = fxij(k,3)

300   continue
500   continue
c
c	Quadrature Method
c
	
      do j = 1,k2
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
      et(j)= 0.5*se
	ei(j) = se/sr-0.5*u(j)*u(j)
	t(j)  = 4.*ei(j)
      p(j)  = r(j)*t(j) 
      end do
	
c
      if (istop.eq.1) go to 2000
      iter = iter + 1
      write(*,15) time, r(k1/2)
15	format ('elapsed time:',F7.4,4X, 'density at x=0.5:',F7.4)
	
      go to 1000

c      
2000  continue
      do 400 j = 2,k2,2
      x(j) = dfloat(j-1)*hdx
      write (8,*) x(j),r(j),T(j),P(j),u(j)
400   continue
      stop
      end
