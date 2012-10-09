program SCI_BGK_2D
implicit real(8) (a-h,o-z)
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
real, dimension (igh,igh,imaxgridx,imaxgridy) :: f,feq
real, dimension (imaxgridx,imaxgridy)         :: ux,uy,z,t,r,et,p 
real, dimension (imaxgridx)					  :: x,y  !/imaxgridy assuming symm. domain
real, dimension (igh)						  :: c,v
!	A Direct solver authored by Bagus p. Muljadi (d97543016@ntu.edu.tw)
!	For solving 2D-Riemann problems of semiclassical Boltzmann-BGK Equation
!	using TVD and Discrete Ordinate Method
mx			= 100
my			= 100
nx          = mx+3	
ny          = my+3	
cfl         = 0.8
outtime     = 0.2
theta       = 0.	! maxwellian = 0., fermion = 1., boson = -1.
isolver		= 0		! euler = 0, navier-stokes & beyond = 1
r_time		= 0.01	! relaxation time (tau)
nv			= 20 	! discretized velocity points (uniform for vx and vy)
open (unit = 10, file = 'abscissas.tec', status = 'unknown')
open (unit = 20, file = 'sv_point(1).tec', status = 'unknown')
open (unit = 30, file = 'sv_point(2).tec', status = 'unknown')
open (unit = 60, file = 'results.tec', status = 'unknown')
    write(20,888)
    write(30,888)
	write(60,888)
	write(20,*) 'zone t="sv_point(1)",i=',nx-2,',j=',ny-2
    write(30,*) 'zone t="sv_point(2)",i=',nx-2,',j=',ny-2
    write(60,*) 'zone t="results",i=',nx-2,',j=',ny-2
888 format ('variables = "x","y","n","p","z"')    
call weightings_gauss_hermite (nv,c,v)
call quadrants_vals(z1,ux1,uy1,t1,z2,ux2,uy2,t2,z3,ux3,uy3,t3,z4,ux4,uy4,t4)
call grid2 (nx,ny,dx,dy,x,y)
call initialization (nv,nx,ny,x,y,theta,z1,ux1,uy1,t1,z2,ux2,uy2,t2,z3,ux3,uy3,t3,z4,ux4,uy4,t4,v,z,ux,uy,t,f)
iter  = 1
time  = 0
istop = 0
1000 continue
call caldt (nv,dx,dy,cfl,istop,time,outtime,v,dt,dtdx,dtdy)
call equilibrium (nv,nx,ny,theta,v,z,ux,uy,t,feq)
	if (isolver.eq.0) then
	call eu_projection (nx,ny,nv,f,feq)
	end if
call iteration (nv,nx,ny,dtdx,dtdy,dt,r_time,feq,v,f)
call caldom (nx,ny,nv,c,v,f,r,ux,uy,et)
	if (theta .eq. 0.) go to 1100
call bisect_ztp (nx,ny,theta,r,ux,uy,et,z,t,p)
go to 1300
1100 continue
call maxwellian_ztp (nx,ny,r,ux,uy,et,z,t,p)
   
1300 continue
call wrt_savingpoints (nx,ny,iter,outtime,dt,x,y,r,p,z)
write(*,777) time, r(nx/2,ny/2)
777	format (1X,'Elapsed time:',F7.4,4X, 'Density at x=0.5,y=0.5:',F7.4)
if (istop .eq. 1) goto 2000
iter = iter + 1
goto 1000
2000 continue

call wrt_results (nx,ny,x,y,r,p,z)
stop
end program    
!==================!				   
!	SUBROUTINES	   !
!==================!
subroutine quadrants_vals(z1,ux1,uy1,t1,z2,ux2,uy2,t2,z3,ux3,uy3,t3,z4,ux4,uy4,t4)
implicit none
real(8) z1,ux1,uy1,t1,z2,ux2,uy2,t2,z3,ux3,uy3,t3,z4,ux4,uy4,t4
	Z1   = 0.142
    UX1  = -0.75
    UY1  = -0.5
    T1   = 2.078
    Z2   = 0.4253
    UX2  = -0.75
    UY2  = 0.5
    T2   = 1.1494
    Z3   = 0.142
    UX3  = 0.75
    UY3  = 0.5
    T3   = 2.078
    Z4   = 0.6635
    UX4  = 0.75
    UY4  = -0.5
    T4   = 0.87685
return
end subroutine
!
subroutine grid2 (nx,ny,dx,dy,x,y)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer i,j,nx,ny
real(8) dx,dy
real, dimension (imaxgridx) :: x,y
dx = 1./dfloat(nx-3)
dy = 1./dfloat(ny-3)

x(1)=0.-dx
y(1)=0.-dy
do i = 2,nx
	x(i) = x(i-1) + dx
end do
do j = 2,ny
	y(j) = y(j-1) + dy
end do
return
end subroutine
!
subroutine weightings_gauss_hermite (nv,c,v)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer i,nv
real,dimension (igh):: c,v,gh,w
gh(:) =(/-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843,-2.25497400209,-1.73853771212,-1.2340762154,-0.737473728545,-0.245340708301,0.245340708301,0.737473728545,1.2340762154,1.73853771212,2.25497400209,2.78880605843,3.34785456738,3.94476404012,4.60368244955,5.38748089001/)
w(:)  =(/0.898591961453,0.704332961176,0.62227869619,0.575262442852,0.544851742366,0.524080350949,0.509679027117,0.499920871336,0.493843385272,0.490921500667,0.490921500667,0.493843385272,0.499920871336,0.509679027117,0.524080350949,0.544851742366,0.575262442852,0.62227869619,0.704332961176,0.898591961453/)
!
do i = 1, nv
    c(i)    = w(i)
    v(i)    = gh(i)
    write (10,*) c(i),v(i)    
end do
return
end subroutine
!
subroutine initialization (nv,nx,ny,x,y,theta,z1,ux1,uy1,t1,z2,ux2,uy2,t2,z3,ux3,uy3,t3,z4,ux4,uy4,t4,v,z,ux,uy,t,f)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer k,l,i,j,nx,ny,nv
real(8) z1,ux1,uy1,t1,z2,ux2,uy2,t2,z3,ux3,uy3,t3,z4,ux4,uy4,t4,pp,theta
real, dimension (igh,igh,imaxgridx,imaxgridy) :: f
real, dimension (imaxgridx,imaxgridy) :: z,ux,uy,t
real, dimension (imaxgridx) :: x,y
real, dimension (igh) :: v 
do i = 1,nx
  do j = 1,ny
	if ((y(j) .le. 0.5) .and. (x(i) .ge. 0.5))  then
	  z(i,j) = z4
	  ux(i,j) = ux4
	  uy(i,j) = uy4
	  t(i,j) = t4
	 else if ((y(j) .ge. 0.5) .and. (x(i) .le. 0.5))	then
	  z(i,j) = z2
	  ux(i,j) = ux2
	  uy(i,j) = uy2
	  t(i,j) = t2
	 else if ((y(j) .ge. 0.5) .and. (x(i) .ge. 0.5)) then
	  z(i,j) = z1
	  ux(i,j) = ux1
	  uy(i,j) = uy1
	  t(i,j) = t1
	 else
	  z(i,j) = z3
	  ux(i,j) = ux3
	  uy(i,j) = uy3
	  t(i,j) = t3
	end if
	do k = 1,nv
	  do l = 1,nv
		pp = ( (v(k)-ux(i,j))**2 + (v(l)-uy(i,j))**2 ) / t(i,j)
		f(k,l,i,j) = 1/((exp(pp)/z(i,j))+theta)
	  end do
	end do
  end do
end do
return
end subroutine
!
subroutine equilibrium (nv,nx,ny,theta,v,z,ux,uy,t,f0)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
real,dimension (igh,igh,imaxgridx,imaxgridy) :: f0
real,dimension (imaxgridx,imaxgridy) :: z,ux,uy,t
real,dimension (igh) :: v
integer nv,nx,ny,k,l,i,j
real(8) pp,theta
do i = 1, nx
  do j = 1, ny
	do k = 1,nv
	  do l = 1,nv
		pp = ( (v(k)-ux(i,j))**2 + (v(l)-uy(i,j))**2 ) / t(i,j)
		f0(k,l,i,j) = 1/((exp(pp)/z(i,j))+theta)
	  end do
	end do
  end do
end do
return
end subroutine
!
subroutine eu_projection (nx,ny,nv,f,f0)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer k,l,i,j,nx,ny,nv
real,dimension (igh,igh,imaxgridx,imaxgridy) :: f,f0
do i = 1, nx
  do j = 1, ny
    do k = 1, nv
      do l = 1, nv
        f(k,l,i,j) = f0(k,l,i,j)	
      end do
    end do  
  end do
end do
return
end subroutine
!
subroutine caldt (nv,dx,dy,cfl,istop,time,outtime,v,dt,dtdx,dtdy)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer nv,istop
real(8) dx,dy,cfl,time,outtime,dt,dtdx,dtdy,dtcfl
real,dimension (igh) :: v
dt = min(dx,dy) * cfl/v(nv)
time = time + dt
dtdx = dt/dx
dtdy = dt/dy
!    
if (time .gt. outtime) then
    dtcfl = outtime - (time - dtcfl)  
    time = outtime
    dt = dtcfl/cfl
    dtdx = dtcfl / dx
    dtdy = dtcfl / dy
    istop = 1
end if
return
end subroutine
!
subroutine caldom (nx,ny,nv,c,v,f,r,ux,uy,et)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer k,l,i,j,nx,ny,nv
real(8) sr,sux,suy,se
real, dimension (igh) :: c,v
real, dimension (imaxgridx,imaxgridy) :: r,ux,uy,et
real, dimension (igh,igh,imaxgridx,imaxgridy) :: f
do i = 1, nx
  do j = 1, ny
    sr  = 0
    sux = 0
    suy = 0
    se  = 0
      do k = 1, nv
        do l = 1, nv
         sr  = sr + c(k)*c(l) * f(k,l,i,j)
         sux = sux + c(k)*c(l) * f(k,l,i,j) * v(k)
         suy = suy + c(k)*c(l) * f(k,l,i,j) * v(l)
         se  = se + c(k)*c(l) * f(k,l,i,j) * (0.5 * (v(k)*v(k) + v(l)*v(l)))
        end do
      end do
    r(i,j)    = sr
    ux(i,j)   = sux/sr 
    uy(i,j)   = suy/sr
    et(i,j)   = se       
  end do
end do
return
end subroutine
!
subroutine bisect_ztp (nx,ny,theta,r,ux,uy,et,z,t,p)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer l,i,j,nx,ny
real(8) pi,za,zb,zc,ga1,gb1,gc1,ga2,gb2,gc2,psia,psib,psic,theta
real, dimension (imaxgridx,imaxgridy)::r,ux,uy,et,z,t,p
pi = atan2(1.,1.)*4.
do i = 1, nx
 do j = 1, ny
  za = 0.001
  zb = 0.9
  do while (abs(za-zb) .gt. 0.0001)
    ga1 = 0
    gb1 = 0
    ga2 = 0
    gb2 = 0
        do l = 1, 50
            if (theta .eq. 1.) then      
            ga1 = ga1 + (za**l) * (-1)**(l-1)/l
            gb1 = gb1 + (zb**l) * (-1)**(l-1)/l
            ga2 = ga2 + (za**l) * (-1)**(l-1)/(l**2)
            gb2 = gb2 + (zb**l) * (-1)**(l-1)/(l**2)
            else
            ga1 = ga1 + (za**l) /l
            gb1 = gb1 + (zb**l) /l
            ga2 = ga2 + (za**l) /(l**2)
            gb2 = gb2 + (zb**l) /(l**2)
            end if    
        end do
    psia = 2*et(i,j) - (ga2*(r(i,j)/ga1)**2)/pi - r(i,j)*(ux(i,j)*ux(i,j)+uy(i,j)*uy(i,j))
    psib = 2*et(i,j) - (gb2*(r(i,j)/gb1)**2)/pi - r(i,j)*(ux(i,j)*ux(i,j)+uy(i,j)*uy(i,j))
    zc = (za+zb)/2
    gc1 = 0
    gc2 = 0
        do l = 1, 50
            if (theta .eq. 1.) then
            gc1 = gc1 + (zc**l) * (-1)**(l-1)/l
            gc2 = gc2 + (zc**l) * (-1)**(l-1)/(l**2)
            else 
            gc1 = gc1 + (zc**l)/l
            gc2 = gc2 + (zc**l)/(l**2)
            end if
        end do
    psic = 2*et(i,j) - (gc2*(r(i,j)/gc1)**2)/pi - r(i,j)*(ux(i,j)*ux(i,j)+uy(i,j)*uy(i,j))
    if ((psia*psic) .lt. 0) then
        zb = zc
    else
        za = zc
    end if
  end do
  z(i,j) = zc
  t(i,j) = r(i,j)/(pi*gc1)
  p(i,j) = et(i,j) - 0.5*r(i,j)*(ux(i,j)**2+uy(i,j)**2)
 end do
end do
return
end subroutine
!
subroutine maxwellian_ztp (nx,ny,r,ux,uy,et,z,t,p)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer	nx,ny,i,j
real(8) pi
real,dimension (imaxgridx,imaxgridy) :: r,ux,uy,et,z,t,p
pi = atan2(1.,1.)*4
do i = 1, nx
 do j = 1, ny
  t(i,j) = (2* et(i,j)/r(i,j))-(ux(i,j)**2+uy(i,j)**2)
  z(i,j) = r(i,j)/(pi*t(i,j))
  p(i,j) = r(i,j)*t(i,j)/2
 end do
end do
return
end subroutine
!
subroutine wrt_savingpoints (nx,ny,iter,outtime,dt,x,y,r,p,z)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer iter,i,j,nx,ny
real(8)	dt,sv_point1,sv_point2,outtime
real,dimension (imaxgridx) :: x,y
real,dimension (imaxgridx,imaxgridy) :: r,p,z
sv_point1 = (1./3.)*outtime   
sv_point2 = (2./3.)*outtime
do i = 2, nx-1
  do j = 2, ny-1
   if (iter.eq.floor(sv_point1/dt)) then
     write (20,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   else if (iter.eq.floor(sv_point2/dt)) then
     write (30,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
   end if
  end do 
end do
return
end subroutine
!
subroutine iteration (nv,nx,ny,dtdx,dtdy,dt,r_time,feq,v,f)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer k,l,i,j,nv,nx,ny
real(8)	vxp,vxm,vyp,vym,dtdx,dtdy,dt,r_time
real,dimension (igh) :: v
real,dimension (imaxgridx,imaxgridy) ::fn,fl,fh,gl,gh2
real,dimension (igh,igh,imaxgridx,imaxgridy) :: f,feq
do k = 1, nv
  do l = 1, nv
    vxp = max(v(k),0.)
    vxm = min(v(k),0.) 
    vyp = max(v(l),0.)
    vym = min(v(l),0.)
    call tvdfluxes (nx,ny,k,l,vxp,vxm,vyp,vym,dtdx,dtdy,v,f,fl,fh,gl,gh2)
!	call explicit (k,l,nx,ny,fl,fh,gl,gh2,dt,dtdx,dtdy,feq,f,r_time,fn)
	call implicit_mtd (nx,ny,k,l,vxp,vyp,vxm,vym,fl,fh,gl,gh2,dt,dtdx,dtdy,f,fn)
	do i = 2,nx-1
	  do j = 2,ny-1
		f(k,l,i,j) = fn(i,j)
	    call b_cond (k,l,i,j,nx,ny,fn,f)
	  end do
	end do	  
  end do
end do
return
end subroutine 
!
subroutine tvdfluxes (nx,ny,k,l,vxp,vxm,vyp,vym,dtdx,dtdy,v,f,fl,fh,gl,gh2)
implicit none		 
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
real,dimension (igh) :: v
real,dimension (imaxgridx,imaxgridy) :: tha,thb,phia,phib
real,dimension (imaxgridx,imaxgridy) ::fll,flh,fl,fhl,fhh,fh,gll,glh,gl,ghl,ghh,gh2
real,dimension (igh,igh,imaxgridx,imaxgridy) :: f
integer k,l,i,j,nx,ny
real(8) vxp,vxm,vyp,vym,dtdx,dtdy 
! theta
    do i = 2, nx-1
      do j = 2, ny-1
        if (f(k,l,i+1,j) .eq. f(k,l,i,j)) then
            tha(i,j) =  0
        else
            tha(i,j) = (f(k,l,i-sign(1.,v(k))+1,j) - f(k,l,i-sign(1.,v(k)),j))/(f(k,l,i+1,j)-f(k,l,i,j))
        end if 
        if (f(k,l,i,j+1) .eq. f(k,l,i,j)) then
            thb(i,j) =  0
        else
            thb(i,j) = (f(k,l,i,j-sign(1.,v(l))+1) - f(k,l,i,j-sign(1.,v(l))))/(f(k,l,i,j+1)-f(k,l,i,j))
        end if         
      end do
    end do 
    do i = 1, nx
        thb(i,1)=1
        thb(i,ny)=1
    end do
    do j = 1, ny
        tha(1,j)=1
        tha(nx,j)=1    
    end do
! van leer limiter
    do i = 1,nx
      do j = 1,ny
        if (tha(i,j) .le. 0) then
            phia(i,j)=0
        else
            phia(i,j)= (abs(tha(i,j))+tha(i,j))/(1+abs(tha(i,j)))
        end if
        if (thb(i,j) .le. 0) then
            phib(i,j)=0
        else
            phib(i,j)= (abs(thb(i,j))+thb(i,j))/(1+abs(thb(i,j)))
        end if    
      end do
    end do
    do i = 2, nx-1
      do j = 2, ny-1 
        fll(i,j) = vxp*f(k,l,i-1,j)+vxm*f(k,l,i,j)
        flh(i,j) = 0.5d0*v(k)*(f(k,l,i-1,j)+f(k,l,i,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i,j)-f(k,l,i-1,j)) 
        fl(i,j)  = fll(i,j)+ (phia(i-1,j)*(flh(i,j)-fll(i,j)))  !F(i-1/2,j)       
        fhl(i,j) = vxp*f(k,l,i,j)+vxm*f(k,l,i+1,j)
        fhh(i,j) = 0.5d0*v(k)*(f(k,l,i,j)+f(k,l,i+1,j))-0.5d0*dtdx*v(k)**2*(f(k,l,i+1,j)-f(k,l,i,j)) 
        fh(i,j)  = fhl(i,j)+ (phia(i,j)*(fhh(i,j)-fhl(i,j)))   !F(i+1/2,j) 
        gll(i,j) = vyp*f(k,l,i,j-1)+vym*f(k,l,i,j)
        glh(i,j) = 0.5d0*v(l)*(f(k,l,i,j-1)+f(k,l,i,j))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j)-f(k,l,i,j-1))
        gl(i,j)  = gll(i,j)+ (phib(i,j-1)*(glh(i,j)-gll(i,j))) !F(i,j-1/2)   
        ghl(i,j) = vyp*f(k,l,i,j)+vym*f(k,l,i,j+1)
        ghh(i,j) = 0.5d0*v(l)*(f(k,l,i,j)+f(k,l,i,j+1))-0.5d0*dtdy*v(l)**2*(f(k,l,i,j+1)-f(k,l,i,j))
        gh2(i,j) = ghl(i,j)+ (phib(i,j)*(ghh(i,j)-ghl(i,j)))  !F(i,j+1/2)
      end do
    end do
return
end subroutine
! 
subroutine explicit (k,l,nx,ny,fl,fh,gl,gh2,dt,dtdx,dtdy,feq,f,r_time,fn)
implicit none		 
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer k,l,i,j,nx,ny
real(8) dt,dtdx,dtdy,r_time
real, dimension (imaxgridx,imaxgridy) :: fl,fh,gl,gh2,fn
real, dimension (igh,igh,imaxgridx,imaxgridy) :: f,feq
 do i = 2, nx-1
    do j = 2, ny-1
	  fn(i,j) = f(k,l,i,j) - dtdx*(fh(i,j)-fl(i,j)) - dtdy*(gh2(i,j)-gl(i,j)) + dt*(feq(k,l,i,j)-f(k,l,i,j))/r_time
    end do
 end do
return
end subroutine
!
subroutine implicit_mtd (nx,ny,k,l,vxp,vyp,vxm,vym,fl,fh,gl,gh2,dt,dtdx,dtdy,f,fn)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer	k,l,i,j,nx,ny
real(8) dt,dtdx,dtdy,vxp,vyp,vxm,vym
real, dimension (imaxgridx,imaxgridy) :: fn,rhs,invA,invB,dq1,dq,fh,fl,gl,gh2
real, dimension (igh,igh,imaxgridx,imaxgridy) :: f
call right_hand_side (nx,ny,fh,fl,gh2,gl,dtdx,dtdy,rhs)
call inv_A (nx,ny,vxp,vyp,dt,invA)
call inv_B (nx,ny,vxm,vym,dt,invB)
call multiply (nx,invA,rhs,dq1)
call multiply (nx,invB,dq1,dq)
do i = 2, nx-1
  do j = 2, ny-1
	fn(i,j) = f(k,l,i,j) + dq(i,j)
  end do
end do
return
end subroutine
!
subroutine right_hand_side (nx,ny,fh,fl,gh2,gl,dtdx,dtdy,rhs)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer i,j,nx,ny
real(8) dtdx,dtdy
real,dimension (imaxgridx,imaxgridy) :: rhs,fh,fl,gl,gh2
do i = 2, nx-1
  do j = 2, ny-1
	rhs(i,j) =  - dtdx*(fh(i,j)-fl(i,j)) - dtdy*(gh2(i,j)-gl(i,j))
  end do
end do
return
end subroutine
!
subroutine inv_A (nx,ny,vxp,vyp,dt,u)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer i,j,k,nx,ny
real(8) vxp,vyp,dt,bet
real,dimension (imaxgridx) :: a, b, c, gam
real,dimension (imaxgridx,imaxgridy) :: u, r 
do i = 3,nx-1
	a(i) = dt*(-vxp-vyp)
end do
c(2)=a(3)
do i = 2,nx-1
	b(i) = 1+(dt*(vxp+vyp))
end do
do	i = 3,nx-2
	c(i) = 0
end do
do i = 2,nx-1
	do j = 2,ny-1
		r(i,j) = 0
	end do
	r(i,i) = 1
end do
do  k = 2,ny-1
	bet=b(2)
	u(1,k) = r(1,k)/bet
	do  j = 2,nx-1
		gam(j) = c(j-1)/bet
		bet=b(j)-a(j)*gam(j)
		u(j,k) = (r(j,k)-a(j)*u(j-1,k))/bet
	end do 
	do  j = nx-2,1,-1
		u(j,k) = u(j,k) - gam(j+1)*u(j+1,k)
	end do 
end do
return
end subroutine
!
subroutine inv_B (nx,ny,vxm,vym,dt,u)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer i,j,k,nx,ny
real(8) vxm,vym,dt,bet
real,dimension (imaxgridx) :: a, b, c, gam
real,dimension (imaxgridx,imaxgridy) :: u, r 
do i = 3,nx-1
	a(i) = 0.
end do
do i = 2,nx-1
	b(i) = 1+dt*(-vxm-vym)
end do
do	i = 2,nx-2
	c(i) = dt*(vxm+vym)
end do
a(nx-1)=c(nx-2)
do i = 2,nx-1
	do j = 2,ny-1
		r(i,j) = 0
	end do
	r(i,i) = 1
end do
do  k = 2,ny-1
	bet=b(2)
	u(1,k) = r(1,k)/bet
	do  j = 2,nx-1
		gam(j) = c(j-1)/bet
		bet=b(j)-a(j)*gam(j)
		u(j,k) = (r(j,k)-a(j)*u(j-1,k))/bet
	end do 
	do  j = nx-2,1,-1
		u(j,k) = u(j,k) - gam(j+1)*u(j+1,k)
	end do 
end do
return
end subroutine
!
subroutine multiply (nx,a,b,c)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer i,j,k,nx
real,dimension (imaxgridx,imaxgridy) ::	a,b,c
do i = 2,nx-1
  do j = 2,nx-1
	  c(i,j) = 0.
  end do
end do
do i = 2,nx-1
  do j = 2,nx-1
    do k = 2,nx-1
	  c(i,j)=c(i,j) + a(i,k)*b(k,j)
    end do
  end do
end do
return
end subroutine
!
subroutine b_cond(k,l,i,j,nx,ny,fn,f)
implicit none		  !	This is a sub-subroutine within tvdfluxes
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer k,l,i,j,nx,ny
real,dimension (imaxgridx,imaxgridy) :: fn
real,dimension (igh,igh,imaxgridx,imaxgridy) :: f
	f(k,l,i,1) = fn(i,3)
    f(k,l,i,ny)= fn(i,ny-2)
    f(k,l,1,j) =fn(3,j)
    f(k,l,nx,j)=fn(nx-2,j)
return
end subroutine
!
subroutine wrt_results (nx,ny,x,y,r,p,z)
implicit none
integer igh,imaxgridx,imaxgridy
parameter (igh=20,imaxgridx=400,imaxgridy=400)
integer i,j,nx,ny
real,dimension (imaxgridx) :: x,y
real,dimension (imaxgridx,imaxgridy) :: r,p,z
do i = 2, nx-1
do j = 2, ny-1
   write (60,*) x(i), y(j), r(i,j), p(i,j), z(i,j)
end do 
end do
return
end subroutine