program main
implicit real (a-h,o-z)

real, dimension (1000,1000) :: f
real, dimension (1000) :: x,z,c,v,r,p,u,t,et,ei,fn,fa,fb,fc,av,am,bm,cm,delta_f,gam
real, dimension (20)   :: gh,w

open (unit = 10, file = 'postiter.tec',status = 'unknown')
write(10,*) 'variables = "x","r","p","z","t" '
nx          = 100
nxp1        = nx
ghnc        = 1
cfl         = 4
outtime     = 0.1
it          = 1

if (ghnc .eq. 0) then
nv = 20
else 
nv = 201
end if

pi = atan2(1.,1.)*4.
gh(:) =(/-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843,-2.25497400209,-1.73853771212,-1.2340762154,-0.737473728545,-0.245340708301,0.245340708301,0.737473728545,1.2340762154,1.73853771212,2.25497400209,2.78880605843,3.34785456738,3.94476404012,4.60368244955,5.38748089001/)
w(:)  =(/0.898591961453,0.704332961176,0.62227869619,0.575262442852,0.544851742366,0.524080350949,0.509679027117,0.499920871336,0.493843385272,0.490921500667,0.490921500667,0.493843385272,0.499920871336,0.509679027117,0.524080350949,0.544851742366,0.575262442852,0.62227869619,0.704332961176,0.898591961453/)
if (ghnc .eq. 1) go to 10
    do i = 1, nv
    c(i)    = w(i)
    v(i)    = gh(i)  
    end do
    go to 30
10  continue
    v1  = -20
    v2  = 20
    dv  = (v2-v1)/(nv-1.)
    do i = 1,nv
    v(i) = v1 + dv * (i-1.)
    end do
    
    do i = 2, (nv-1)
    c(i) = 64./45. * dv
    if(mod(i,4).eq.1) c(i) = 28./45. * dv
    if(mod(i,4).eq.3) c(i) = 24./45. * dv
    end do
    c(1)    = 14./45. * dv
    c(nv)   = c(1)
30  continue
    
    ul  = 0.
    tl  = 4.38385
    zl  = 0.2253353
    ur  = 0.
    tr  = 8.972544
    zr  = 0.1204582
    
    dx = 1./dfloat(Nx-1)
    x(1)  = -0.5 *dx
    do j = 2,nxp1
    x(j) = x(j-1) + dx
    end do
    
    do j = 1, nxp1
     if (x(j).le.0.5)then
     u(j) = ul
     t(j) = tl
     z(j) = zl
     else
     u(j) = ur
     t(j) = tr
     z(j) = zr
     end if   
        do k = 1, nv
         pp       = (v(k)-u(j))*(v(k)-u(j))/t(j)
         f(k,j)   = 1/((exp(pp)/z(j)) + it) 
        end do       
    end do

iter  = 1
time  = 0
istop = 0

1000 continue
    
    dt = dx * cfl/v(nv)
    time = time + dt
    dtdx = dt/dx
    
    if (time .gt. outtime) then
    dtcfl = outtime - (time - dtcfl)  
    time = outtime
    dt = dtcfl/cfl
    dtdx = dtcfl / dx
    istop = 1
    end if  
         
    do j = 1, nxp1
    do k = 1, nv
        pp      = (v(k)-u(j)) * (v(k)-u(j)) /t(j)
        f(k,j)   = 1/((exp(pp)/z(j)) + it )
    end do    
    end do

    do k = 1, nv  
             vxp = max(v(k),0.)
             vxm = min(v(k),0.)             
             do j = 2, nxp1-1
             fn(j)  =  - vxp * dtdx * (f(k,j)-f(k,j-1)) - vxm * dtdx * (f(k,j+1)-f(k,j))
             end do 
             fn(1) = fn(2)
             fn(nxp1) = fn(nxp1-1)      
     
        ! [i + dt.dx.a]
        ! the tridiagonal matrix is represented by three vectors of am(nxp1),bm(nxp1) and cm(nxp1); n = discrete space points
        do i = 1,nxp1
            am(i) = -vxp * dt
            bm(i) = 1 + ((-vxm+vxp) * dt)
            cm(i) = vxm * dt
        end do  
        bet = bm(1)
        delta_f(1) = fn(1)/bet
        !decomposition and forward substitution
        do i = 2,nxp1
            gam(i) = cm(i-1)/bet
            bet = bm(i)-am(i)*gam(i)
            delta_f(i) = (fn(i)-am(i)*delta_f(i-1))/bet 
        end do
        !backsubstitution
        do i = nxp1-1,1,-1
            delta_f(i) = delta_f(i)-gam(i+1)*delta_f(i+1)
        end do
        do j = 1,nxp1   
            f(k,j) = f(k,j) + delta_f(j)
        end do
 
    end do    
    do j = 1, nxp1
        sr = 0
        su = 0
        se = 0
        sav= 0
    do k = 1, nv
        sr = sr + c(k) * f(k,j)
        su = su + c(k) * f(k,j) * v(k)
        se = se + c(k) * f(k,j) * (0.5 * v(k) * v(k))
        sav = sav + c(k) * f(k,j) * abs(v(k))
    end do
        r(j)    = sr
        u(j)    = (su/sr)  
        et(j)   = se  
        av(j)   = sav   
    end do   
    if (it.eq.0) go to 1111
    do j = 1, nxp1
        za = 0.0001
        zb = 0.99
    do while (abs(za-zb) .gt. 0.00001)
        ga12 = 0
        gb12 = 0
        ga32 = 0
        gb32 = 0
        
        do l = 1, 50
            if (it.eq.1) then 
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
            if (it.eq.1) then
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
        t(j) = r(j)**2 / (pi*gc12**2 )  
        p(j) = et(j) - 0.5 * r(j) * u(j)**2 
    end do  
   
    go to 1112
1111 continue
      do j = 1, nxp1
            t(j)    = 4*et(j)/r(j) - 2*u(j)*u(j) 
            z(j)    = r(j) / sqrt(pi * t(j))    
            p(j) = et(j) - 0.5 * r(j) * u(j)**2  
      end do  
1112 continue        
        
        write(*,*) iter, x(nxp1/2), r(nxp1/2), z(nxp1/2)
        if (istop .eq. 1) go to 2000
        iter = iter + 1
        go to 1000
2000 continue

    do j = 1, nxp1
    write (10,*) x(j), r(j), p(j), z(j), t(j)
    end do
stop    
end program




