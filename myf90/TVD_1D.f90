	program tvd
	implicit real (a-h,o-z)
!	to solve du/dt + a du/dx = 0
!	TVD scheme
	real, dimension (1000) :: x,u,un,r,theta,phi,fll,flh,fl,frl,frh,fr
	open (unit=10,file = 'u.tec',status = 'unknown')
	write(10,*) 'variables = "x","u","r"'
	cfl     = 0.2
	nx      = 200
	outtime = 0.5
!
!	Grid generation (dx: constant)
	dx = 1./dfloat(nx-1)
	x(1) = 0.
	do i = 2,nx
		x(i) = x(i-1) + dx
	end do
!
!	initial condition (u: data)
	do i = 1,nx
		if (x(i).le.0.5) then
			u(i) = 1.
		else
			u(i) = 0.5
		end if
	end do
!
!	Advection Speed
	a  = 0.5
	ap = max(a,0.)
	am = min(a,0.)	
!
	iter  = 1
	time  = 0
	istop = 0
!
!	Time step (dt: depends on "a" value)
	dt = dx * cfl/abs(a)
	dtdx = dt/dx ! value of charactesristics 
	v=a*dtdx
!
1000 continue
!	Time
	time = time + dt
!	but, 
	if (time.ge.outtime) then
		dtn = outtime - time
		time = outtime
		dt = dtn
		istop = 1
	end if
!
!   Flux Limiter: Van Leer limiter
!
!	Step 1: Ratio of consecutive gradients (theta)
    do i = 2, nx-1	 !98 iterations
        if (u(i+1) .eq. u(i)) then
            theta(i) = 0 !
        else
            theta(i) = (u(i-sign(1.,a)+1.)-u(i-sign(1.,a)))/(u(i+1.)-u(i))
        end if
	end do
    theta(1)=0	! the missing 2	iterations
    theta(nx)=0
!    
!	Step 2: Van Leer Limiter Function (phi)
	do i = 3,nx-2
        if (theta(i) .le. 0) then
            phi(i) = 0
        else
            phi(i) = (abs(theta(i))+theta(i))/(1+abs(theta(i)))
        end if
    end do
!
!	Fluxes Definition:
    do i = 2, nx-1
		!right
		frl(i) = ap*u(i) + am*u(i+1)
        frh(i) = 0.5d0*a*(u(i) + u(i+1))-0.5d0*dtdx*a**2*(u(i+1) - u(i))
		fr(i)  = frl(i) + (phi(i)*(frh(i) - frl(i)))  !F(i+1/2)
		!left
        fll(i) = ap*u(i-1) + am*u(i)
        flh(i) = 0.5d0*a*(u(i-1) + u(i))-0.5d0*dtdx*a**2*(u(i) - u(i-1)) 
		fl(i)  = fll(i) + (phi(i-1)*(flh(i) - fll(i))) !F(i-1/2)
		!new data
		un(i) = u(i) - dtdx*(fr(i) - fl(i))
		u(i)  = un(i)
		u(1)  = un(2)
		u(nx) = un(nx-1)
	end do
	!
	if (istop.eq.1) goto 2000
		iter = iter + 1
	goto 1000
!
2000	continue
!
!	Compute real solution for this case:
	do i = 1,nx
		if (x(i).le.(0.5+a*time)) then
			r(i) = 1.
		else
			r(i) = 0.5
		end if
	end do
!
!	Write to data file	
	do j = 1, nx
		write (10,*) x(j),u(j),r(j),theta(j)
	end do 
!
!	Ending program
	stop
	end program