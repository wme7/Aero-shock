	program upwind
	implicit real (a-h,o-z)
!	 to solve du/dt + a du/dx = 0
!	 upwind scheme
	real, dimension (1000) :: x,u,un
	open (unit=10,file = 'u.tec',status = 'unknown')
	write(10,*) 'varibles = "x","u"'
	a = 0.5
	ap = max(a,0.)
	am = min(a,0.)	
	cfl = 0.5
	nx = 200
	outtime = 0.2 
!	 grid generation
	dx = 1./dfloat(nx-1)
	x(1) = 0.
	do i = 2,nx
	x(i) = x(i-1) + dx
	end do
!	initial condition
	do i = 1,nx
	if (x(i).le.0.5) then
	u(i) = 1
	else
	u(i) = 0.5
	end if
	end do
!
	iter = 1
	time = 0
	istop = 0
!
1000 continue
	dt = dx * cfl/abs(a)
!
	time = time + dt
	dtdx = dt/dx
!
	if (time.ge.outtime) then
	dtn = outtime - time
	time = outtime
	dt = dtn
	istop = 1
	end if
!	
	do i = 2,nx-1
	un(i) = u(i) - ( dtdx*ap*(u(i)-u(i-1)) + dtdx*am*(u(i+1)-u(i)) )
	end do
	do i = 1,nx
	u(i) = un(i)
	end do
!	boundary condition
	u(1) = un(2)
	u(nx) = un(nx-1)
!
	if (istop.eq.1) goto 2000
	iter = iter + 1
	goto 1000
!
2000	continue
!	
	do j = 1, nx
	write (10,*) x(j),u(j)
	end do 
!
	stop
	end program