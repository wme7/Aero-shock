	program tvd
	implicit real (a-h,o-z)
!	to solve du/dt + a du/dx = 0
!	TVD scheme
	real, dimension (1000) :: x,u,r,theta,phi,FN,F 
!	real, dimension (1000) :: fll,flh,fl
!	real, dimension (1000) :: frl,frh,fr
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
			f(i) = 1.
		else
			f(i) = 0.5
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
       
	       ANU  = a*DTDX 
             DO j = 2,NX-1
                IF (F(J+1).EQ.F(J)) THEN
                THETA(j) = 0.
                ELSE
                THETA(j) = (F(J-SIGN(1.,ANU)+1)-F(J-SIGN(1.,ANU)))/(F(J+1)-F(J))
             END IF        
             END DO
             THETA(1) = 0
             THETA(NX) = 0
             DO j = 1,nx
                if (THETA(j) .LE. 0) then
                phi(j) = 0          
                else
                phi(J) = (ABS(THETA(J))+THETA(J))/(1+ABS(THETA(J)))
                end if
             END DO
             DO J=3,(NX-2)
             FN(J)  =  F(J) - DTDX * ( 0.5D0*a*( ( F(J) + F(J+1) ) - ( F(J-1)+F(J) ) ) - 0.5D0 * SIGN(1.,ANU)*a * ( ( F(J+1)-F(J) ) - ( F(J)-F(J-1) ) )+ 0.5D0 * ( (phi(J) * (SIGN(1.,ANU)-ANU) *a*(F(J+1)-F(J)) ) - (phi(J-1) * (SIGN(1.,ANU)-ANU)*a*( F(J)-F(J-1) ) ) ) ) !+ (DT/VIS(J))*(FEQ(J)-F(J))
            end do
			Do J=1,nx
			F(J) = FN(J)  
            F(1) = fn(4)
            F(2) = fn(3) 
            F(NX) = fn(nx-3)      
            F(NX-1) = fn(nx-2) 
			end do

	if (istop.eq.1) goto 2000
	iter = iter + 1
	goto 1000
!
2000	continue
!
!	Compute real solution for this case:
	do i = 1,nx
		if (x(i).le.(0.5d0+a*time)) then
		r(i) = 1.
		else
		r(i) = 0.5
		end if
	end do
!
!	Write to data file	
	do j = 1, nx
		write (10,*) x(j),f(j),r(j)!,theta(j)
	end do 
!
!	Ending program
	stop
	end program