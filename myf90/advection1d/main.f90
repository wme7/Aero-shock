program main
use extramodule

implicit none
 !
 ! we will try to solve u_t +a u_x = 0
 ! with IC u(0) = sin( x/(2 pi ) )
 !
 real :: u0(10), u(10), u_new(10), x(10)
 integer :: i,j,k
 real, parameter :: pi = 4*atan(1.0)
 real :: dt,dx,tEnd,c
 real, parameter :: cfl = 0.9
 integer :: time,steps

 ! create x domain
 x = linspace(0.0,1.0,10)
 dx = 0.1

 ! load initial condition into domain
 u0 = sin(2*pi*x)

 ! using upwind scheeme to evolve IC
 u = u0 ! load IC into u

 c = 2.0
 tEnd = 1.2
 dt = real(cfl*dx/c) 
 steps = ceiling(tEnd/dt)

do time  = 0,steps
  do i = 2,10
   u_new(i) = u(i) - c*dt/dx*(u(i)-u(i-1))
  end do
 u_new(1) = u(10) ! periodic BC
 u = u_new        ! update info in domain
 print *, 'step: ',time
end do

! write answer
open(unit=12,file='output.plt')
write(12,*) 'variables =x, u'
do i = 1,10
write(12,*) x(i), u(i)
end do  

 print *, u



end program main
