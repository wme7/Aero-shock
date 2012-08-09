program upwind
  !
  ! 1D Linear Advection Equation.
  !
  ! du/dt + a du/dx = 0
  !
  ! Subroutine for solving using Upwind Method.
  ! by Manuel Diaz, manuel.ade'at'gmail.com 
  ! Institute of Applied Mechanics, 2012.08.08


implicit real (a-h, o-z)
real, dimension (500) :: x,t,u_0,u,u_next
integer t_step,nx,i,j,k
real a,dt,dx,dtdx,cfL

! Parameters
a = 0.6
nx = 100
dx = (2.0-1.0)/nx
cfl = 0.5
t_step = 20
dt = cfl*dx/a
dtdx = dt/dx

! Discretize the Domain
x(1) = 1
do i = 1,nx
   x(i+1) = x(i)+dx
end do

t(1) = 0
do i = 1,t_step
   t(i+1) = t(i)+dt
end do 

! The initial condition
do i =1,nx
  if (x(i).le.1.5) then
    u_0(i) = 2.0
  else
    u_0(i) = 1.0
  end if
end do

! Main Loop
do i = 1,nx
  u(i) = u_0(i)
end do

do k = 1,t_step
  do j = 2,nx
     u_next(j) = u(j) - dtdx*a*(u(j)-u(j-1))
  end do
  ! BC
  u_next(1) = u_next(2)
  ! Update info
  do i = 1,nx  
     u(i) = u_next(i)
  end do
end do

! Write down our result
open (10, file='result.dat', status='new') 
  write(10,*) 'Variables x, u'
  do i = 1,nx
     write (10,*) x(i), u(i)
  end do
close(10)

end program upwind
