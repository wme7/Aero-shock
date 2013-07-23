program projectile
  implicit none

  ! define constants
  real, parameter :: g = 9.81
  real, parameter :: pi = 3.1415927

  real :: a, t, u, x, y
  real :: theta, v, vx, vy

  ! Read values for a, t and u from terminal
  read(*,*) a, t, u

  ! Convert angle to radians
  a = a*pi/180.0
  
  x = u*cos(a)*t
  y = u*sin(a)*t-0.5*g*t*t

  vx=u*cos(a)
  vy=u*sin(a)-g*t
  v=sqrt(vx*vx+vy*vy)
  theta=atan(vy/vx)*180.0/pi

  write(*,*) 'x: ',x,' y: ',y
  write(*,*) 'v: ',v,' theta: ',theta

end program projectile
