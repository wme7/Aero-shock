program vertical
  !
  ! Vertical motion under gravity
  !
  implicit none
   
  ! Acceleration due to gravity
  real, parameter :: g = 9.81
  
  ! Variables
  real :: s ! displacement
  real :: t ! time
  real :: u ! initial speed (m/s)

  ! Set values of variables
  t = 6.0
  u = 60.0

  ! Calculate displacement
  s = u*t-g*(t**2)/2

  ! Output results
  write(*,*) 'time=',t,'   Displacement=',s

end program vertical
