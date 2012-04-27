program tripos1
  implicit none

  real :: p1, p2, p3 ,maths
  real :: av1, av2

  ! read int the marks
  read(*,*) p1, p2, p3, maths
 
  !work out two averages
  av1 = p1 + p2 + p3
  av2 = av1 + maths
  av1 = av1 / 3.0; av2 = av2 / 4.0
 
  !use an if statement
  if (av2>av1) then
     write(*,*) 'final average = ',av2
  else
     write(*,*) 'Final average = ',av1
  end if

end program tripos1
