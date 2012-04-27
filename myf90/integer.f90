program arithmetic
  implicit none

  ! Define real and integer varaibles
  real :: d, r, rres, k
  integer :: i, j, ires, m

  ! Assign some values
  d = 2.0 ; r = 3.0 ; k = 10.0
  i = 2 ; j = 3 ; m = 10

  ! Now the examples 
  rres = r/d
 
  ! Print result both text and a value.
  ! Note how teh text and value are seprated by a coma:
  write(*,*) 'rres = r/d: ',rres 

  ! now some more examples
  ires = j/i; write(*,*) 'ires = j/i: ',ires
  ires = r/i; write(*,*) 'ires = r/i: ',ires
  rres = r/i; write(*,*) 'rres = r/i: ',rres
  rres = k/r; write(*,*) 'rres = k/r: ',rres
  rres = m/r; write(*,*) 'rres = m/r: ',rres
  ires = k/j; write(*,*) 'ires = k/j: ',ires
  rres = m/j; write(*,*) 'rres = m/j: ',rres

end program arithmetic
