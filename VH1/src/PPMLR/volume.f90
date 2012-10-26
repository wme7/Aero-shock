subroutine volume(nmin, nmax, ngeom, radius, xa, dx, dvol)

! Calculate zone volume and average area based on geometry of sweep
!---------------------------------------------------------------
! GLOBALS
use sweepsize

IMPLICIT NONE

! LOCALS
INTEGER :: n, nmin, nmax, ngeom
REAL :: radius
REAL, DIMENSION(maxsweep) :: xa, dx, dvol
REAL, PARAMETER :: third = 1.0/3.0

!------------------------------------------------------------------------

select case (ngeom)
  case (0) 
    dvol = dx
  case (1)
    dvol = dx*(xa+0.5*dx)
  case (2)
    dvol = dx*(xa*(xa+dx)+dx*dx*third)
  case (3) 
    dvol = dx*radius
  case (4) 
    do n = nmin-3, nmax+4
      dvol (n) = (cos(xa(n))-cos(xa(n+1)))*radius
    enddo
  case (5)
    dvol = dx*radius
  case default 
    write(*,*) 'Geometry', ngeom, ' not implemented.'
    stop
end select               

return
end
