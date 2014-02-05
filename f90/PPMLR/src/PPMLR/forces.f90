subroutine forces(xf,grav,fict)
  
! Calculate both real (ie, gravity and whatnot) 
! and fictitious (ie, coriolis and centrifugal) forces.
!-------------------------------------------------------------
! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: n
REAL :: sinxf0
REAL, DIMENSION(maxsweep) :: grav, fict, xf

!--------------------------------------------------------------
 
if(sweep=='x') then

  if(ngeom==0) then  ! CARTESIAN
  
    do n = nmin-4, nmax+5
      grav(n) = 0.
      fict(n) = 0.
    enddo

  else if(ngeom==1) then  ! CYLINDRICAL R 

    do n = nmin-4, nmax+5
      grav(n) = 0.
      fict(n) = v(n)*v(n)/xf(n)
    enddo

  else if(ngeom==2) then  ! SPHERICAL R

    do n = nmin-4, nmax+5
      grav(n) = 0.0
      if (xf(n)/=0.0) then
        fict(n) = (w(n)*w(n)+v(n)*v(n))/xf(n)
      else
        fict(n) = 0.0
      endif
    enddo

  endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
else if(sweep=='y') then

  if(ngeom==0) then  ! CARTESIAN

    do n = nmin-4, nmax+5
      grav(n) = 0.
      fict(n) = 0.
    enddo

  else if(ngeom==1) then  ! CYLINDRICAL R (2D)

    do n = nmin-4, nmax+5
      grav(n) = 0.
      if (xf(n)/=0.0) then
        fict(n) = v(n)*v(n)/xf(n)
      else
        fict(n) = 0.0
      endif
    enddo

  else if(ngeom==3) then  ! CYLINDRICAL THETA

    do n = nmin-4, nmax+5
      grav(n) = 0.
      fict(n) = -u(n)*w(n) / radius
    enddo

  else if(ngeom==4) then  ! SPHERICAL THETA

    do n = nmin-4, nmax+5
      fict(n) = -u(n)*w(n) / radius
      sinxf0 = sin(xf(n))
      if(abs(sinxf0)>1.0e-5) then
        fict(n) = fict(n)+v(n)*v(n)/radius*cos(xf(n))/sinxf0
      endif
      grav(n) = 0.
   enddo

  endif

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
else

  if(ngeom==0) then  ! CARTESIAN OR CYLINDRICAL

     do n = nmin-4, nmax+5
       grav(n) = 0.
       fict(n) = 0.
     enddo

  else if(ngeom==5) then  ! SPHERICAL PHI

    do n = nmin-4, nmax+5
      grav(n) = 0.0
      fict(n) = -u(n)*v(n)/radius*stheta - u(n)*w(n)/radius*cos(theta) 
    enddo

  endif

endif

return
end

