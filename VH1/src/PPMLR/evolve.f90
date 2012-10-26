subroutine evolve( umid, pmid )

! Use umid and pmid from Riemann solver to update velocity, density, and total energy.
! Physical zones are from nmin to nmax.  Zone boundary numbers run from nmin to nmax+1
!-----------------------------------------------------------------------
! GLOBALS
use global
use sweeps
      
IMPLICIT NONE

! LOCALS
INTEGER :: n
REAL :: dtheta
REAL, DIMENSION(maxsweep) :: umid, pmid, amid, uold, xa1, dvol1, upmid, dm, dtbdm
REAL, DIMENSION(maxsweep) :: grav0, grav1, xa2, fict0, fict1, xa3
REAL, PARAMETER :: third = 1.0 / 3.0

!------------------------------------------------------------------------

call volume (nmin, nmax, ngeom, radius, xa , dx , dvol1 )

! grid position evolution

do n = nmin-3, nmax + 4
  dm   (n) = r(n) * dvol1(n)
  dtbdm(n) = dt / dm(n)
  xa1  (n) = xa(n)
  xa   (n) = xa(n) + dt * umid(n) / radius
  upmid(n) = umid(n) * pmid(n)
enddo
 
xa1(nmin-4) = xa(nmin-4)
xa1(nmax+5) = xa(nmax+5)
      
do n = nmin-4, nmax+5
  xa2(n)   = xa1(n) + 0.5*dx(n)
  dx (n)   = xa(n+1) - xa(n)
  xa3(n)   = xa (n) + 0.5*dx(n)
enddo

! Calculate forces using coordinates at t and at t+dt, note that
!   fictitious forces at t+dt depend on updated velocity, but we ignore this

call forces(xa2,grav0,fict0)
call forces(xa3,grav1,fict1)

! Calculate dvolume and average area based on geometry of sweep

call volume (nmin, nmax, ngeom, radius, xa , dx , dvol )

select case (ngeom)
  case(1)
    amid = 0.5*(xa + xa1)
  case (2) 
    amid = (xa-xa1)*(third*(xa-xa1)+xa1)+xa1**2
  case (4)
   do n = nmin-3, nmax+4
     dtheta  = xa(n) - xa1(n)
     if(dtheta == 0.0) then
       amid(n) = sin(xa(n))
     else
       amid(n) = (cos(xa1(n))-cos(xa(n)))/dtheta
     endif
   enddo
  case default
   amid = 1.0
end select

do n = nmin-3, nmax+3

! density evolution. lagrangian code, so all we have to do is watch the change in the geometry.

  r(n) = r(n) * ( dvol1(n) / dvol(n) )
  r(n) = max(r(n),smallr)

! velocity evolution due to pressure acceleration and forces.

  uold (n) = u(n)
  u(n) = u(n) - dtbdm(n)*(pmid(n+1)-pmid(n))*0.5*(amid(n+1)+amid(n)) + 0.5*dt*(grav0(n)+fict0(n)+grav1(n)+fict1(n))

! total energy evolution

  e(n) = e(n) - dtbdm(n)*(amid(n+1)*upmid(n+1) - amid(n)*upmid(n)) + 0.5*dt*(uold(n)*grav0(n) + u(n)*grav1(n))
  q(n) = e(n) - 0.5*(u(n)**2+v(n)**2+w(n)**2)
  q(n) = max(q(n),smallp/(gamm*r(n)))

enddo

return
end 
