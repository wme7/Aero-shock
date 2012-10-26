subroutine sweepy

! This subroutine performs sweeps in the Y direction, looping over i by k rows.   
!-----------------------------------------------------------------------

! GLOBALS
use global
use zone
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k, n

!-----------------------------------------------------------------------

sweep  = 'y'
ngeom  = ngeomy
nleft  = nlefty
nright = nrighty
nmin   = 7
nmax   = jmax + 6
radius = 1.0

! Now Loop over each column...

do k = 1, kmax
 do i = 1, imax
   if (ngeom>2) radius = zxc(i)

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do j = 1, jmax
     n = j + 6

     r  (n) = zro(i,j,k)
     p  (n) = zpr(i,j,k)
     u  (n) = zuy(i,j,k)
     v  (n) = zuz(i,j,k)
     w  (n) = zux(i,j,k)
     f  (n) = zfl(i,j,k)

     xa (n) = zya(j)
     dx (n) = zdy(j)
     xa0(n) = zya(j)
     dx0(n) = zdy(j)

     p  (n) = max(smallp,p(n))
     e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! Perform 1D hydrodynamic update using PPMLR algorithm
   call ppmlr

   ! Put updated values into 2D arrays, dropping ghost zones

   do j = 1, jmax
      n = j + 6
      zro(i,j,k) = r(n)
      zpr(i,j,k) = p(n)
      zuy(i,j,k) = u(n)
      zuz(i,j,k) = v(n)
      zux(i,j,k) = w(n)
      zfl(i,j,k) = f(n)
   enddo

 enddo
enddo

return
end

