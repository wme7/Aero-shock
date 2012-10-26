subroutine sweepz

! This subroutine performs 1D hydro sweeps in the Z direction, looping over i by j rows.
!-----------------------------------------------------------------------

! GLOBALS
use global
use zone
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k, n

!-----------------------------------------------------------------------

sweep  = 'z'
ngeom  = ngeomz
nleft  = nleftz
nright = nrightz
nmin   = 7
nmax   = kmax + 6

! Now Loop over each row...

do j = 1, jmax
 do i = 1, imax
   radius = zxc(i)
   theta  = zyc(j)
   stheta = sin(theta)
   radius = radius * stheta

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do k = 1,kmax
     n = k + 6

     r  (n) = zro(i,j,k)
     p  (n) = zpr(i,j,k)
     u  (n) = zuz(i,j,k)
     v  (n) = zux(i,j,k)
     w  (n) = zuy(i,j,k)
     f  (n) = zfl(i,j,k)

     xa (n) = zza(k)
     dx (n) = zdz(k)
     xa0(n) = zza(k)
     dx0(n) = zdz(k)

     p  (n) = max(smallp,p(n))
     e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! Perform 1D hydrodynamic update using PPMLR algorithm
   call ppmlr

   ! Put updated values back into 3D arrays, dropping ghost zones
   do k = 1, kmax
     n = k + 6
     zro(i,j,k) = r(n)
     zpr(i,j,k) = p(n)
     zuz(i,j,k) = u(n)
     zux(i,j,k) = v(n)
     zuy(i,j,k) = w(n)
     zfl(i,j,k) = f(n)
   enddo

 enddo
enddo
 
return
end
