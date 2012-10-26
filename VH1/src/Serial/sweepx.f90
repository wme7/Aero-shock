subroutine sweepx

! This subroutine performs 1D hydro sweeps in the X direction, looping over j by k rows.
!-----------------------------------------------------------------------

! GLOBALS
use zone
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k, n

!-----------------------------------------------------------------------

sweep  = 'x'
ngeom  = ngeomx
nleft  = nleftx
nright = nrightx
nmin   = 7
nmax   = imax + 6

! Now Loop over each row...

do k = 1, kmax
 do j = 1, jmax

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do i = 1,imax
     n = i + 6

     r  (n) = zro(i,j,k)
     p  (n) = zpr(i,j,k)
     u  (n) = zux(i,j,k)
     v  (n) = zuy(i,j,k)
     w  (n) = zuz(i,j,k)
     f  (n) = zfl(i,j,k)

     xa0(n) = zxa(i)
     dx0(n) = zdx(i)
     xa (n) = zxa(i)
     dx (n) = zdx(i)

     p  (n) = max(smallp,p(n))
     e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! Perform 1D hydrodynamics evolution using PPMLR algorithm
   call ppmlr
   
   ! Put updated values back into 3D arrays, dropping ghost zones
   do i = 1, imax
     n = i + 6
     zro(i,j,k) = r(n)
     zpr(i,j,k) = p(n)
     zux(i,j,k) = u(n)
     zuy(i,j,k) = v(n)
     zuz(i,j,k) = w(n)
     zfl(i,j,k) = f(n)
   enddo

 enddo
enddo
       
return
end

