subroutine sweepx2

! This subroutine performs 1D hydro sweeps in the X direction, 
! looping over js by ks rows.  Data is pulled from the RECV buffer
! used in the sweepy ALLTOALL call.
!-----------------------------------------------------------------------

! GLOBALS
use zone
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k, n, m

!-----------------------------------------------------------------------

sweep  = 'x'
ngeom  = ngeomx
nleft  = nleftx
nright = nrightx
nmin   = 7
nmax   = imax + 6

! Now Loop over each row...

do k = 1, ks
 do j = 1, js

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do m = 1, npey
    do i = 1, isy
      n = i + isy*(m-1) + 6
      r(n) = recv2(1,k,i,j,m)
      p(n) = recv2(2,k,i,j,m)
      u(n) = recv2(3,k,i,j,m)
      v(n) = recv2(4,k,i,j,m)
      w(n) = recv2(5,k,i,j,m)
      f(n) = recv2(6,k,i,j,m)
    enddo
   enddo

   do i = 1,imax
     n = i + 6
     xa0(n) = zxa(i)
     dx0(n) = zdx(i)
     xa (n) = zxa(i)
     dx (n) = zdx(i)
     p  (n) = max(smallp,p(n))
     e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! Do 1D hydro update using PPMLR
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

