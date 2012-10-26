subroutine sweepz

! This subroutine performs sweeps in the Y direction, looping over i by k rows.  
!  First data is copied from RECV buffer from sweepy ALLTOALL (where it was 
!     transposed back to the original layout) into send buffer sweepz ALLTOALL
!  The data is then transposed such that an entire row in k is local to each process
!  Finally, data is read in from RECV into 1D arrays for TWO hydro updates.
!-----------------------------------------------------------------------

! GLOBALS
use global
use zone
use sweeps
use mpi

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k, n, ntot, m, mpierr

!-----------------------------------------------------------------------
do j = 1, js
 do m = 1, npey
  do i = 1, isy
   n = i + isy*(m-1)
   do k = 1, ks
     send3(1,j,k,n) = recv2(1,k,i,j,m)
     send3(2,j,k,n) = recv2(2,k,i,j,m)
     send3(3,j,k,n) = recv2(3,k,i,j,m)
     send3(4,j,k,n) = recv2(4,k,i,j,m)
     send3(5,j,k,n) = recv2(5,k,i,j,m)
     send3(6,j,k,n) = recv2(6,k,i,j,m)
   enddo
  enddo
 enddo
enddo

call MPI_ALLTOALL(send3,Za2abuff_size,VH1_DATATYPE,recv3,Za2abuff_size,VH1_DATATYPE,MPI_COMM_COL,mpierr)  

sweep  = 'z'
ngeom  = ngeomz
nleft  = nleftz
nright = nrightz
nmin   = 7
nmax   = kmax + 6

! Now Loop over each column...

do j = 1, js
 do i = 1, isz
   radius = zxc(i+mypez*isz)
   theta  = zyc(j+mypey*js)
   stheta = sin(theta)
   radius = radius * stheta

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do m = 1, npez
    do k = 1, ks
     n = k + ks*(m-1) + 6
     r(n) = recv3(1,j,k,i,m)
     p(n) = recv3(2,j,k,i,m)
     u(n) = recv3(5,j,k,i,m)
     v(n) = recv3(3,j,k,i,m)
     w(n) = recv3(4,j,k,i,m)
     f(n) = recv3(6,j,k,i,m)
    enddo
   enddo

   do k = 1, kmax
     n = k + 6
     xa (n) = zza(k)
     dx (n) = zdz(k)
     xa0(n) = zza(k)
     dx0(n) = zdz(k)
     e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! 1D Hydro update using PPMLR
   call ppmlr

! ##################   REPEAT HYDRO EVOLUTION A SECOND TIME
   do k = 1, kmax
     n = k + 6
     xa (n) = zza(k)
     dx (n) = zdz(k)
     xa0(n) = zza(k)
     dx0(n) = zdz(k)
     e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! 1D Hydro update using PPMLR
   call ppmlr
! ##########################################################

   ! Put updated values into 2D arrays, dropping ghost zones

   do k = 1, kmax
     n = k + 6
     send4(1,j,i,k) = r(n)
     send4(2,j,i,k) = p(n)
     send4(3,j,i,k) = v(n)
     send4(4,j,i,k) = w(n)
     send4(5,j,i,k) = u(n)
     send4(6,j,i,k) = f(n)
   enddo

 enddo
enddo

call MPI_ALLTOALL(send4,Za2abuff_size,VH1_DATATYPE,recv4,Za2abuff_size,VH1_DATATYPE,MPI_COMM_COL,mpierr)  

do k = 1, ks
 do i = 1, isz
  do m = 1, npez
   n = i + isz*(m-1)
   do j = 1, js
     send1(1,k,j,n) = recv4(1,j,i,k,m)
     send1(2,k,j,n) = recv4(2,j,i,k,m)
     send1(3,k,j,n) = recv4(3,j,i,k,m)
     send1(4,k,j,n) = recv4(4,j,i,k,m)
     send1(5,k,j,n) = recv4(5,j,i,k,m)
     send1(6,k,j,n) = recv4(6,j,i,k,m)
   enddo
  enddo
 enddo
enddo

return
end

