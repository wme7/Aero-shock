subroutine sweepy

! This subroutine performs sweeps in the Y direction, looping over is by ks rows.   
!   After call to MPI_ALLTOALL, data is read out of recv buffer into 1D arrays for hydro
!   If only two dimensions, perform a second hydro update
!   After hydro update, data is put back into a send buffer to transpose back to original grid.
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

call MPI_ALLTOALL(send1,Ya2abuff_size,VH1_DATATYPE,recv1,Ya2abuff_size,VH1_DATATYPE,MPI_COMM_ROW,mpierr)  

sweep  = 'y'
ngeom  = ngeomy
nleft  = nlefty
nright = nrighty
nmin   = 7
nmax   = jmax + 6

! Now Loop over each column...

do k = 1, ks
 do i = 1, isy
   radius = zxc(i+mypey*isy)

   ! Put state variables into 1D arrays, padding with 6 ghost zones
   do m = 1, npey
    do j = 1, js
     n = j + js*(m-1) + 6
     r(n) = recv1(1,k,j,i,m)
     p(n) = recv1(2,k,j,i,m)
     u(n) = recv1(4,k,j,i,m)
     v(n) = recv1(5,k,j,i,m)
     w(n) = recv1(3,k,j,i,m)
     f(n) = recv1(6,k,j,i,m)
    enddo
   enddo

   do j = 1, jmax
     n = j + 6
     xa (n) = zya(j)
     dx (n) = zdy(j)
     xa0(n) = zya(j)
     dx0(n) = zdy(j)
     e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
   enddo

   ! 1D Hydrodynamic update using PPMLR
   call ppmlr

   if (ndim==2) then ! #########   REPEAT HYDRO EVOLUTION A SECOND TIME

    do j = 1, jmax
      n = j + 6
      xa (n) = zya(j)
      dx (n) = zdy(j)
      xa0(n) = zya(j)
      dx0(n) = zdy(j)
      e  (n) = p(n)/(r(n)*gamm)+0.5*(u(n)**2+v(n)**2+w(n)**2)
    enddo
    call ppmlr

   endif    ! ##########################################################

   ! Put updated values into 2D arrays, dropping ghost zones

   do j = 1, jmax
     n = j + 6
     send2(1,k,i,j) = r(n)
     send2(2,k,i,j) = p(n)
     send2(3,k,i,j) = w(n)
     send2(4,k,i,j) = u(n)
     send2(5,k,i,j) = v(n)
     send2(6,k,i,j) = f(n)
   enddo

 enddo
enddo

call MPI_ALLTOALL(send2,Ya2abuff_size,VH1_DATATYPE,recv2,Ya2abuff_size,VH1_DATATYPE,MPI_COMM_ROW,mpierr)  

return
end

