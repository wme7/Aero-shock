subroutine dtcon

! set the timestep using various constraints.
! NB - additional constraints may be required depending on the problem!
!-----------------------------------------------------------------------

! GLOBALS
use global
use zone
use mpi

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k, mpierr
REAL ::  ridt, rodt, dtx, dt3, xvel, yvel, zvel, width, widthy, widthz, sinthe

!------------------------------------------------------------------------
!   Hydro constraint on timestep.
ridt = 0.

do k = 1, ks
 do j = 1, js
  sinthe = sin(zyc(j+mypey*js))
  do i = 1, imax
    widthy = zdy(j+mypey*js)
    widthz = zdz(k+mypez*ks)
    if(ngeomy >  2) widthy = widthy*zxc(i)
    if(ngeomz >  2) widthz = widthz*zxc(i)
    if(ngeomz == 5) widthz = widthz*sinthe
    xvel = abs(zux(i,j,k)) / zdx(i)
    yvel = abs(zuy(i,j,k)) / widthy
    zvel = abs(zuz(i,j,k)) / widthz
    ridt = max(xvel,yvel,zvel,ridt)
  enddo
 enddo
enddo

ridt = max(svel,ridt)
call MPI_ALLREDUCE( ridt, rodt, 1, VH1_DATATYPE, MPI_MAX, MPI_COMM_WORLD, mpierr )

dtx  = courant / rodt     ! global time constraint for given courant parameter
dt3  = 1.1 * dt           ! limiting constraint on rate of increase of dt              
dt   = min( dt3, dtx )    ! use smallest required timestep                                        
      
if (time/dt > 1.e8) then   ! if timestep becomes too small, stop the program!
  if (mype == 0) then
    write(*,*) 'Timestep has become too small: dt = ',dt
    write(*,*) '                             time = ',time
  endif
  call prin
  call MPI_FINALIZE(mpierr)
  stop
endif

return
end

