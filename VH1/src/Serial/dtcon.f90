subroutine dtcon

! set the timestep using various constraints.
! NB - additional constraints may be required depending on the problem!
!-----------------------------------------------------------------------

! GLOBALS
use global
use zone

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k
REAL ::  ridt, dtx, dt3, xvel, yvel, zvel, widthy, widthz

!------------------------------------------------------------------------
!   Hydro constraint on timestep.  Use R*d(theta) if y geometry is angular
ridt = 0.

if(ndim==1) then
  do i = 1, imax
    xvel = abs(zux(i,1,1)) / zdx(i)
    ridt = max(xvel,ridt)
  enddo
else if(ndim==2) then
  do j = 1, jmax
   do i = 1, imax
     widthy = zdy(j)
     if(ngeomy > 2) widthy = widthy*zxc(i)
     xvel = abs(zux(i,j,1)) / zdx(i)
     yvel = abs(zuy(i,j,1)) / widthy
     ridt = max(xvel,yvel,ridt)
   enddo
  enddo
else if(ndim==3) then 
  do k = 1, kmax
   do j = 1, jmax
    do i = 1, imax
      widthy = zdy(j)
      widthz = zdz(k)
      if(ngeomy >  2) widthy = widthy*zxc(i)
      if(ngeomz >  2) widthz = widthz*zxc(i)
      if(ngeomz == 5) widthz = widthz*sin(zyc(j))
      xvel = abs(zux(i,j,k)) / zdx(i)
      yvel = abs(zuy(i,j,k)) / widthy
      zvel = abs(zuz(i,j,k)) / widthz
      ridt = max(xvel,yvel,zvel,ridt)
    enddo
   enddo
  enddo
endif

ridt = max(svel,ridt)
dtx  = courant / ridt     ! global time constraint for given courant parameter
dt3  = 1.1 * dt           ! limiting constraint on rate of increase of dt              
dt   = min( dt3, dtx )    ! use smallest required timestep                                        
      
if (time/dt > 1.e7) then   ! if timestep becomes too small, stop the program!
  write(*,*) 'Timestep has become too small: dt = ',dt
  write(*,*) '                             time = ',time
  call prin('ABORT')
  stop
endif

return
end
