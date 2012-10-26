subroutine init

! Sod shock tube problem (a whimpy test) in 2 or 3 dimensions
! 24jan92 blondin
!=======================================================================

! GLOBALS
use global
use zone
use mpi

IMPLICIT NONE

! LOCALS
INTEGER :: i, j, k, mpierr
REAL :: xmin, xmax, ymin, ymax, zmin, zmax
REAL :: rodt, ridt, xvel, yvel, zvel, width, widthz, widthy
REAL :: dleft, pleft, dright, pright, plane

!--------------------------------------------------------------------------------
! Set up geometry and boundary conditions of grid
!
! Boundary condition flags : nleft, nright
!   = 0  :  reflecting boundary condition
!   = 1  :  inflow/outflow boundary condition
!   = 2  :  fixed inflow boundary condition
!   = 3  :  periodic
! Geometry flag : ngeom                         |  Cartesian:
!   = 0  :  planar                              |    gx = 0, gy = 0, gz = 0
!   = 1  :  cylindrical radial                  |  Cylindrical:
!   = 2  :  spherical   radial             3D= {     gx = 1, gy = 3, gz = 0
!   = 3  :  cylindrical angle                   |
!   = 4  :  spherical polar angle (theta)       |  Spherical:
!   = 5  :  spherical azimu angle (phi)         |    gx = 2, gy = 4, gz = 5

! Define the problem...

ngeomx = 0
ngeomy = 0
ngeomz = 0
nleftx = 0
nrightx= 0
nlefty = 0
nrighty= 0
nleftz = 0
nrightz= 0
xmin   = 0.
xmax   = 1.0
ymin   = 0.
ymax   = 1.0
zmin   = 0.
zmax   = 1.

! If any dimension is angular, multiply coordinates by pi...
if(ngeomy > 2) then
   ymin = ymin * pi
   ymax = ymax * pi
endif
if(ngeomz > 2) then
   zmin = zmin * pi
   zmax = zmax * pi
endif

!======================================================================
! Set up parameters from the problem (Sod shock tube)

pright = 0.1 
dright = 0.125 
pleft  = 1.0 
dleft  = 1.0 
gam    = 5. / 3.
gamm   = gam - 1.0

!=======================================================================
! set time and cycle counters

time   = 0.0
timep  = 0.0
timem  = 0.0
ncycle = 0
ncycp  = 0
ncycd  = 0
ncycm  = 0
nfile  = 1000

! Set up grid coordinates 

call grid(imax,xmin,xmax,zxa,zxc,zdx)
call grid(jmax,ymin,ymax,zya,zyc,zdy)
call grid(kmax,zmin,zmax,zza,zzc,zdz)

if (ndim == 2) zzc(1) = 0.0

!=======================================================================
! Log parameters of problem in history file

if (mype == 0) then
  write (8,*) 'Oblique Sod shock tube '
  if (ndim == 3) then
   write (8,"(' Grid dimensions: ',i4,' x ',i4,' x ',i4)") imax,jmax,kmax
  else
   write (8,"(' Grid dimensions: ',i4,' x ',i4)") imax,jmax
  endif
  write (8,*) 
  write (8,*) 'Adiabatic index, gamma = ', gam
  write (8,*) 'Pressure ratio is ', pright
  write (8,*) 'Density ratio is ', dright
  write (8,*) 
endif

! initialize grid:

do k = 1, ks
 do j = 1, js
  do i = 1, imax
    plane = zxc(i)+zyc(mypey*js+j)+zzc(mypez*ks+k)
    if(plane >= 0.5) then
      zro(i,j,k) = dright
      zpr(i,j,k) = pright
    else
      zro(i,j,k) = dleft
      zpr(i,j,k) = pleft
    endif
    zux(i,j,k) = 0.
    zuy(i,j,k) = 0.
    zuz(i,j,k) = 0.
    zfl(i,j,k) = 0.
  enddo
 enddo
enddo

!########################################################################
! Compute Courant-limited timestep

ridt = 0.

if (ndim == 2) then

 do j = 1, js
  do i = 1, imax
    widthy = zdy(j+mypey*js)
    if(ngeomy > 2) widthy = widthy*zxc(i)
    width  = min(zdx(i),widthy)
    svel = sqrt(gam*zpr(i,j,1)/zro(i,j,1))/width
    xvel = abs(zux(i,j,1)) / zdx(i)
    yvel = abs(zuy(i,j,1)) / widthy
    ridt = max(xvel,yvel,svel,ridt)
  enddo
 enddo

else

 do k = 1, ks
  do j = 1, js
   do i = 1, imax
     widthy = zdy(j+mypey*js)
     widthz = zdz(k+mypez*ks)
     if(ngeomy.gt.2) widthy = widthy*zxc(i)
     if(ngeomz.gt.2) widthz = widthz*zxc(i)
     if(ngeomz.eq.5) widthz = widthz*sin(zyc(j+mypey*js))
     width  = min(zdx(i),widthy,widthz)
     svel = sqrt(gam*zpr(i,j,k)/zro(i,j,k))/width
     xvel = abs(zux(i,j,k)) / zdx(i)
     yvel = abs(zuy(i,j,k)) / widthy
     zvel = abs(zuz(i,j,k)) / widthz
     ridt = max(xvel,yvel,zvel,svel,ridt)
   enddo
  enddo
 enddo

endif

call MPI_ALLREDUCE( ridt, rodt, 1, VH1_DATATYPE, MPI_MAX, MPI_COMM_WORLD, mpierr )
dt = courant / rodt

return
end

!#####################################################################

subroutine grid( nzones, xmin, xmax, xa, xc, dx )

! Create grid to cover physical size from xmin to xmax
! number of physical grid zones is nzones
!
! xa(1) is left boundary location - at xmin
! xa(nzones+1) is right boundary location - at xmax
!----------------------------------------------------------------------
IMPLICIT NONE

! LOCALS
integer :: nzones, n
real, dimension(nzones) :: xa, dx, xc
real :: dxfac, xmin, xmax

!=======================================================================

dxfac = (xmax - xmin) / float(nzones)
do n = 1, nzones
  xa(n) = xmin + (n-1)*dxfac
  dx(n) = dxfac
  xc(n) = xa(n) + 0.5*dx(n)
enddo

return
end
