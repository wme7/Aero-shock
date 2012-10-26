program VHone

! This is a 1D version of VH-1, stripped down to the bare bones.
!--------------------------------------------------------------------------
!
!          VV              VV   HH       HH         111
!           VV            VV    HH       HH        1111
!            VV          VV     HH       HH          11
!             VV        VV      HH       HH          11
!              VV      VV       HHHHHHHHHHH   ==     11
!               VV    VV        HH       HH          11
!                VV  VV         HH       HH          11
!                 VVVV          HH       HH          11
!                  VV           HH       HH        111111
!
! 
!                        VIRGINIA HYDRODYNAMICS #1
!
!--------------------------------------------------------------------------
!  The Virginia Numerical Bull Session ideal hydrodynamics PPMLR
!
!  Version 3.1 July 2008
!    Cleaned up code, minor changes to algorithm, 2/3D output
!  Version 3.0 c. 2003
!    Converted to F90 free-format, use netcdf for binary output  
!  Version 2.0 July 1997
!    Many small modifications, mainly brought BCs into only one call to sweepbc 
!  Version 1.0 October 1991 
!    This version was removed from the ``loop'' in the continued development
!    of VH-1 on Sep 1 1990, but has followed a similar evolution until this
!    release date.
! 
!-------------------------------------------------------------------------------
! GLOBAL VARIABLES declared in modules

use global
use sweeps

IMPLICIT NONE

! LOCAL VARIABLES used only in this main program

CHARACTER(len=8) :: todayis 
CHARACTER(len=2) :: rstrt 
CHARACTER(len=4) :: tmp1
CHARACTER(len=50):: hstfile
CHARACTER(len=50):: filename

INTEGER :: ncycend, nprin, nmovie, ndump, imax, n
REAL :: endtime, tprin, tmovie, dxfac, xmin, xmax, ridt, xvel, dt3, dtx

! ADDITIONAL VARIABLES NEEDED FOR CURRENT SIMULATION (Sod Shock)
INTEGER :: nmid
REAL :: pright, pleft, dright, dleft

NAMELIST / hinput / rstrt, prefix, ncycend, ndump, nprin, nmovie, endtime, tprin, tmovie

!-----------------------------------------------------------------------------------------
! Begin by reading input file (namelist) for job control values

open (unit=15,file='indat',status='old',form='formatted')
 read (15,nml=hinput)
close(15)

! Name and open a file for recording the simulation metadata

hstfile  = 'output/' // trim(prefix) // '.hst'
open (unit=8,file=hstfile,form='formatted')
call date_and_time(todayis)
write (8,*) 'History File for VH-1 simulation run on ', todayis(5:6), ' / ', todayis(7:8), ' / ', todayis(1:4)
write (8,*) 

!============================================================
! CREATE GRID

! Set some flags for geometry and boundary conditions

sweep  = 'x'   ! default is x-sweep (only option for 1D)
ngeom  = 0     ! Cartesian geometry
nleft  = 0     ! Reflecting at xmin
nright = 0     ! Reflecting at xmax

! Create a grid of imax zones, making room for 6 'ghost zones' on each end

imax = 100           ! total number of real zones on grid
nmin = 7             ! first real zone
nmax = imax + 6      ! last real zone
xmin = 0.            ! x value at left edge of grid
xmax = 1.0           ! x value at right edge of grid

! Initialize grid coordinates: xa0(n) is coordinate of left edge of zone n

dxfac = (xmax - xmin) / float(imax)  ! width of each zone
do n = nmin, nmax
  xa0(n) = xmin + (n-nmin)*dxfac 
  dx0(n) = dxfac
enddo

!============================================================
! INITIAL CONDITIONS

! Set up parameters for the problem; in this case a Sod shock tube

pright = 0.1 
dright = 0.125 
pleft  = 1.0 
dleft  = 1.0 
gam    = 1.4   ! We always need a ratio of specific heats, gamma

gamm   = gam - 1.0 

! initialize grid for Sod shock problem:

nmid = imax/2+6
do n = nmin, nmid      ! left half of grid has density, pressure given by dleft, pleft
 r(n) = dleft
 p(n) = pleft
 u(n) = 0.0            ! velocity is zero everywhere
 v(n) = 0.0            ! note that we have to carry around the transverse
 w(n) = 0.0            ! velocities even though this is a 1D code
 f(n) = 0.0            ! set initial flattening to zero
enddo

do n = nmid+1, nmax    ! right half of grid has density, pressure given by dright, pright
 r(n) = dright
 p(n) = pright
 u(n) = 0.0            ! velocity is zero everywhere
 v(n) = 0.0            ! note that we have to carry around the transverse
 w(n) = 0.0            ! velocities even though this is a 1D code
 f(n) = 0.0            ! set initial flattening to zero
enddo

! Write out initial data to a file
nfile = 0
write(tmp1,"(i4)") nfile + 1000 ; nfile = nfile + 1
filename = 'output/' // trim(prefix) // tmp1 // '.dat'
open(unit=3,file=filename,form='formatted')
do n = nmin, nmax
  write(3, 1003) xa0(n), r(n), p(n), u(n)
enddo
close(3)

! Log parameters of problem in history file

write (8,*) 'Sod shock tube in 1 dimension'
write (8,"('Grid dimensions: ',i4)") imax
write (8,*) 
write (8,"('Ratio of specific heats = ',f5.3)") gam
write (8,"('Pressure Ratio = ',f5.3)") pright/pleft
write (8,"('Density Ratio = ',f5.3)")  dright/dleft
write (8,*) 

! Compute initial time step based on Courant condition to ensure stability

ridt = 0.
do n = nmin, nmax
 svel = sqrt(gam*p(n)/r(n))/dx0(n)
 xvel = abs(u(n)) / dx0(n)
 ridt = max(xvel,svel,ridt)
enddo
dt = courant / ridt

! set time and cycle counters

time   = 0.0
timep  = 0.
ncycle = 0
ncycp  = 0

!############################################################################
!  MAIN COMPUTATIONAL LOOP

do while (ncycle < ncycend)

  if ( time + dt  >  endtime ) then 
    dt = endtime - time     ! set dt so time lands on endtime
    ncycend = ncycle-1      ! set ncycend so simulation stops after this step
  endif

  ! construct total energy; reset Lagrangian coordinates to Eulerian grid
  !    The Eulerian coordinates, xa0(n), will stay fixed; 
  !    The Lagrangian coordinates, xa(n), will move with the flow each time step
  do n = nmin, nmax
    e (n) = p(n)/(r(n)*gamm) + 0.5*(u(n)**2+v(n)**2+w(n)**2)
    xa(n) = xa0(n)
    dx(n) = dx0(n)
  enddo
  
  ! perform 1D hydro update with PPMLR algorithm
  call ppmlr

  ! update time and cycle counters
  time   = time  + dt
  timep  = timep + dt
  ncycle = ncycle + 1
  ncycp  = ncycp  + 1

  !============================================================
  ! TIME STEP CONTROL
  
  ridt = 0.         ! searching for maximum, so start out with zero
  do n = nmin, nmax
   svel = sqrt(gam*p(n)/r(n))/dx0(n)
   xvel = abs(u(n)) / dx0(n)
   ridt = max(xvel,svel,ridt)
  enddo
  dtx  = courant / ridt     ! global time constraint for given courant parameter
  dt3  = 1.1 * dt           ! limiting constraint on rate of increase of dt              
  dt   = min( dt3, dtx )    ! use smallest required timestep                                        

  !============================================================
  ! DATA OUTPUT
  
  if ( ncycp >= nprin .or. timep >= tprin ) then    ! Print out a data file

    write(tmp1,"(i4)") nfile + 1000 ; nfile = nfile + 1
    filename = 'output/' // trim(prefix) // tmp1 // '.dat'
    write(8,6001) trim(prefix) // tmp1, time, ncycle

    open(unit=3,file=filename,form='formatted')
     do n = nmin, nmax
       write(3, 1003) xa0(n), r(n), p(n), u(n)
     enddo
    close(3)
  
    timep = 0.
    ncycp = 0
    
  endif

enddo
!                           END OF MAIN LOOP
!#########################################################################

close( 8 )

 1003 format(1pe13.5,4e13.5)
 6001 format('Wrote ',a18,' to disk at time =',1pe12.5,' (ncycle =', i7,')')

stop      
end 
