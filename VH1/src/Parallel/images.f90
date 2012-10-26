subroutine imagesxy

! gather data from processors to make a 2D slice, append data to a netCDF file
! use ncview to animate by scanning through the time dimension
! this subroutine should only be called by pes with one value of mypez (eg, if mypez==0)
!-----------------------------------------------------------------------------

! GLOBALS

use NETCDF
use global
use zone
use mpi

IMPLICIT NONE

! LOCALS

CHARACTER(LEN=10) :: coord, varname
CHARACTER(LEN=50) :: filename
INTEGER :: i, j, mpierr
INTEGER :: ncstat, ncid, frame_number, rank, nv, natt, icheck, jcheck, count, kmid
INTEGER :: xDimID, yDimID, tDimID, varID, XScale_varID, YScale_varID, TScale_varID
INTEGER, DIMENSION(3) :: start

REAL(4), DIMENSION(imax,jmax/pey) :: slab
REAL(4), DIMENSION(imax,jmax) :: var

!======================================================================
!On each node, Construct variable to be plotted...

varname = 'Density'
kmid = 1
do j = 1, js
 do i = 1, imax
   slab(i,j) = zro(i,j,kmid)
 enddo
enddo

! collect the pieces of the image from all processors
count = imax*js
call MPI_GATHER(slab,count,MPI_REAL,var,count,MPI_REAL,0,MPI_COMM_ROW, mpierr)

if (mypey == 0) then  ! Since there is a unique value of mypez here,
                      !  there is only one pe that opens a file

! create the netCDF filename and try to open it
filename = 'output/' // trim(prefix) // 'XY.nc'
ncstat = nf90_open(filename, NF90_WRITE, ncid)

if (ncstat == 2) then   ! file does not exist, so create it and start definitions

  ! create the file and define coordinates/variables
  ncstat = nf90_create(filename, NF90_Clobber, ncid)
  if (ncstat /= nf90_NoErr ) print *, ncstat, NF90_STRERROR(ncstat)

  ! Initialize Dimensions
  ncstat = nf90_def_dim(ncid, "x", imax, xDimID)      
  ncstat = nf90_def_dim(ncid, "y", jmax, yDimID)
  ncstat = nf90_def_dim(ncid, "time" , NF90_UNLIMITED, tDimID)
  if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  ! define the coordinate scales as 1D variables
  ncstat = nf90_def_var(ncid, "x", nf90_float,(/ xDimID /), XScale_varID)
  ncstat = nf90_def_var(ncid, "y", nf90_float,(/ yDimID /), YScale_varID)
  ncstat = nf90_def_var(ncid, "time", nf90_float,(/ tDimID /), TScale_varID)
  if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  ! define the simulation variables
  ncstat = nf90_def_var(ncid,varname,nf90_float,(/ xDimID, yDimID, tDimID /), varID)
  if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  ! take dataset out of definition mode
  ncstat = NF90_ENDDEF(ncid)

  ! write out coordinate arrays
  ncstat = nf90_put_var(ncid, XScale_varID, zxc)
  ncstat = nf90_put_var(ncid, YScale_varID, zyc)
  if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  frame_number = 1

else

  ! inquire about dimensions to make sure it fits current run
  ncstat = nf90_inquire(ncid, rank, nv, natt, tDimID)
  ncstat = nf90_inquire_dimension(ncid, 1, coord, icheck)
  if(icheck /= imax) print *, 'Dim 1 does not equal imax', icheck, imax
  ncstat = nf90_inquire_dimension(ncid, 2, coord, jcheck)
  if(jcheck /= jmax) print *, 'Dim 2 does not equal jmax', jcheck, jmax
  ncstat = nf90_inquire_dimension(ncid, tDimID, coord, frame_number)
  frame_number = frame_number + 1

endif

start(1) = 1
start(2) = 1
start(3) = frame_number
ncstat = nf90_put_var(ncid, 4, var, start)
if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

start(1) = frame_number
ncstat = nf90_put_var(ncid, 3, time, start)
if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

ncstat = nf90_close(ncid)

endif

return
end

!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------

subroutine imagesxz

! gather data from processors to make a 2D slice, append data to a netCDF file
! use ncview to animate by scanning through the time dimension
! this subroutine should only be called by pes with one value of mypey (eg, if mypey==0)
!-----------------------------------------------------------------------------

! GLOBALS

use NETCDF
use global
use zone
use mpi

IMPLICIT NONE

! LOCALS

CHARACTER(LEN=10) :: coord, varname
CHARACTER(LEN=50) :: filename
INTEGER :: i, k, mpierr
INTEGER :: ncstat, ncid, frame_number, rank, nv, natt, icheck, kcheck, count, jmid
INTEGER :: xDimID, yDimID, tDimID, varID, XScale_varID, YScale_varID, TScale_varID
INTEGER, DIMENSION(3) :: start

REAL(4), DIMENSION(imax,kmax/pez) :: slab
REAL(4), DIMENSION(imax,kmax) :: var

!======================================================================
!On each node, Construct variable to be plotted...

varname = 'Density'
jmid = 1
do k = 1, ks
 do i = 1, imax
   slab(i,k) = zro(i,jmid,k)
 enddo
enddo

! collect the pieces of the image from all processors
count = imax*ks
call MPI_GATHER(slab,count,MPI_REAL,var,count,MPI_REAL,0,MPI_COMM_COL, mpierr)

if (mypez == 0) then  ! since unique value of mypey, there is only one pe to open a file

! create the netCDF filename and try to open it
filename = 'output/' // trim(prefix) // 'XZ.nc'
ncstat = nf90_open(filename, NF90_WRITE, ncid)

if (ncstat == 2) then   ! file does not exist, so create it and start definitions

  ! create the file and define coordinates/variables
  ncstat = nf90_create(filename, NF90_Clobber, ncid)
  if (ncstat /= nf90_NoErr ) print *, ncstat, NF90_STRERROR(ncstat)

  ! Initialize Dimensions
  ncstat = nf90_def_dim(ncid, "x", imax, xDimID)      
  ncstat = nf90_def_dim(ncid, "z", kmax, yDimID)
  ncstat = nf90_def_dim(ncid, "time" , NF90_UNLIMITED, tDimID)
  if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  ! define the coordinate scales as 1D variables
  ncstat = nf90_def_var(ncid, "x", nf90_float,(/ xDimID /), XScale_varID)
  ncstat = nf90_def_var(ncid, "z", nf90_float,(/ yDimID /), YScale_varID)
  ncstat = nf90_def_var(ncid, "time", nf90_float,(/ tDimID /), TScale_varID)
  if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  ! define the simulation variables
  ncstat = nf90_def_var(ncid,varname,nf90_float,(/ xDimID, yDimID, tDimID /), varID)
  if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  ! take dataset out of definition mode
  ncstat = NF90_ENDDEF(ncid)

  ! write out coordinate arrays
  ncstat = nf90_put_var(ncid, XScale_varID, zxc)
  ncstat = nf90_put_var(ncid, YScale_varID, zzc)
  if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

  frame_number = 1

else

  ! inquire about dimensions to make sure it fits current run
  ncstat = nf90_inquire(ncid, rank, nv, natt, tDimID)
  ncstat = nf90_inquire_dimension(ncid, 1, coord, icheck)
  if(icheck /= imax) print *, 'Dim 1 does not equal imax', icheck, imax
  ncstat = nf90_inquire_dimension(ncid, 2, coord, kcheck)
  if(kcheck /= kmax) print *, 'Dim 2 does not equal kmax', kcheck, kmax
  ncstat = nf90_inquire_dimension(ncid, tDimID, coord, frame_number)
  frame_number = frame_number + 1

endif

start(1) = 1
start(2) = 1
start(3) = frame_number
ncstat = nf90_put_var(ncid, 4, var, start)
if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

start(1) = frame_number
ncstat = nf90_put_var(ncid, 3, time, start)
if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

ncstat = nf90_close(ncid)

endif


return
end


