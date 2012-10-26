subroutine images

! Produce a 3D file containing 2D maps of the density over time.
! Data is stored in floating point format using netcdf.  The third 
! dimension uses the unlimited dimension capability of netcdf to 
! grow the third dimension each time a new time step is output to the movie.
!=================================================================================
use global
use zone
use netcdf

IMPLICIT NONE

! LOCALS
CHARACTER(LEN=8) :: coord
CHARACTER(LEN=50):: filename
INTEGER :: nc_stat, ncid, frame_number, rank, nvar, natt, icheck, jcheck
INTEGER :: xDimID, yDimID, tDimID, varID, XScale_varID, YScale_varID, TScale_varID
INTEGER, DIMENSION(3) :: start, count

!======================================================================

! create the netCDF filename and try to open it
filename = 'output/' // trim(prefix) // 'M.nc'
nc_stat = nf90_open(filename, NF90_WRITE, ncid)

if (nc_stat == 2) then   ! file does not exist, so create it and start definitions

  ! create the file and define coordinates/variables
  nc_stat = nf90_create(filename, NF90_Clobber, ncid)
  if (nc_stat /= nf90_NoErr ) print *, nc_stat, NF90_STRERROR(nc_stat)

  ! Initialize Dimensions
  nc_stat = nf90_def_dim(ncid, "x", imax, xDimID)      
  nc_stat = nf90_def_dim(ncid, "y", jmax, yDimID)
  nc_stat = nf90_def_dim(ncid, "time" , NF90_UNLIMITED, tDimID)
  if (nc_stat /= nf90_NoErr ) print *, NF90_STRERROR(nc_stat)

  ! define the coordinate scales as 1D variables
  nc_stat = nf90_def_var(ncid, "x", nf90_float,(/ xDimID /), XScale_varID)
  nc_stat = nf90_def_var(ncid, "y", nf90_float,(/ yDimID /), YScale_varID)
  nc_stat = nf90_def_var(ncid, "time", nf90_float,(/ tDimID /), TScale_varID)
  if (nc_stat /= nf90_NoErr ) print *, NF90_STRERROR(nc_stat)

  ! define the simulation variables
  nc_stat = nf90_def_var(ncid, "Density",nf90_float,(/ xDimID, yDimID, tDimID /), varID)
  if (nc_stat /= nf90_NoErr ) print *, NF90_STRERROR(nc_stat)

  ! take dataset out of definition mode
  nc_stat = NF90_ENDDEF(ncid)

  ! write out coordinate arrays
  nc_stat = nf90_put_var(ncid, XScale_varID, zxc)
  nc_stat = nf90_put_var(ncid, YScale_varID, zyc)
  if (nc_stat /= nf90_NoErr ) print *, NF90_STRERROR(nc_stat)

  frame_number = 1

else

  ! inquire about dimensions to make sure it fits current run
  nc_stat = nf90_inquire(ncid, rank, nvar, natt, tDimID)
  nc_stat = nf90_inquire_dimension(ncid, 1, coord, icheck)
  if(icheck.ne.imax) print *, 'Dim 1 does not equal imax', icheck, imax
  nc_stat = nf90_inquire_dimension(ncid, 2, coord, jcheck)
  if(jcheck.ne.jmax) print *, 'Dim 2 does not equal jmax', jcheck, jmax
  nc_stat = nf90_inquire_dimension(ncid, tDimID, coord, frame_number)
  frame_number = frame_number + 1

endif


start(1) = 1
start(2) = 1
start(3) = frame_number
count(1) = imax
count(2) = jmax
count(3) = 1
nc_stat = nf90_put_var(ncid, 4, zro, start, count)
if (nc_stat /= nf90_NoErr ) print *, NF90_STRERROR(nc_stat)

start(1) = frame_number
nc_stat = nf90_put_var(ncid, 3, time, start)
if (nc_stat /= nf90_NoErr ) print *, NF90_STRERROR(nc_stat)

nc_stat = nf90_close(ncid)

return
end


