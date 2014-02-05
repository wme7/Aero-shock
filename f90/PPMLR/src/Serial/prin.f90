subroutine prin

! For 1D, write out separate ascii file for each call to prin
! For 2D/3D, write out one netcdf data file containing all variables, with time as
!  the last (expandable) dimension
!---------------------------------------------------------------------------
use NETCDF
use global
use zone

IMPLICIT NONE

! LOCALS
character(LEN=4) :: tmp1
character(LEN=8) :: coord
character(LEN=50):: filename, filename1d
 
integer :: i, nc_stat, ncid, frame_number, natt, nvar, rank
integer :: xDimID, yDimID, zDimID, tDimID
integer :: density_varID, pressure_varID, XVelocity_varID, YVelocity_varID, ZVelocity_varID
integer :: XScale_varID, YScale_varID, ZScale_varID, TScale_varID
integer, dimension(ndim+1) :: start, count, DimIDs

!------------------------------------------------------------------------------

if (ndim == 1) then         ! Keep 1D output simple, just write out in ascii...

  ! Create filename from integer nfile (in global.mod) and prefix such that filename
  ! looks like prefx1000.dat where 1000 is the value of nfile
  write(tmp1,"(i4)") nfile + 1000
  filename1d = 'output/' // trim(prefix) // tmp1 // '.dat'

  open(unit=3,file=filename1d,form='formatted')
   do i = 1, imax
     write(3, 1003) zxa(i), zro(i,1,1),zpr(i,1,1), zux(i,1,1)
   enddo
  close(3)
  1003 format(1pe13.5,4e13.5)

  write(8,6001) trim(prefix)//tmp1//'.dat', time, ncycle

!#################################################################################

else    ! use netCDF calls to write out 2D or 3D arrays

 ! create the netCDF file
 filename = 'output/' // trim(prefix) // '.nc'
 nc_stat = nf90_open(filename, NF90_WRITE, ncid)
 
 if (nc_stat == 2) then   !  file does not exist, so create it and define variables

 nc_stat = nf90_create(filename, nf90_Clobber, ncid); if (nc_stat /= nf90_NoErr ) print *, NF90_STRERROR(nc_stat)

 ! Initialize Dimensions
             nc_stat = nf90_def_dim(ncid, "x", imax, xDimID)      
             nc_stat = nf90_def_dim(ncid, "y", jmax, yDimID)
 if(ndim==3) nc_stat = nf90_def_dim(ncid, "z", kmax, zDimID)
             nc_stat = nf90_def_dim(ncid, "time" , NF90_UNLIMITED, tDimID)

 ! define the coordinate scales as 1D variables
             nc_stat = nf90_def_var(ncid, "x",    nf90_float, xDimID, XScale_varID)
             nc_stat = nf90_def_var(ncid, "y",    nf90_float, yDimID, YScale_varID)
 if(ndim==3) nc_stat = nf90_def_var(ncid, "z",    nf90_float, zDimID, ZScale_varID)
             nc_stat = nf90_def_var(ncid, "time", nf90_float, tDimID, TScale_varID)
 
             DimIDs(1) = xDimID
             DimIDs(2) = yDimID
 if(ndim==3) DimIDs(ndim) = zDimID
             DimIDs(ndim+1) = tDimID
 
 ! define the simulation variables
              nc_stat = nf90_def_var(ncid, "Density",   nf90_float, DimIDs, density_varID)
              nc_stat = nf90_def_var(ncid, "Pressure",  nf90_float, DimIDs, pressure_varID)
              nc_stat = nf90_def_var(ncid, "XVelocity", nf90_float, DimIDs, XVelocity_varID)
              nc_stat = nf90_def_var(ncid, "YVelocity", nf90_float, DimIDs, YVelocity_varID)
  if(ndim==3) nc_stat = nf90_def_var(ncid, "ZVelocity", nf90_float, DimIDs, ZVelocity_varID)

 ! label the velocity components as vector using a variable attribute
              nc_stat = nf90_put_att(ncid, XVelocity_varID, "vector", "x")
              nc_stat = nf90_put_att(ncid, YVelocity_varID, "vector", "y")
  if(ndim==3) nc_stat = nf90_put_att(ncid, ZVelocity_varID, "vector", "z")

 ! define some global attributes?

 nc_stat = nf90_put_att(ncid, nf90_global, "coordinatesystem", ngeomx) 

 nc_stat = NF90_ENDDEF(ncid) ! take dataset out of definition mode
 
             nc_stat = nf90_put_var(ncid, XScale_varID, zxc)
             nc_stat = nf90_put_var(ncid, YScale_varID, zyc)
 if(ndim==3) nc_stat = nf90_put_var(ncid, ZScale_varID, zzc)

 frame_number = 1

else

  nc_stat = nf90_inquire(ncid, rank, nvar, natt, tDimID)
  nc_stat = nf90_inquire_dimension(ncid, tDimID, coord, frame_number)

              nc_stat = nf90_inq_varid(ncid, "Density",   density_varID)
              nc_stat = nf90_inq_varid(ncid, "Pressure",  pressure_varID)
              nc_stat = nf90_inq_varid(ncid, "XVelocity", XVelocity_varID)
              nc_stat = nf90_inq_varid(ncid, "YVelocity", YVelocity_varID)
  if(ndim==3) nc_stat = nf90_inq_varid(ncid, "ZVelocity", ZVelocity_varID)

  frame_number = frame_number + 1

endif

start(1) = 1
start(2) = 1
count(1) = imax
count(2) = jmax
start(ndim+1) = frame_number
count(ndim+1) = 1

if (ndim==3) then
  start(ndim) = 1
  count(ndim) = kmax
endif

 ! Populate variables (dataset)
             nc_stat = nf90_put_var(ncid, density_varID,   zro, start, count)
             nc_stat = nf90_put_var(ncid, pressure_varID,  zpr, start, count)
             nc_stat = nf90_put_var(ncid, XVelocity_varID, zux, start, count)
             nc_stat = nf90_put_var(ncid, YVelocity_varID, zuy, start, count)
 if(ndim==3) nc_stat = nf90_put_var(ncid, ZVelocity_varID, zuz, start, count)
 
 start(1) = frame_number
 nc_stat = nf90_put_var(ncid, ndim+1, time, start)

 nc_stat = nf90_close(ncid)
 write(8,6000) trim(prefix) // '.nc', time, ncycle

endif

nfile = nfile + 1


 6000 format('Wrote time slice to ',a20,' at time =',1pe12.5,' (ncycle =', i7,')')
 6001 format('Wrote ',a16,' to disk at time =',1pe12.5,' (ncycle =', i7,')')

return
end

