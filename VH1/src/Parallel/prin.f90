subroutine prin

! collect all rows of data onto mypey=0, and then 
! write out pez netcdf data files containing all variables (plus time as an attribute).
! These nc data files can then be merged later using mergeslabs.f90.
!---------------------------------------------------------------------------

use NETCDF
use global
use zone
use mpi

IMPLICIT NONE

! LOCALS
CHARACTER(LEN=1) :: char
CHARACTER(LEN=4) :: tmp1, tmp2
CHARACTER(LEN=50) :: filename
CHARACTER(LEN=16), DIMENSION(nvar) :: varname
 
INTEGER, PARAMETER :: nvarout = 5   ! # of arrays to write into netcdf file
INTEGER :: i,j,k,ncstat, ncid, gathbuffer_size, mpierr, m, jsk, nv, nvars
INTEGER :: xDimID, yDimID, zDimID
INTEGER, DIMENSION(ndim) :: DimIDs
INTEGER, DIMENSION(nvar) :: varID
INTEGER :: XScale_varID, YScale_varID, ZScale_varID

REAL(4), DIMENSION(kmax/pez) :: zshort
REAL(4), DIMENSION(imax,jmax,kmax/pez) :: prin_buff
REAL(4), DIMENSION(imax,nvarout,jmax/pey,kmax/pez) :: send_buff
REAL(4), DIMENSION(imax,nvarout,jmax/pey,kmax/pez,pey) :: recv_buff

!------------------------------------------------------------------------------
! everybody loads up a send buffer; data is gathered on mypey=0 procs.

varname(1) = 'Density'
varname(2) = 'Pressure'
varname(3) = 'Xvelocity'
varname(4) = 'Yvelocity'
varname(5) = 'Zvelocity'

do k = 1, ks
 do j = 1, js
  do i = 1, imax
    send_buff(i,1,j,k) = zro(i,j,k)
    send_buff(i,2,j,k) = zpr(i,j,k)
    send_buff(i,3,j,k) = zux(i,j,k)
    send_buff(i,4,j,k) = zuy(i,j,k)
    send_buff(i,5,j,k) = zuz(i,j,k)
  enddo
 enddo
enddo
    
gathbuffer_size = imax * js * ks * nvarout
call MPI_GATHER(send_buff, gathbuffer_size, MPI_REAL, recv_buff, gathbuffer_size, MPI_REAL, 0, MPI_COMM_ROW, mpierr)

!--------------------------------------------------------------------------------
! only the mypey=0 procs create nc files, unload receive buffer and write data to disk

if (mypey == 0) then

! Create filename from integers nfile and mypez and prefix such that filename
! looks like prefx_1000.0000.nc where 1000 is the value of nfile and 0000 is mypez
! For 2D create filename from integer nfile and prefix such that filename looks like prefx1000.nc

write(tmp1,"(i4)") nfile
nfile = nfile + 1

if (ndim==2) then
 filename = 'output/' // trim(prefix) // tmp1 // '.nc'
else
 write(tmp2,"(i4)") mypez
 do i = 1, 4
  if ((tmp2(i:i)) == ' ') tmp2(i:i) = '0'
 enddo
 filename = 'output/' // trim(prefix) // '_' // tmp1 // '.' // tmp2 // '.nc'
endif

ncstat = nf90_create(filename, nf90_Clobber, ncid)
   if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

if (mype==0) write(8,6001) trim(prefix)//tmp1//'.nc', time, ncycle

! Initialize Dimensions
ncstat = nf90_def_dim(ncid, "x", imax, xDimID); DimIDs(1) = xDimID
   if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)
ncstat = nf90_def_var(ncid, "x", nf90_float,(/ xDimID /), XScale_varID)
   if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

ncstat = nf90_def_dim(ncid, "y", jmax, yDimID); DimIDs(2) = yDimID
   if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)
ncstat = nf90_def_var(ncid, "y", nf90_float,(/ yDimID /), YScale_varID)
   if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

if (ndim==3) then
 ncstat = nf90_def_dim(ncid, "z", ks  , zDimID); DimIDs(ndim) = zDimID
 ncstat = nf90_def_var(ncid, "z", nf90_float,(/ zDimID /), ZScale_varID)
 nvars = nvarout
else
 nvars = nvarout - 1  ! drop zuz output if only two dimensions
endif

! define the simulation variables
do nv = 1, nvars
 ncstat = nf90_def_var(ncid, varname(nv), nf90_float, DimIDs, varID(nv))
   if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)
enddo

! Label velocity components as vectors
ncstat = nf90_put_att(ncid, varID(3), "vector", "x")
ncstat = nf90_put_att(ncid, varID(4), "vector", "y")
if (ndim==3) ncstat = nf90_put_att(ncid, varID(5), "vector", "z")

! perhaps add variable attributes like units, min/max...

! define some global attributes
ncstat = nf90_put_att(ncid, nf90_global, "coordinatesystem", ngeomx)
ncstat = nf90_put_att(ncid, nf90_global, "time", time)
ncstat = nf90_put_att(ncid, nf90_global, "gamma", gam)
 
! take dataset out of definition mode
ncstat = NF90_ENDDEF(ncid)
   if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

do nv = 1, nvars

  do m = 1, npey
   do k = 1, ks
    do j = 1, js
     jsk = (m-1)*js + j
     do i = 1, imax
       prin_buff(i,jsk,k) = recv_buff(i,nv,j,k,m)
     enddo
    enddo
   enddo
  enddo
  ncstat = nf90_put_var(ncid, varID(nv), prin_buff)
   if (ncstat /= nf90_NoErr ) print *, NF90_STRERROR(ncstat)

enddo

ncstat = nf90_put_var(ncid, XScale_varID, zxc)
ncstat = nf90_put_var(ncid, YScale_varID, zyc)

if (ndim==3) then
 do k = 1, ks
  zshort(k) = zzc(mypez*ks+k)
 enddo
 ncstat = nf90_put_var(ncid, ZScale_varID, zshort)
endif

! always a necessary statment to flush output 
ncstat = nf90_close(ncid)

endif

6001 format('Wrote files for ',a12,' to disk at time =',1pe12.5,' (ncycle =', i8,')')

return
end



