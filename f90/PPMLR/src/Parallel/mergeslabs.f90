program mergevh1data

! read in MVH-1 data in netCDF format, sliced in j; write out concatenated netCDF datasets
! compile with:
! ifort -o mergeslabs mergeslabs.f90 -I/share/apps/netcdf/include -L/share/apps/netcdf/lib -lnetcdf

      
use netcdf ! netcdf modules
      
IMPLICIT NONE

LOGICAL, DIMENSION(10) :: isvector

CHARACTER(LEN=4) :: tmp, tmp1, tmp2
CHARACTER(LEN=50) :: base
CHARACTER(LEN=50) :: filename
CHARACTER(LEN=50) :: masterfile
CHARACTER(LEN=16) :: coord1, coord2, coord3
CHARACTER(LEN=16), DIMENSION(9) :: varname, vector

INTEGER, DIMENSION(9) :: varID
INTEGER :: xDimID, yDimID, zDimID, x_varID, y_varID, z_varID
INTEGER :: status, ncid, ncid_slab, ncstat

INTEGER, DIMENSION(3) :: dims, start, mdims
INTEGER :: imax, jmax, kmax, ks, nv, nvar

INTEGER :: n, i, j, k, rank, npe, nstart, nstop, nfile, nsteps, error, existing, npey, ngeom

REAL(KIND=4) :: time
REAL(KIND=4), DIMENSION(:,:,:), ALLOCATABLE :: slab
REAL(KIND=4), DIMENSION(:), ALLOCATABLE :: zxc, zyc, zzc, zzs

!######################################################################################################
open(27, file="postprocess", status="OLD", iostat=existing)
if (existing > 0) then
  write(*,*) 'Input name of file to paste'
  read(5,*) base
  write(*,*) 'Input first file number'
  read(5,*) nstart
  write(*,*) 'Input number of time steps to merge'
  read(5,*) nsteps
  write(*,*) 'Input number of mpi tasks used in simulation'
  read(5,*) npe
else
  read (27,*) base
  read (27,*) nstart
  read (27,*) nsteps
  read (27,*) npe
  close(27)
endif

n = 0
nstop = nstart + nsteps - 1

! create a filename for the first slab of the first data set
write(tmp1,"(i4)") nstart
write(tmp2,"(i4)") n
do i = 1, 4
  if ((tmp1(i:i)) .eq. ' ') tmp1(i:i) = '0'
  if ((tmp2(i:i)) .eq. ' ') tmp2(i:i) = '0'
enddo
filename = trim(base) // '_' // tmp1 // '.' // tmp2 // '.nc'

! open up first data file to retrieve dimensions and coordinates
ncstat = nf90_open(filename, nf90_nowrite, ncid_slab)
ncstat = nf90_inquire(ncid_slab, rank, nvar)
ncstat = nf90_inquire_dimension(ncid_slab, 1, coord1, imax)
ncstat = nf90_inquire_dimension(ncid_slab, 2, coord2, jmax)
ncstat = nf90_inquire_dimension(ncid_slab, 3, coord3, ks)

nvar = nvar - rank  ! first 'rank' variables are assumed to be the coordinate arrays
do nv = 1, nvar
  ncstat = nf90_inquire_variable(ncid_slab, nv+rank, varname(nv))
  write(*,*) nv, varname(nv), ncstat
  ncstat = nf90_inquire_attribute(ncid_slab, nv+rank, "vector")
  if (ncstat==nf90_noerr) then
    isvector(nv) = .true.
    ncstat = nf90_get_att(ncid_slab,nv+rank,"vector",vector(nv))
  else
    isvector(nv) = .false.
  endif
enddo
ncstat = nf90_get_att(ncid_slab, nf90_global, 'time', time)
ncstat = nf90_get_att(ncid_slab, nf90_global, 'coordinatesystem', ngeom)
write(*,*) 'time = ', time
ncstat = nf90_close(ncid_slab)

dims(1) = imax
dims(2) = jmax
dims(3) = ks
kmax = ks * npe

! allocate data array
ALLOCATE( slab(dims(1),dims(2),dims(3)) )
ALLOCATE( zxc(dims(1)) )
ALLOCATE( zyc(dims(2)) )
ALLOCATE( zzc(dims(2)) )
ALLOCATE( zzs(dims(3)) )

! mdims is the shape of the full dataset
mdims(1) = imax
mdims(2) = jmax
mdims(3) = kmax

! start is the pointer to the start of the hyperslab in the full array
start(1) = 1
start(2) = 1
start(3) = 1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

do nfile = nstart, nstop

 n = 0
 write(tmp1,"(i4)") nfile
 write(tmp2,"(i4)") n
 do i = 1, 4
   if ((tmp1(i:i)) .eq. ' ') tmp1(i:i) = '0'
   if ((tmp2(i:i)) .eq. ' ') tmp2(i:i) = '0'
 enddo
 filename   = trim(base) // '_' // tmp1 // '.' // tmp2 // '.nc'
 masterfile = trim(base) // tmp1 // '.nc'
 

 ! open up master data file
 status = nf90_create(masterfile, nf90_Clobber, ncid)
 
 ! define the dimensions
status = nf90_def_dim(ncid, coord1, mdims(1), xDimID)
status = nf90_def_dim(ncid, coord2, mdims(2), yDimID)
status = nf90_def_dim(ncid, coord3, mdims(3), zDimID)

 ! define the coordinate variables
status = nf90_def_var(ncid, coord1, nf90_float, (/ xDimID /), x_varID)
status = nf90_def_var(ncid, coord2, nf90_float, (/ yDimID /), y_varID)
status = nf90_def_var(ncid, coord3, nf90_float, (/ zDimID /), z_varID)

 ! define the simulation variables
do nv = 1, nvar
 status = nf90_def_var(ncid, varname(nv), nf90_float, (/ xDimID, yDimID, zDimID /), varID(nv) )
 if (isvector(nv)) status = nf90_put_att(ncid, varID(nv), "vector", vector(nv))
enddo
 
 ! initialize any attributes
 status = nf90_put_att(ncid, nf90_global, "time", 0.0)
 status = nf90_put_att(ncid, nf90_global, "coordinatesystem", ngeom)

 ! take dataset out of definition mode
 status = NF90_ENDDEF(ncid)
 write(*,*) status, 'finished definition mode'
 if(status.ne.0) write(*,*) NF90_STRERROR(status)

 !#############################  LOOP OVER SLABS ###########################
 do n = 0, npe-1

   write(tmp2,"(i4)") n
   do i = 1, 4
     if ((tmp2(i:i)) .eq. ' ') tmp2(i:i) = '0'
   enddo
   filename   = trim(base) // '_' // tmp1 // '.' // tmp2 // '.nc'
   write(*,*) filename
   
   ncstat = nf90_open(filename, nf90_nowrite, ncid_slab)
   
     ncstat = nf90_get_att(ncid_slab, nf90_global, 'time', time)
     ncstat = nf90_get_var(ncid_slab, 1, zxc) 
     ncstat = nf90_get_var(ncid_slab, 2, zyc)
     ncstat = nf90_get_var(ncid_slab, 3, zzs)
              do k = 1, ks
                zzc(n*ks+k) = zzs(k)
              enddo
   
     start(3) = ks * n + 1
     do nv = 1, nvar
      ncstat = nf90_get_var(ncid_slab, nv+rank, slab)
      ncstat = nf90_put_var(ncid, varID(nv), slab, start, dims)
     enddo

   ncstat = nf90_close(ncid_slab)
      
 enddo
 !############################################################################

 
 ! write coordinates to master data file
 status = nf90_put_var(ncid, x_varID, zxc)
 status = nf90_put_var(ncid, y_varID, zyc)
 status = nf90_put_var(ncid, z_varID, zzc)
 status = nf90_put_att(ncid, nf90_global, "time", time)

 status = nf90_close(ncid)
 if(status.ne.0) write(*,*) 'File did not close properly: ',NF90_STRERROR(status)
 
enddo
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


END PROGRAM mergevh1data
