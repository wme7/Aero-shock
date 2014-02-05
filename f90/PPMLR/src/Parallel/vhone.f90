program vhone

!--------------------------------------------------------------------------
!
!    MM        MM    VV              VV   HH       HH       3333333
!    MMMM    MMMM     VV            VV    HH       HH      333    333
!    MM MM  MM MM      VV          VV     HH       HH              333
!    MM  MMMM  MM       VV        VV      HH       HH            333
!    MM   MM   MM        VV      VV       HHHHHHHHHHH          333
!    MM        MM         VV    VV        HH       HH            333    
!    MM        MM          VV  VV         HH       HH              333
!    MM        MM           VVVV          HH       HH      333    333
!    MM        MM            VV           HH       HH       3333333
!
! 
!                   MPI - VIRGINIA HYDRODYNAMICS - 3D
!
!--------------------------------------------------------------------------
!  The Virginia Numerical Bull Session ideal hydrodynamics PPMLR
!  
!  Version 3.1 July 2008
!    Cleaned up code, minor changes to algorithm, 2/3D output
!    Diced MPI version uses a 'julienne' domain decomposition for better scaling
!  Version 3.0 c. 2003
!    Converted to F90 free-format, use netcdf for binary output
!  Version 2.0 July 1997
!    Many small modifications, mainly brought BCs into only one call to sweepbc 
!  Version 1.0 October 1991 
!    This version was removed from the ``loop'' in the continued development
!    of VH-1 on Sep 1 1990, but has followed a similar evolution until this
!    release date.
! 
! Output channels:
!                 4 = aaaaadaa.####,   binary dump files for restarting
!                 8 = aaaaa.hst,   history file
!                     aaaaaXY.nc,   netCDF movie of XY slice
!                     aaaaa_1000.NNN, netCDF output of dataset (all variables) for slab NNN
!-------------------------------------------------------------------------------

! GLOBALS
use global
use zone
use sweepsize
use mpi

IMPLICIT NONE

! LOCALS
CHARACTER(len=50) :: hstfile
CHARACTER(len=8)  :: dmpfile
CHARACTER(len=4)  :: tmp4 
CHARACTER(len=2)  :: rstrt 
CHARACTER(len=8)  :: todayis
CHARACTER(len=MPI_MAX_PROCESSOR_NAME) :: myname

INTEGER :: jcol, krow, n, nfile_start
INTEGER :: ncycend, nprin, ndump, nmovie, tmpcyc, mpierr,start_cycl, namelength
REAL :: endtime,tprin,tmovie, olddt
REAL(8) :: start_time,end_time,run_time,zones

NAMELIST / hinput / rstrt, prefix, ncycend, ndump, nprin, nmovie, endtime, tprin, tmovie

!-----------------------------------------------------------------------------------------
! Initialize MPI
call MPI_INIT(mpierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,  npe, mpierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, mype, mpierr)

! Set VH1_DATATYPE for MPI calls
if (SIZEOF(endtime)==8) then
  VH1_DATATYPE = MPI_DOUBLE_PRECISION
elseif (SIZEOF(endtime)==4) then
  VH1_DATATYPE = MPI_REAL
else
  if(mype==0) write(*,*) 'Mode not recognized. MPI commands will fail.'
  stop
endif

! set dimensions of local simulation grid  (imax X js X ks)
js   = jmax / pey
ks   = kmax / pez
isy  = imax / pey
isz  = imax / pez
ndim = 3; if (kmax==1) ndim = 2

! set size of send/receive buffers for alltoall calls
Za2abuff_size = ks * isz * js * nvar
Ya2abuff_size = ks * isy * js * nvar

jcol = mod(mype,pey)   ! processor rank in the Y direction
krow = mype / pey      ! processor rank in the Z direction

! Create communicators for each value of krow (for XY transpose)
call MPI_COMM_SPLIT(MPI_COMM_WORLD, krow, mype, MPI_COMM_ROW, mpierr)
  call MPI_COMM_RANK(MPI_COMM_ROW, mypey, mpierr)
  call MPI_COMM_SIZE(MPI_COMM_ROW, npey, mpierr)

! Create communicators for each value of jcol (for XZ transpose)
call MPI_COMM_SPLIT(MPI_COMM_WORLD, jcol, mype, MPI_COMM_COL, mpierr)
  call MPI_COMM_SIZE(MPI_COMM_COL, npez, mpierr)
  call MPI_COMM_RANK(MPI_COMM_COL, mypez, mpierr)

! should have npey = pey, npez = pez, mypey = jcol, mypez = krow

! Check that code is running with desired number of PE's
if (pey*pez /= npe) then
  write(*,*) 'sorry, I was compiled for ',pey*pez,' PEs, but npe= ', npe
  stop
endif

! Check that arrays are large enough for desired number of physical zones
if (max(imax,jmax,kmax)+12 > maxsweep) then
  write(*,*) 'maxsweep too small'
  stop
endif

! Check that imax and jmax are multiples of npey
if (isy*npey /= imax) then
  write(*,*) 'imax is not an integer multiple of pey', imax, pey
  stop
endif
if (js*npey /= jmax) then
  write(*,*) 'jmax is not an integer multiple of pey', jmax, pey
  stop
endif

! Check that imax and kmax are multiples of npez
if (isz*npez /= imax) then
  write(*,*) 'imax is not an integer multiple of pez', imax, pez
  stop
endif
if (ks*npez /= kmax) then
  write(*,*) 'kmax is not an integer multiple of pez', kmax, pez
  stop
endif

! Begin by reading input deck, and defining file names

open (unit=15,file='indat',status='old',form='formatted')
read (15,nml=hinput)
close(15)

write(tmp4,"(i4)") mype
do n = 1, 4
  if ((tmp4(n:n)) == ' ') tmp4(n:n) = '0'
enddo
dmpfile = 'daa.' // tmp4

if ((rstrt == 'NO').or.(rstrt == 'no')) then ! Initialize variables for new problem; 

  if (mype == 0) then

   hstfile = 'output/' // trim(prefix) // '.hst'
   open (unit=8,file=hstfile,form='formatted')
   call date_and_time(todayis)
   write (8,*) 'History File for VH1-MPI simulation run on ', todayis(5:6), ' / ', todayis(7:8), ' / ', todayis(1:4)

   call MPI_GET_PROCESSOR_NAME(myname, namelength, mpierr)
   write (8,"(' machine id: ',a25)") myname
   if (ndim==3) then
    write (8,"(' Running 3 dimensions on ',i4,' processors, split ',i3, ' X ', i3)") npe,npey,npez
   else
    write (8,"(' Running 2 dimensions on ',i3,' processors' )") npe
   endif

   if(VH1_DATATYPE==MPI_DOUBLE_PRECISION) then
     write(8,*) 'Running in double precision mode'
   else
     write(8,*) 'Running in single precision mode'
   endif
   write (8,*)

  endif

  call init
  nfile_start = nfile
  call prin

else ! Restart from old dumpfile...
  
  dmpfile = 'd' // rstrt(1:2) // dmpfile(4:8)
  if (mype == 0) then
   hstfile = 'output/' // trim(prefix) // '.hst'
   open (unit=8,file=hstfile,form='formatted',position='append')
   call date_and_time(todayis)
   write (8,*) 'Restarting run on ', todayis(5:6), ' / ', todayis(7:8), ' / ', todayis(1:4)
  endif
  call undump(dmpfile)
  if(timep >= tprin)  timep = tprin -2.0*dt    ! adjust so that dt is not set to zero on first step
  if(timem >= tmovie) timem = tmovie-2.0*dt
  nfile_start = nfile

endif

if(mype == 0) then
  start_time = MPI_WTIME()
  start_cycl = ncycle
  write(8,*) 'Starting on cycle number ',ncycle
endif

!############################################################################
!                         MAIN COMPUTATIONAL LOOP

do while (ncycle < ncycend)

!  if(mype == 0) write(*,*) 'STEP = ',ncycle

  ncycle = ncycle + 2
  ncycp  = ncycp  + 2
  ncycd  = ncycd  + 2
  ncycm  = ncycm  + 2
  olddt  = dt
  svel   = 0.

  if ( time + 2*dt > endtime ) then ! set dt to land on endtime
    if(mype==0) write(8,*) 'cutting to the end...', ncycle, ncycend
    dt = 0.5*(endtime - time)
    ncycend = ncycle-1
    ncycp   = nprin
    ncycd   = ndump
  else if ( timep+2.0*dt > tprin ) then ! set dt to land on tprin
    dt = 0.5*(tprin - timep)
    ncycp = nprin
  else if ( timem+2*dt > tmovie ) then ! set dt to land on tmovie
    dt = 0.5*(tmovie - timem)
    ncycm = nmovie
  endif

! Alternate sweeps to approximate 2nd order operator splitting

  call sweepx1
  call sweepy

  if (ndim==3) then
   call sweepz  
   call sweepy
  endif

  call sweepx2

  time  = time  + 2.0*dt
  timep = timep + 2.0*dt
  timem = timem + 2.0*dt
  dt = olddt

  call dtcon   ! Check constraints on the timestep

  ! output movie images/datasets/dumpfiles based on simulation time or cycle #

  if ( ncycm >= nmovie .or. timem >= tmovie ) then
    if (mypez == 0) call imagesxy
    if ((mypey == 0).and.(ndim==3)) call imagesxz
    timem = 0.0
    ncycm = 0
  endif

  if ( ncycp >= nprin .or. timep >= tprin ) then
    call prin
    timep = 0.
    ncycp = 0
  endif

  if ( ncycd >= ndump ) then
    ncycd = 0
    call dump(dmpfile)
  endif

enddo

!                           END OF MAIN LOOP
!#########################################################################
      
if(mype == 0) then

  end_time = MPI_WTIME()
  zones = (imax*js*ks)*(ncycle - start_cycl)
  run_time = end_time - start_time
  write(8,*) 'successful stop at cycle number', ncycle
  if (run_time>4000.0) then
    write(8,*) 'elapsed time = ',run_time/3600., ' hours'
  else
    write(8,*) 'elapsed time = ',run_time/60., ' minutes'
  endif
  write(8,"('speed = ',f5.1,' kz/s/pe')") 1.0e-3*zones/run_time
  close( 8 )

  if (ndim==3) then ! open up a file to store data needed for post-processing
   open(9,file='output/postprocess')
     write(9,*) trim(prefix)
     write(9,*) nfile_start
     write(9,*) nfile - nfile_start
     write(9,*) npez
   close(9)
  endif

endif
      
call MPI_FINALIZE(mpierr)

stop      
end 
