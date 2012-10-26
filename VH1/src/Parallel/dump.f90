subroutine dump(dumpfile)

! Write out a raw binary file containing all variables needed to continue computation
!-------------------------------------------------------------------------------------
! GLOBALS     
use global
use zone

IMPLICIT NONE

! LOCALS
CHARACTER(len=8)  :: dumpfile
CHARACTER(len=50) :: filename
CHARACTER(len=1)  :: sf1,sf2,char
INTEGER :: isf1, isf2

!----------------------------------------------------------------------

filename = 'output/' // trim(prefix) // dumpfile
open(unit=4,file=filename,status='unknown',form='unformatted')
write(4) zro,zpr,zux,zuy,zuz,zfl,zxa,zxc,zdx,zya,zyc,zdy,zza,zzc,zdz,           &
  gam,gamm,time,dt,timem,timep,uinflo,vinflo,winflo,dinflo,pinflo,einflo,       &
  uotflo,votflo,wotflo,dotflo,potflo,eotflo,                                    &
  ncycle,ncycp,ncycm,ncycd,ngeomx,ngeomy,ngeomz,nleftx,nlefty,nleftz,nrightx,   &
  nrighty,nrightz,nfile
close(4)

if (mype == 0) write(8,*) 'Dumped ', trim(prefix) // dumpfile(1:3),' at cycle ', ncycle

! Increment file name by one letter in the suffix
if (dumpfile(3:3) == 'z' .or. dumpfile(3:3) == 'Z') then
  sf1  = dumpfile(2:2)
  isf1 = ichar(sf1)
  sf2  = dumpfile(3:3)
  isf2 = ichar(sf2)
  isf1 = isf1 + 1
  isf2 = isf2 - 25
  dumpfile(2:2) = char(isf1)
  dumpfile(3:3) = char(isf2)
else
  sf2  = dumpfile(3:3)
  isf2 = ichar(sf2)
  isf2 = isf2 + 1
  dumpfile(3:3) = char(isf2)
endif

return
end




subroutine undump(dumpfile)

! Write out a raw binary file containing all variables needed to continue computation
!-------------------------------------------------------------------------------------
! GLOBALS     
use global
use zone

IMPLICIT NONE

! LOCALS
CHARACTER(len=8)  :: dumpfile, todayis
CHARACTER(len=50) :: filename
CHARACTER(len=1)  :: sf1,sf2,char
INTEGER :: isf1, isf2

!----------------------------------------------------------------------

if (mype == 0) then
  call date_and_time(todayis)
  write(8,*) '-------------------------------------------------'
  write(8,*) 'Restarting from ',trim(prefix)//dumpfile(1:3), ' on ',todayis(5:6),' / ',todayis(7:8),' / ',todayis(1:4)
endif

filename = 'output/' // trim(prefix) // dumpfile
open(unit=4,file=filename,status='old',form='unformatted')
read(4) zro,zpr,zux,zuy,zuz,zfl,zxa,zxc,zdx,zya,zyc,zdy,zza,zzc,zdz,            &
  gam,gamm,time,dt,timem,timep,uinflo,vinflo,winflo,dinflo,pinflo,einflo,       &
  uotflo,votflo,wotflo,dotflo,potflo,eotflo,                                    &
  ncycle,ncycp,ncycm,ncycd,ngeomx,ngeomy,ngeomz,nleftx,nlefty,nleftz,           &
  nrightx,nrighty,nrightz,nfile
close(4)
  
if(dumpfile(3:3) == 'z' .or. dumpfile(3:3) == 'Z') then ! Increment dump filename
   sf1  = dumpfile(2:2)
   isf1 = ichar(sf1) + 1
   sf2  = dumpfile(3:3)
   isf2 = ichar(sf2) - 25
   dumpfile(2:2) = char(isf1)
   dumpfile(3:3) = char(isf2)
else 
   sf2  = dumpfile(3:3)
   isf2 = ichar(sf2) + 1
   dumpfile(3:3) = char(isf2)
endif

return
end


