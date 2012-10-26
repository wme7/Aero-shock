module zone 
!=======================================================================
! (formerly zone.h) global (3D) data arrays
!======================================================================= 
 
 INTEGER, PARAMETER :: imax=100, jmax=100, kmax=1   ! Memory dimensions

 REAL, DIMENSION(imax,jmax,kmax) :: zro, zpr, zux, zuy, zuz, zfl
 
 REAL, DIMENSION(imax) :: zxa, zdx, zxc
 REAL, DIMENSION(jmax) :: zya, zdy, zyc
 REAL, DIMENSION(kmax) :: zza, zdz, zzc
 
 INTEGER :: ngeomx, ngeomy, ngeomz       ! XYZ Geometry flag
 INTEGER :: nleftx, nlefty, nleftz       ! XYZ Lower Boundary Condition
 INTEGER :: nrightx,nrighty,nrightz      ! XYZ Upper Boundary Condition

end module zone
