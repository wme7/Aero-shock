module global
!=======================================================================
! global variables accessible by (most) anything
!-----------------------------------------------------------------------

 character(len=50) :: prefix              ! prefix for output filenames

 integer :: ncycle, ncycp, ncycm, ncycd  ! cycle number
 integer :: nfile                        ! output file numbers
 integer :: ndim

 real :: time, dt, timem, timep, svel 
 real :: gam, gamm

 real, parameter :: courant = 0.5           ! timestep fraction of courant limit
 real, parameter :: pi = 3.1415926535897931 ! shouldn't computers know this?
 real, parameter :: xwig = 0.00             ! fraction of a zone to wiggle grid for dissipation
 real, parameter :: smallp = 1.0e-15        ! Set small values to prevent divide by zero
 real, parameter :: smallr = 1.0e-15
 real, parameter :: small  = 1.0e-15

 real :: uinflo, dinflo, vinflo, winflo, pinflo, einflo 
 real :: uotflo, dotflo, votflo, wotflo, potflo, eotflo
      
end module global

module sweepsize
!=======================================================================
! Dimension of 1D sweeps.  maxsweep must be as long as the longest of the 
! 3D arrays PLUS the ghost zones: 
! ie, must have maxsweep >= max(imax,jmax,kmax) + 12
!----------------------------------------------------------------------

integer, parameter :: maxsweep=1036 

end module sweepsize

module sweeps      
!=======================================================================
! Data structures used in 1D sweeps, dimensioned maxsweep  (set in sweepsize.mod)
!----------------------------------------------------------------------

use sweepsize

character(len=1) :: sweep                                    ! direction of sweep: x,y,z
integer :: nmin, nmax, ngeom, nleft, nright                  ! number of first and last real zone  
real, dimension(maxsweep) :: r, p, e, q, u, v, w             ! fluid variables
real, dimension(maxsweep) :: xa, xa0, dx, dx0, dvol          ! coordinate values
real, dimension(maxsweep) :: f, flat                         ! flattening parameter
real, dimension(maxsweep,5) :: para                          ! parabolic interpolation coefficients
real :: radius, theta, stheta

end module sweeps

