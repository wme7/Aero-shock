subroutine states( pl, ul, rl, p6, u6, r6, dp, du, dr, plft, ulft, rlft, prgh, urgh, rrgh )

! This subroutine takes the values of rho, u, and P at the left hand
! side of the zone, the change accross the zone, and the parabolic 
! coefficients, p6, u6, and rho6, and computes the left and right states
! (integrated over the charachteristics) for each variable for input 
! to the Riemann solver.
!-----------------------------------------------------------------------

! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: n, np
REAL, DIMENSION(maxsweep) :: plft, ulft, rlft, prgh, urgh, rrgh, dp, du, dr
REAL, DIMENSION(maxsweep) :: pl, ul, rl, p6, u6, r6, grav, fict, Cdtdx, fCdtdx
REAL :: hdt

REAL, PARAMETER :: fourthd = 4.0 / 3.0

! Input variables are:
!    pl    =  pressure at left edge of zone
!    dp    =  pressure slope
!    p6    =  pressure parabola coeffecient
!    ul    =  velocity at left edge of zone
!    du    =  velocity slope
!    u6    =  velocity parabola coeffecient
!    rl    =  density at left edge of zone
!    dr    =  density slope
!    r6    =  density parabola coeffecient
!
! Output variables are:
!    plft =  average pressure in the + wave state (left of interface)
!    prgh =  average pressure in the - wave state (right of interface)
!    ulft =  average velocity in the + wave state (left of interface)
!    urgh =  average velocity in the - wave state (right of interface)
!    rlft =  average density in the  + wave state (left of interface)
!    rrgh =  average density in the  - wave state (right of interface)
!
!--------------------------------------------------------------------------
! Calculate the domain of dependence along space coordinate,
! C*dt, by multiplying the Eulerian sound speed in 
! each zone by dt.  Divide by two to save flops later.
! If angular coordinate, then divide out by radius to get physical dimensions

hdt   = 0.5*dt
do n = nmin-4, nmax+4
  Cdtdx (n) = sqrt(gam*p(n)/r(n))/(dx(n)*radius)
  svel      = max(svel,Cdtdx(n))
  Cdtdx (n) = Cdtdx(n)*hdt
  fCdtdx(n) = 1. - fourthd*Cdtdx(n)
enddo

! Include gravitational and ficticious forces
call forces(xa,grav,fict)

! Obtain averages of rho, u, and P over the domain (+/-)Cdt
!                    lft is the + wave on the left  side of the boundary
!                    rgh is the - wave on the right side of the boundary

do n = nmin-4, nmax+4
  np = n+1
  plft(np) = pl(n)+dp(n)-Cdtdx(n)*(dp(n)-fCdtdx(n)*p6(n))
  ulft(np) = ul(n)+du(n)-Cdtdx(n)*(du(n)-fCdtdx(n)*u6(n))
  rlft(np) = rl(n)+dr(n)-Cdtdx(n)*(dr(n)-fCdtdx(n)*r6(n))
  plft(np) = max(smallp,plft(np))
  rlft(np) = max(smallr,rlft(np))
  ulft(np) = ulft(np) + hdt*(grav(np)+fict(np))

  prgh(n) = pl(n) + Cdtdx(n)*(dp(n)+fCdtdx(n)*p6(n))
  urgh(n) = ul(n) + Cdtdx(n)*(du(n)+fCdtdx(n)*u6(n))
  rrgh(n) = rl(n) + Cdtdx(n)*(dr(n)+fCdtdx(n)*r6(n))
  prgh(n) = max(smallp,prgh(n))
  rrgh(n) = max(smallr,rrgh(n))
  urgh(n) = urgh(n) + hdt*(grav(n)+fict(n))
enddo

return
end
