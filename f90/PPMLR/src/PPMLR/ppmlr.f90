subroutine ppmlr

! Using the 1D arrays of rho, u, and P, perform a 1D lagrangian hydrodynamics evolution:
!   - set boundary conditions by filling ghost zones
!   - obtain parabolic interpolations of rho, u, P
!   - compute input states from these interpolations for Riemann problem
!   - call the Riemann solver to find the time averages umid, pmid
!   - evolve the finite difference equations to get updated values of rho, u, and E
! Then perform a conservative remap back to the Eulerian grid
!-----------------------------------------------------------------------

! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: n
REAL :: xwag
REAL, DIMENSION(maxsweep) :: dr, du, dp, r6, u6, p6, rl, ul, pl
REAL, DIMENSION(maxsweep) :: rrgh, urgh, prgh, rlft, ulft, plft, umid, pmid
REAL, DIMENSION(maxsweep) :: xaf, dxf

!-----------------------------------------------------------------------
if (ngeom<3) radius = 1.0

! Apply boundary conditions by filling ghost zones
call boundary

! Calculate flattening coefficients for smoothing near shocks
call flatten

! Compute parabolic coefficients 
call paraset( nmin-4, nmax+5, para, dx, xa )

! Interpolate parabolae for fluid variables 
call parabola(nmin-4, nmax+4, para, p, dp, p6, pl, flat)
call parabola(nmin-4, nmax+4, para, r, dr, r6, rl, flat)
call parabola(nmin-4, nmax+4, para, u, du, u6, ul, flat)

! Integrate parabolae to get input states for Riemann problem
call states( pl, ul, rl, p6, u6, r6, dp, du, dr, plft, ulft, rlft, prgh, urgh, rrgh )
   
! Call the Riemann solver to obtain the zone face averages, umid and pmid
call riemann( nmin-3, nmax+4, gam, prgh, urgh, rrgh, plft, ulft, rlft, pmid, umid )

! do lagrangian update using umid and pmid
call evolve( umid, pmid )

!#########################################################################
! EXTRA DISSIPATION TO REDUCE CARBUNCLE NOISE
xwag = sum(flat)
if (xwag*xwig /= 0.0) then ! wiggle grid, remap, then remap back to Eulerian grid

 ! save Eulerian coordinates for second remap
 xaf = xa0
 dxf = dx0

 ! wiggle grid where there is a shock, leaving edges unchanged
 do n = nmin+1, nmax
  if (max(flat(n-1),flat(n)) > 0.0) xa0(n) = xa0(n) + xwig*dx0(n)
  dx0(n-1) = xa0(n) - xa0(n-1)
 enddo
 dx0(nmax) = xa0(nmax+1) - xa0(nmax)

 call remap

 ! put wiggled grid into xa(n), but keep ghost cells untouched
 do n = nmin, nmax+1
  xa(n)   = xa0(n)
  dx(n-1) = xa(n) - xa(n-1)
 enddo
 call volume(nmin, nmax, ngeom, radius, xa, dx, dvol)

 ! put Eulerian grid back into xa0
 xa0 = xaf
 dx0 = dxf

endif
!#########################################################################

! remap onto original Eulerian grid
call remap

return
end

