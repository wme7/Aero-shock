subroutine riemann ( lmin, lmax, gamma, prgh, urgh, vrgh, plft, ulft, vlft, pmid, umid )

! Solve the Riemann shock tube problem for the left and right input states,
! using the Newton interation procedure described in van Leer (1979).
!---------------------------------------------------------------------------------
! GLOBALS
use sweepsize

IMPLICIT NONE

! LOCALS
INTEGER :: l, lmin, lmax, n
REAL :: gamma, gamfac1, gamfac2
REAL, DIMENSION(maxsweep) :: plft, prgh, ulft, urgh, vlft, vrgh, wlft, wrgh, clft, crgh, zlft, zrgh
REAL, DIMENSION(maxsweep) ::  umidl, umidr, pmid, umid, pmold, plfti, prghi

REAL, PARAMETER :: smallp = 1.0e-25
REAL, PARAMETER :: tol = 1.0e-5

! Input variables are:
!
!   lmin   = zone number of first physical zone
!   lmax   = zone number of first ghost zone on right (lmax=nmax+1)
!   gamma    = equation of state gamma
!   prgh   = pressure state on the right side of the boundary
!   plft   = pressure state on the left side of the boundary
!   urgh   = velocity state on the right side of the boundary
!   ulft   = velocity state on the left side of the boundary
!   vrgh   = density state on the right side of the boundary
!   vlft   = density state on the left side of the boundary
!             (vlft and vrgh are inverted to get the specific volume)
!
! Output variables are:
! 
!   umid   = time averaged velocity at the zone interface
!   pmid   = time averaged pressure at the zone interface
!
!----------------------------------------------------------------------
gamfac2 = gamma + 1.0
gamfac1 = 0.5*(gamfac2)/gamma

! Obtain first guess for Pmid by assuming Wlft, Wrgh = Clft, Crgh

do l = lmin, lmax
  clft(l) = sqrt(gamma*plft(l)*vlft(l))
  crgh(l) = sqrt(gamma*prgh(l)*vrgh(l))
  vlft(l) = 1.0/vlft(l)
  vrgh(l) = 1.0/vrgh(l)
  plfti(l)= 1.0/plft(l)
  prghi(l)= 1.0/prgh(l)
  pmid(l) = prgh(l) - plft(l) - crgh(l)*(urgh(l)-ulft(l))     
  pmid(l) = plft(l) + pmid(l) * clft(l)/(clft(l)+crgh(l))      
  pmid(l) = max(smallp,pmid(l)) 
enddo
 
! Iterate up to 8 times using Newton's method to converge on correct Pmid
!     -use previous guess for pmid to get wavespeeds: wlft, wrgh
!     -find the slope in the u-P plane for each state: zlft, zrgh
!     -use the wavespeeds and pmid to guess umid on each side: umidl, umidr
!     -project tangents from (pmid,umidl) and (pmid,umidr) to get new pmid
!     -make sure pmid does not fall below floor value for pressure

do l = lmin, lmax
  do n = 1, 12
    pmold(l) = pmid(l)
    wlft (l) = 1.0 + gamfac1*(pmid(l) - plft(l)) * plfti(l)   
    wrgh (l) = 1.0 + gamfac1*(pmid(l) - prgh(l)) * prghi(l)
    wlft (l) = clft(l) * sqrt(wlft(l))      
    wrgh (l) = crgh(l) * sqrt(wrgh(l))   
    zlft (l) = 4.0 * vlft(l) * wlft(l) * wlft(l)      
    zrgh (l) = 4.0 * vrgh(l) * wrgh(l) * wrgh(l)   
    zlft (l) = -zlft(l) * wlft(l)/(zlft(l) - gamfac2*(pmid(l) - plft(l)))      
    zrgh (l) =  zrgh(l) * wrgh(l)/(zrgh(l) - gamfac2*(pmid(l) - prgh(l)))     
    umidl(l) = ulft(l) - (pmid(l) - plft(l)) / wlft(l)      
    umidr(l) = urgh(l) + (pmid(l) - prgh(l)) / wrgh(l)   
    pmid (l) = pmid(l) + (umidr(l) - umidl(l))*(zlft(l) * zrgh(l)) / (zrgh(l) - zlft(l))   
    pmid (l) = max(smallp,pmid(l))
    if (abs(pmid(l)-pmold(l))/pmid(l) < tol ) exit
  enddo
enddo

! Calculate umid by averaging umidl, umidr based on new pmid
do l = lmin, lmax
  umidl(l) = ulft(l) - (pmid(l) - plft(l)) / wlft(l)      
  umidr(l) = urgh(l) + (pmid(l) - prgh(l)) / wrgh(l)   
  umid (l) = 0.5*(umidl(l) + umidr(l)) 
enddo
 
return     
end 
