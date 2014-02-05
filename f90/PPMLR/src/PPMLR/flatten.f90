subroutine flatten

! Flaten looks for signs of strong shocks and sets the variable flat
! between 0.0 (smooth flow) and 1.0 (strong shock).
! Simplified method of C&W: eqn. A.1 and A.2.
!-----------------------------------------------------------------------

! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: n
REAL, DIMENSION(maxsweep) :: steep
REAL :: delp1, delp2, shock, temp1, temp2, old_flat
                                    ! use the numbers below for more aggresive flattening      
REAL, PARAMETER :: omega1 = 0.75    ! 0.5    
REAL, PARAMETER :: omega2 = 5.0     ! 10.0
REAL, PARAMETER :: epsilon = 0.33   ! 1.00

!--------------------------------------------------------------------------
! Look for presence of a shock using pressure gradient and sign of
! velocity jump:  shock = 1 if there is a shock in the zone, else shock = 0
! Compute steepness parameter based on steepness of pressure jump IF 
! there is a shock.

do n = nmin-4, nmax+4
  delp1 = p(n+1) - p(n-1)
  delp2 = p(n+2) - p(n-2)
  if(abs(delp2) < small) delp2 = small
  shock = abs(delp1)/min(p(n+1),p(n-1))-epsilon
  shock = max(0.0,shock)
  if(shock > 0.0) shock = 1.0
  if(u(n-1) < u(n+1)) shock = 0.0
  temp1 = ( delp1 / delp2 - omega1 ) * omega2
  steep(n) = shock * max( 0., temp1 )
enddo

! Set phony boundary conditions for the steepness parameter

steep(nmin-5) = steep(nmin-4)
steep(nmax+5) = steep(nmax+4)

! Set flatening coefficient based on the steepness in nieghboring zones

flat = 0.0
do n = nmin-4, nmax+4
  temp2   = max( steep(n-1), steep(n), steep(n+1) )
  flat(n) = max( 0.0, min( 0.5, temp2 ) )
enddo

! flat(n) should be set to old_flat if no shock in this direction

do n = nmin-3, nmax+3
  old_flat = f(n) - int(f(n))
  if (flat(n) > 0.0) then
    flat(n) = max(flat(n),old_flat)
    f(n)    = max(flat(n) + (2.0*ndim-3.0), 0.0)   ! f(n)=0 for 1D, 1.f for 2D, 3.f for 3D
   else
    f(n)    = max(0.0, f(n) - 1.0)
    flat(n) = old_flat
  endif
enddo

return
end
