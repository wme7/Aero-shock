subroutine boundary

! Impose boundary conditions on ends of 1D sweep arrays
!-----------------------------------------------------------------------

! GLOBALS
use global
use sweeps

IMPLICIT NONE

! LOCALS
INTEGER :: n

!-----------------------------------------------------------------------
!  Boundary condition flags : nleft, nright
!    = 0 : reflecting
!    = 1 : outflow (zero gradients)
!    = 2 : fixed inflow (eg, uinflo,pinflo,...)
!    = 3 : periodic (eg, u(nmin-1) = u(nmax))

if ( nleft == 0 ) then
  do n = 1, 6
    dx (nmin-n)= dx (nmin+n-1)
    dx0(nmin-n)= dx0(nmin+n-1)
    xa (nmin-n)= xa (nmin-n+1) - dx (nmin-n)
    xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)
    r (nmin-n) = r (nmin+n-1)
    u (nmin-n) = -u(nmin+n-1)
    v (nmin-n) = v (nmin+n-1)
    w (nmin-n) = w (nmin+n-1)
    p (nmin-n) = p (nmin+n-1)
    e (nmin-n) = e (nmin+n-1)
    f (nmin-n) = f (nmin+n-1)
  enddo
else if ( nleft == 1 ) then
  do n = 1, 6
    dx (nmin-n)= dx (nmin)
    dx0(nmin-n)= dx0(nmin)
    xa (nmin-n)= xa (nmin-n+1) - dx (nmin-n)
    xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)
    r (nmin-n) = r (nmin)
    u (nmin-n) = u (nmin)
    v (nmin-n) = v (nmin)
    w (nmin-n) = w (nmin)
    p (nmin-n) = p (nmin)
    e (nmin-n) = e (nmin)
    f (nmin-n) = f (nmin)
  enddo
else if ( nleft == 2 ) then
  do n = 1, 6
    dx (nmin-n)= dx (nmin)
    dx0(nmin-n)= dx0(nmin)
    xa (nmin-n)= xa (nmin-n+1) - dx (nmin-n)
    xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)
    r (nmin-n) = dinflo
    u (nmin-n) = uinflo
    v (nmin-n) = vinflo
    w (nmin-n) = winflo
    p (nmin-n) = pinflo
    e (nmin-n) = pinflo/(dinflo*gamm) + 0.5*( uinflo**2 + vinflo**2 + winflo**2 )
    f (nmin-n) = 0.0
  enddo
else if ( nleft == 3 ) then
  do n = 1, 6
    dx (nmin-n)= dx (nmax+1-n)
    dx0(nmin-n)= dx0(nmax+1-n)
    xa (nmin-n)= xa (nmin-n+1) - dx (nmin-n)
    xa0(nmin-n)= xa0(nmin-n+1) - dx0(nmin-n)
    r (nmin-n) = r (nmax+1-n)
    u (nmin-n) = u (nmax+1-n)
    v (nmin-n) = v (nmax+1-n)
    w (nmin-n) = w (nmax+1-n)
    p (nmin-n) = p (nmax+1-n)
    e (nmin-n) = e (nmax+1-n)
    f (nmin-n) = f (nmax+1-n)
  enddo
endif

if ( nright == 0 ) then     
  do n = 1, 6
    dx (nmax+n)= dx (nmax+1-n)
    dx0(nmax+n)= dx0(nmax+1-n)
    xa (nmax+n)= xa (nmax+n-1) + dx (nmax+n-1)
    xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)
    r (nmax+n) = r (nmax+1-n)
    u (nmax+n) = -u(nmax+1-n)
    v (nmax+n) = v (nmax+1-n)
    w (nmax+n) = w (nmax+1-n)
    p (nmax+n) = p (nmax+1-n)
    e (nmax+n) = e (nmax+1-n)
    f (nmax+n) = f (nmax+1-n)
  enddo
else if ( nright == 1 ) then     
  do n = 1, 6
    dx (nmax+n)= dx (nmax)
    dx0(nmax+n)= dx0(nmax)
    xa (nmax+n)= xa (nmax+n-1) + dx (nmax+n-1)
    xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)
    r (nmax+n) = r (nmax)
    u (nmax+n) = u (nmax)
    v (nmax+n) = v (nmax)
    w (nmax+n) = w (nmax)
    p (nmax+n) = p (nmax)
    e (nmax+n) = e (nmax)
    f (nmax+n) = f (nmax)
  enddo
else if ( nright == 2 ) then     
  do n = 1, 6
    dx (nmax+n)= dx (nmax)
    dx0(nmax+n)= dx0(nmax)
    xa (nmax+n)= xa (nmax+n-1) + dx (nmax+n-1)
    xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)
    r (nmax+n) = dotflo
    u (nmax+n) = uotflo
    v (nmax+n) = votflo
    w (nmax+n) = wotflo
    p (nmax+n) = potflo
    e (nmax+n) = potflo/(dotflo*gamm) + 0.5*( uotflo**2 + votflo**2 + wotflo**2 )
    f (nmax+n) = 0.0
  enddo
else if ( nright == 3 ) then     
  do n = 1, 6
    dx (nmax+n)= dx (nmin+n-1)
    dx0(nmax+n)= dx0(nmin+n-1)
    xa (nmax+n)= xa (nmax+n-1) + dx (nmax+n-1)
    xa0(nmax+n)= xa0(nmax+n-1) + dx0(nmax+n-1)
    r (nmax+n) = r (nmin+n-1)
    u (nmax+n) = u (nmin+n-1)
    v (nmax+n) = v (nmin+n-1)
    w (nmax+n) = w (nmin+n-1)
    p (nmax+n) = p (nmin+n-1)
    e (nmax+n) = e (nmin+n-1)
    f (nmax+n) = f (nmin+n-1)
  enddo
endif

return
end

