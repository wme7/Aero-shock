subroutine parabola( nmin, nmax, para, a, deltaa, a6, al, flat )

! Colella and Woodward, JCompPhys 54, 174-201 (1984) eq 1.5-1.8,1.10
!
! parabola calculates the parabolas themselves. call paraset first
! for a given grid-spacing to set the constants, which can be reused
! each time parabola is called.
!
! flatening coefficients are calculated externally in flaten. 
! nmin/nmax are indicies over which the parabolae are calculated
!-----------------------------------------------------------------------
! GLOBALS
use sweepsize

IMPLICIT NONE

! LOCALS
integer :: n, nmin, nmax
real :: onemfl
real, dimension(maxsweep) :: a, deltaa, a6, al, flat, da, ar, diffa, scrch1, scrch2, scrch3
real, dimension(maxsweep,5) :: para

!----------------------------------------------------------------------
  do n = nmin-2, nmax+1
    diffa(n) = a(n+1) - a(n)
  enddo

!                                                       Equation 1.7 of C&W
!     da(j) = D1 * (a(j+1) - a(j)) + D2 * (a(j) - a(j-1))
  do n = nmin-1, nmax+1
   da(n) = para(n,4) * diffa(n) + para(n,5) * diffa(n-1)
   da(n) = sign( min(abs(da(n)), 2.0*abs(diffa(n-1)), 2.0*abs(diffa(n))), da(n) )
  enddo

!     zero out da(n) if a(n) is a local max/min
  do n = nmin-1, nmax+1
    if(diffa(n-1)*diffa(n) < 0.0) da(n) = 0.0
  enddo

!                                                       Equation 1.6 of C&W
!     a(j+.5) = a(j) + C1 * (a(j+1)-a(j)) + C2 * dma(j+1) + C3 * dma(j)
! MONOT: Limit ar(n) to the range defined by a(n) and a(n+1)

  do n = nmin-1, nmax
    ar(n) = a(n) + para(n,1)*diffa(n) + para(n,2)*da(n+1) + para(n,3)*da(n)
!    ar(n) = max(ar(n),min(a(n),a(n+1)))
!    ar(n) = min(ar(n),max(a(n),a(n+1)))
    al(n+1) = ar(n)
  enddo

! eqn. 4.1 - flaten interpolation in zones with a shock ( flat(n)->1. )

  do n = nmin, nmax
    onemfl= 1.0 - flat(n)
    ar(n) = flat(n) * a(n) + onemfl * ar(n)
    al(n) = flat(n) * a(n) + onemfl * al(n)
  enddo

! MONOTONICITY constraints:

! compute delta_a, a_6
! MONOT: if a is a local max/min, flaten zone structure ar,al -> a.
! MONOT: compute monotonized values using eq. 1.10 of C&W
!        if parabola exceeds al/ar, reset ar/al so that slope -> 0.
! Recalculate delta_a and a_6

do n = nmin, nmax
  deltaa(n) = ar(n) - al(n)
  a6(n)     = 6. * (a(n) - .5 * (al(n) + ar(n)))
  scrch1(n) = (ar(n) - a(n)) * (a(n)-al(n)) 
  scrch2(n) = deltaa(n) * deltaa(n)
  scrch3(n) = deltaa(n) * a6(n)
enddo

do n = nmin, nmax
  if(scrch1(n) <= 0.0) then
    ar(n) = a(n)
    al(n) = a(n)
  endif
  if(scrch2(n) < +scrch3(n)) al(n) = 3. * a(n) - 2. * ar(n)       
  if(scrch2(n) < -scrch3(n)) ar(n) = 3. * a(n) - 2. * al(n)
enddo

do n = nmin, nmax
  deltaa(n)= ar(n) - al(n)
  a6(n) = 6. * (a(n) - .5 * (al(n) + ar(n)))
enddo

return
end

!#######################################################################


subroutine paraset( nmin, nmax, para, dx, xa )

! Colella and Woodward, JCompPhys 54, 174-201 (1984) eq 1.6, 1.7
!
! paraset sets up constants which are re-used each time we want to
! interpolate a parabola on a quantity. First pull out constants
! A, B, and C, and then compute all the equations in terms of those quantities. 
! The quantities calculated here are stored in a array para, to be read by parabola()
! nmin/nmax are index range for which one will calculate parabolae
!-----------------------------------------------------------------------------------
! GLOBALS
use sweepsize

IMPLICIT NONE

! LOCALS
integer :: n, nmin, nmax
real, dimension(maxsweep,5) :: para
real, dimension(maxsweep) :: dx, xa, a, b, c, d, ai, bi, ci

!------------------------------------------------------------------------------

do n = nmin-2, nmax+1
  a (n) = dx(n) + dx(n+1)
  ai(n) = 1.0/a(n)
  b (n) = a(n) + dx(n)
  bi(n) = 1.0/b(n)
  c (n) = a(n) + dx(n+1)
  ci(n) = 1.0/c(n)
enddo

!                                        constants for equation 1.6
!     a(j+.5) = a(j) + C1 * (a(j+1)-a(j)) + C2 * da(j+1) + C3 * da(j)

do n = nmin-1, nmax
  d(n)      = 1. / (a(n-1) + a(n+1))
  para(n,1) = dx(n) * ai(n) + 2. * dx(n+1) * dx(n) * d(n) * ai(n) * ( a(n-1) * bi(n) - a(n+1) * ci(n) )
  para(n,2) = - d(n) * dx(n)   * a(n-1) * bi(n)
  para(n,3) =   d(n) * dx(n+1) * a(n+1) * ci(n)
enddo

!                                        constants for equation 1.7
!     da(j) = D1 * (a(j+1) - a(j)) + D2 * (a(j) - a(j-1))

do n = nmin-1, nmax+1
  d(n) = dx(n) / ( a(n-1) + dx(n+1) )
  para(n,4) = d(n) * b(n-1) * ai(n)
  para(n,5) = d(n) * c(n)   * ai(n-1)
enddo

return
end

