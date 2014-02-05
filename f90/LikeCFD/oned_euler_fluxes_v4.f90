
!*****************************************************************************
!* A collection of one-dimensional Euler numerical fluxes, Version 4 (2010),
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!*------------------------
!* List of Flux Functions:
!
!   Lax-Friedrichs
!   Richtmyer
!   MacCormack
!   Steger-Warming
!   Van Leer
!   AUSM
!   Zha-Bilgen
!   Godunov
!   Osher
!   Roe
!   Rusanov
!   HLL(HLLE)
!   HLLL
!   AUFS
!
!*------------------------
!*
!* These F90 routines were written and made available for download
!* for an educational purpose. Detailed description of each numerical
!* flux can be found in the original paper, or in popular textbooks, or
!* in the (will-be-available) second volume of "I do like CFD".
!*
!* Note that all routines are not efficiently implemented for clarity; you can 
!* improve the efficiency and also you can covert it to double precision 
!* version if you wish.
!*
!* Ver.2 (05-10-2010): Rusanov flux debugged.
!* Ver.3 (05-11-2010): Godunov flux debugged.
!* Ver.4 (05-12-2010): AUFS flux added.
!*
!* This file may be updated in future. (to add Hanel, HLLC, LDFSS, etc.)
!*****************************************************************************


!*****************************************************************************
!* --- Lax-Friedrichs's Flux Function ---
!* 
!* P. D. Lax, Weak Solutions of Nonlinear Hyperbolic Equations and Their
!* Numerical Computation, Commun. Pure and Applied Mathematics, 7, 159-193, 
!* 1954.
!*
!* This requires the function, physical_flux(u), included in this file.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function LaxFriedrichs(uL,uR,dt,dx)
 implicit none
 real :: uL(3), uR(3), dt, dx  !  Input (conservative variables rho*[1, v, E])
 real :: LaxFriedrichs(3)      ! Output (numerical flux across L and R states)

 LaxFriedrichs = 0.5*(physical_flux(uR) + physical_flux(uL) - dx/dt*(uR-uL))

 end function LaxFriedrichs
!-----------------------------------------------------------------------------


!*****************************************************************************
!* --- Richtmyer's Flux Function --- (Two-Step Lax-Wendroff Scheme)
!*
!* R. D. Richtmyer and K. W. Morton, `Difference Methods for Initial Value 
!* Problems', second edition, John Wiley & Sons, 1967.
!*
!* This requires the function, physical_flux(u), included in this file.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function Richtmyer(uL,uR,dt,dx)
 implicit none
 real :: uL(3), uR(3), dt, dx  !  Input (conservative variables rho*[1, v, E])
 real :: Richtmyer(3)          ! Output (numerical flux across L and R states)
!Local variables
 real :: us(3)

!Step 1. Compute the intermediate state, us(3).
    us = 0.5 * ( (uL + uR) - dt/dx*( physical_flux(uR) - physical_flux(uL) ) )

!Step 2. Compute the numerical flux, Richtmyer(3).
    Richtmyer = physical_flux(us)

 end function Richtmyer
!-----------------------------------------------------------------------------


!*****************************************************************************
!* --- MacCormack's Flux Function --- (Two-Step Lax-Wendroff Scheme)
!*
!* R. W. MacCormack, The Effect of Viscosity in Hypervelocity Impact Cratering,
!* AIAA Paper 69-354, 1967.
!*
!* This requires the function, physical_flux(u), included in this file.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function MacCormack(uL,uR,dt,dx)
 implicit none
 real :: uL(3), uR(3), dt, dx  !  Input (conservative variables rho*[1, v, E])
 real :: MacCormack(3)         ! Output (numerical flux across L and R states)
!Local variables
 real :: us(3)                 ! Intermediate solution

!Step 1. Compute the intermediate state, us(3).

     us = uL - dt/dx * ( physical_flux(uR) - physical_flux(uL) )

!Step 2. Compute the numerical flux, MacCormack(3).

    MacCormack = 0.5 * ( physical_flux(uR) + physical_flux(us) )

 end function MacCormack
!-----------------------------------------------------------------------------


!*****************************************************************************
!* --- Steger-Warming's Flux-Vector Splitting ---
!*
!* J. L. Steger and R. F. Warming, Flux Vector Splitting of the Inviscid 
!* Gas Dynamics Equations with Application to Finite Difference Methods,
!* Journal of Computational Physics, 40, 263-293, 1981.
!*
!* This requires the function, physical_flux(u), included in this file.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function StegerWarming(uL,uR)
 implicit none
 real :: uL(3), uR(3)     !  Input (conservative variables rho*[1, v, E])
 real :: StegerWarming(3) ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, half, one, two, three  ! Numbers
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: aL, aR, ML, MR               ! Speeds of sound and Mach numbers.
 real :: Fp(3), Fm(3)                 ! F_plus and F_minus
 real :: Z, Q                         ! Temporary variables

!Constants.
     gamma = 1.4
      zero = 0.0
      half = 0.5
       one = 1.0
       two = 2.0
     three = 3.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      ML = vL/aL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      MR = vR/aR

!Positive Part of Flux evaluated in the left cell.
 if     (ML < -one) then
     Fp(1) = zero
     Fp(2) = zero
     Fp(3) = zero
 elseif (ML >  one) then
     Fp    = physical_flux(uL)
 else
         Z = half * rhoL * aL / gamma * (ML + one)
         Q = rhoL * aL * (gamma - one) / gamma
     Fp(1) = Z             + Q       * max(zero, ML)
     Fp(2) = Z*aL*(ML+one) + Q*aL*ML * max(zero, ML)
     Fp(3) = Z*half*aL*aL*((ML+one)*(ML+one)+(three-gamma)/(gamma-one)) &
                           + Q*half*aL*aL*ML*ML * max(zero, ML)
 endif

!Negative Part of Flux evaluated in the right cell.
 if     (MR >  one) then
     Fm(1) = zero
     Fm(2) = zero
     Fm(3) = zero
 elseif (MR < -one) then
     Fm    = physical_flux(uR)
 else
         Z = half * rhoR * aR / gamma * (MR-one)
         Q = rhoR * aR * (gamma-one) / gamma
     Fm(1) = Z             + Q       * min(zero, MR)
     Fm(2) = Z*aR*(MR-one) + Q*aR*MR * min(zero, MR)
     Fm(3) = Z*half*aR*aR*((MR-one)*(MR-one)+(three-gamma)/(gamma-one)) &
                           + Q*half*aR*aR*MR*MR * min(zero, MR)
 endif  

!Compute the flux: Fp(uL)+Fm(uR).
   StegerWarming = Fp + Fm

 end function StegerWarming
!-----------------------------------------------------------------------------


!*****************************************************************************
!* --- Van Leer's Flux-Vector Splitting ---
!*
!* B. van Leer, Flux Vector Splitting for the Euler Equations, Proc. 8th 
!* International Conference on Numerical Methods in Fluid Dynamics, Berlin, 
!* Springer Verlag, 1982.
!* 
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function VanLeer(uL,uR)
 implicit none
 real :: uL(3), uR(3) !  Input (conservative variables rho*[1, v, E])
 real :: VanLeer(3)   ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, quarter, half, one, two
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: aL, aR, ML, MR               ! Speeds of sound and Mach numbers.
 real :: Fp(3), Fm(3)                 ! F_plus and F_minus

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0
       two = 2.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      ML = vL/aL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      MR = vR/aR

!Positive Part of Flux evaluated in the left cell.
 Fp(1) =   quarter*rhoL*aL*(ML+one)*(ML+one)
 Fp(2) =   Fp(1)*two*aL*(one+half*(gamma-one)*ML)/gamma
 Fp(3) =   Fp(1)*two*aL*aL*(one+half*(gamma-one)*ML)**2/(gamma*gamma-one)

!Negative Part of Flux evaluated in the right cell.
 Fm(1) = - quarter*rhoR*aR*(MR-one)*(MR-one)
 Fm(2) =   Fm(1)*two*aR*(-one+half*(gamma-one)*MR)/gamma
 Fm(3) =   Fm(1)*two*aR*aR*(one-half*(gamma-one)*MR)**2/(gamma*gamma-one)

!Compute the flux: Fp(uL)+Fm(uR).
   VanLeer = Fp + Fm

 end function VanLeer
!-----------------------------------------------------------------------------


!*****************************************************************************
!* --- AUSM Flux Function ---
!*
!* M.-S. Liou and C. J. Steffen, A New Flux Splitting Scheme, Journal of 
!* Computational Physics, 107, pp. 23-39, 1993.
!*
!* NB: This may require a low CFL number to get started.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function AUSM(uL,uR)
 implicit none
 real :: uL(3), uR(3)  !  Input (conservative variables rho*[1, v, E])
 real :: AUSM(3)       ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, quarter, half, one, two
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: HL, HR                       ! Specific enthaply (per unit mass).
 real :: aL, aR, ML, MR               ! Speeds of sound and Mach numbers.
 real :: Pp, Pm, Mp, Mm
 real :: Fp(3), Fm(3)                 ! F_plus and F_minus

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0
       two = 2.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      ML = vL/aL
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      MR = vR/aR
      HR = ( uR(3) + pR ) / rhoR

!Positive M and p in the LEFT cell.
 if (ML <= -one) then
   Mp = zero
   Pp = zero
 elseif (ML < one) then
   Mp = quarter*(ML+one)*(ML+one)
   Pp = quarter*PL*(one+ML)*(one+ML)*(two-ML) ! or use Pp = half*(one+ML)*pL
 else
   Mp = ML
   Pp = PL
 endif

!Negative M and p in the RIGHT cell.
 if   (MR <= -one) then
   Mm = MR
   Pm = PR;
 elseif (MR < one) then
   Mm = -quarter*(MR-one)*(MR-one)
   Pm =  quarter*pR*(one-MR)*(one-MR)*(two+MR) ! or use Pm = half*(one-MR)*pR
 else
   Mm = zero
   Pm = zero
 endif

!Positive Part of Flux evaluated in the left cell.
 Fp(1) = max(zero,Mp+Mm)*aL * rhoL
 Fp(2) = max(zero,Mp+Mm)*aL * rhoL*vL  + Pp
 Fp(3) = max(zero,Mp+Mm)*aL * rhoL*HL

!Negative Part of Flux evaluated in the right cell.
 Fm(1) = min(zero,Mp+Mm)*aR * rhoR
 Fm(2) = min(zero,Mp+Mm)*aR * rhoR*vR  + Pm
 Fm(3) = min(zero,Mp+Mm)*aR * rhoR*HR

!Compute the flux: Fp(uL)+Fm(uR).
    AUSM = Fp + Fm

 end function AUSM
!-----------------------------------------------------------------------------


!*****************************************************************************
!* --- ZHA-BILGEN Flux Function ---
!*
!* G.-C. Zha and E. Bilgen, Numerical Solutions of Euler Equations by Using
!* a New Flux Vector Splitting Schemes, International Journal for Numerical
!* Methods in Fluids, 17, pp. 11-144, 1993.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function ZhaBilgen(uL,uR)
 implicit none
 real :: uL(3), uR(3)   !  Input (conservative variables rho*[1, v, E])
 real :: ZhaBilgen(3)   ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, quarter, half, one, two
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: EL, ER                       ! Specific total energy (per unit mass).
 real :: aL, aR, ML, MR               ! Speeds of sound and Mach numbers.
 real :: Pp, Pm, Mp, Mm, Pup, Pum
 real :: Fp(3), Fm(3)                 ! F_plus and F_minus

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0
       two = 2.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      ML = vL/aL
      EL = uL(3) / rhoL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      MR = vR/aR
      ER = uR(3) / rhoR

!Positive p and p*u in the LEFT cell.
 if (ML <= -one) then
   Pp  = zero
   Pup = zero
 elseif (ML < one) then
   Pp = quarter*PL*(one+ML)*(one+ML)*(two-ML) ! or use Pp = half*(one+ML)*pL
  Pup = half*(vL+aL)*pL
 else
   Pp = PL
  Pup = pL*vL
 endif

!Negative p and p*u in the RIGHT cell.
 if   (MR <= -one) then
   Pm = pR
  Pum = pR*vR
 elseif (MR < one) then
   Pm =  quarter*pR*(one-MR)*(one-MR)*(two+MR) ! or use Pm = half*(one-MR)*pR
  Pum =  half*(vR-aR)*pR
 else
   Pm = zero
  Pum = zero
 endif

!Positive Part of Flux evaluated in the left cell.
 Fp(1)=max(zero,vL)*rhoL
 Fp(2)=max(zero,vL)*rhoL*vL  + Pp
 Fp(3)=max(zero,vL)*rhoL*EL  + Pup

!Negative Part of Flux evaluated in the right cell.
 Fm(1)=min(zero,vR)*rhoR
 Fm(2)=min(zero,vR)*rhoR*vR  + Pm
 Fm(3)=min(zero,vR)*rhoR*ER  + Pum

!Compute the flux: Fp(uL)+Fm(uR).
  ZhaBilgen = Fp + Fm

 end function ZhaBilgen
!-----------------------------------------------------------------------------


!*****************************************************************************
!* --- The Godunov Flux ---
!*
!* S. K. Godunov, A Difference Scheme for Numerical Computation of 
!* Disconsinuous Solution of Hydrodynamic Equations, Math. Sbornik, 47,
!* pp. 271-306, 1959 (in Russian). Translated US Joint Publ. Res. Service,
!* JPRS 7226 (1969).
!*
!* See textbooks such as Toro's or Hirsh's books, or "I do Like CFD, VOL.1"
!* for details on the algorithm.
!*
!* This requires the functions, massflux(), sonic(), and physical_flux(u).
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function Godunov(uL,uR)
 implicit none
 real :: uL(3), uR(3) !  Input (conservative variables rho*[1, v, E])
 real :: Godunov(3)   ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, quarter, half, one, two
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: aL, aR                       ! Speeds of sound.
 real :: Fp(3), Fm(3)                 ! F_plus and F_minus

 real    :: gam,gam2,tol,pm1,pm2,mL,mR,vm,p(2),rm(2),r(2),rmI,amL,amR
 real    :: SmL,SmR,Um2,Um3
 real    :: um(3)
 integer :: k, kmax

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0
       two = 2.0
       tol = 1.0D-05
      gam2 = half*(gamma-one)/gamma

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
 
   gam = (gamma+one)/(gamma-one)
  r(1) = rhoL
  r(2) = rhoR
  P(1) = pL
  P(2) = pR;

!*** Supersonic Flow to the RIGHT
 if     (vL/aL >=  one) then

  Godunov = physical_flux(uL)
  return

!*** Supersonic Flow to the LEFT
 elseif (vR/aR <= -one) then 

  Godunov = physical_flux(uR)
  return

!*** Otherwise
 else

! Initial solution: (intersection of two linearized integral curves,
!                    which is actually the upper bound of the solution.)
   pm1 = (  (half*(vL-vR)*(gamma-one)+aL+aR)/( &
             aL*pL**((one-gamma)/gamma*half)   &
           + aR*pR**((one-gamma)/gamma*half) )  )**(two*gamma/(gamma-one))

! Fixed-point iteration to find the pressure and velocity in the middle.
! (i.e., find the intersection of two nonlinear integral curves.)

     k = 0
  kmax = 100

  do

    mL = massflux(rhoL,aL,pL,pm1)
    mR = massflux(rhoR,aR,pR,pm1)
   pm2 = (mL*pR+mR*pL-mL*mR*(vR-vL))/(mL+mR)
     k = k + 1
   pm1 = pm2

   if (abs(pm2-pm1) < tol) exit
   if (k > kmax) then
    write(*,*) " Fixed-point iteration did not converge. Stop."
    stop
   endif

  end do

  mL = massflux(rhoL,aL,pL,pm2)
  mR = massflux(rhoR,aR,pR,pm2)
  vm = (mL*vL+mR*vR-(pR-pL))/(mL+mR)

!Density in the middle.
 do k = 1, 2
  if (pm2/p(K) >= one) then
    rm(k) = r(k)*(one+gam*pm2/p(k))/(gam+pm2/p(k))
  else                    
    rm(k) = r(k)*( pm2/p(k) )**(one/gamma)
  endif
 end do

!Contact wave to the right or left?
 if (vm >= 0) then
  rmI = rm(1)
 else
  rmI = rm(2) 
 endif

!Wave speeds at the interface, x/t = 0.
    amL = sqrt(gamma*pm2/rm(1))
    amR = sqrt(gamma*pm2/rm(2))
    SmL = vm - amL
    SmR = vm + amR

!Sonic case.
   if (SmL<=zero .and. SmR>=zero) then
    Um2 = rmI*vm
    Um3 = Pm2/(gamma-one)+half*rmI*vm*vm
   elseif (SmL>zero .and. vL-aL<zero) then
    call sonic(rmI,Um2,Um3,vL,aL,PL,vm,amL,vL-aL,SmL)
   elseif (SmR<zero .and. vR+aR>zero) then
    call sonic(rmI,Um2,Um3,vR,aR,PR,vm,amR,vR+aR,SmR)
   endif

!Compute the flux: evaluate the physical flux at the interface (middle).
  um(1) = rmI
  um(2) = Um2
  um(3) = Um3
  Godunov = physical_flux(um)

 endif

 end function Godunov
!-----------------------------------------------------------------------------

!*****************************************************************************
!* --- Osher's Flux Function ---
!*
!* S. Osher, Riemann Solvers, the Entropy Condition and Difference 
!* Approximations, SIAM Journal of Numerical Analysis, 21, pp. 217-235, 1984.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function Osher(uL,uR)
 implicit none
 real :: uL(3), uR(3) !  Input (conservative variables rho*[1, v, E])
 real :: Osher(3)     ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, quarter, half, one, two
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: aL, aR                       ! Speeds of sound
 real :: pLm,rLm,uLm,aLm
 real :: pRm,uRm,rRm,aRm
 real :: sgnU,sgnUmC,sgnUpC
 real :: uLma(3),SL,SLm
 real :: uRma(3),SR,SRm
 real :: us(3)

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0
       two = 2.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
 
!Compute the intermediate states.
 !*** Lm (or j+1/3)****
   pLm = ( (half*(gamma-one)*(vR-vL)+(aR+aL)) &
       &   /( aL*(one+sqrt(rhoL/rhoR*(pR/pL)**(one/gamma)))) )**(two*gamma/(gamma-one))*pL
   rLm = (pLm/pL)**(one/gamma)*rhoL
   aLm = sqrt(gamma*pLm/rLm)
   uLm = vL-2/(gamma-1)*(aL-aLm)
  uLma(1) = rLm
  uLma(2) = rLm*uLm
  uLma(3) = pLm/(gamma-one) + half*rLm*uLm*uLm
 !*** Rm (or j+2/3)****
      pRm = pLm
      uRm = uLm
      rRm = (pRm/pR)**(one/gamma)*rhoR
      aRm = sqrt(gamma*pRm/rRm)
  uRma(1) = rRm
  uRma(2) = rRm*uRm
  uRma(3) = pRm/(gamma-one)+half*rRm*uRm*uRm

!Centered flux.
   Osher = half * ( physical_flux(uL) + physical_flux(uR) )

!Wave speeds.

   if     (uRm < zero) then
      sgnU = -one
   elseif (uRm > zero) then
      sgnU = one
   else
      sgnU = zero
   endif

   if     (vR-aR < zero) then
     sgnUmC = -one
   elseif (vR-aR > zero) then
     sgnUmC =  one
   else
     sgnUmC = zero
   endif

   if     (vL+aL < zero) then
     sgnUpC = -one
   elseif (vL+aL > zero) then
     sgnUpC =  one
   else
     sgnUpC = zero
   endif

!*
  Osher = Osher - half*(  -sgnUpC*physical_flux(uL) + sgnUmC*physical_flux(uR) )
  Osher = Osher - half*( ( sgnU - sgnUmC) * physical_flux(uRma) )
  Osher = Osher - half*( (-sgnU + sgnUpC) * physical_flux(uLma) )

!Take into account sonic cases.

 !*** The first Path: u+c    FUpC
     SLm = uLm+aLm
      SL = vL+aL  !*** wave speeds0.25*WL1*CL*(ML+1)*(ML+1)
  if (SL*SLm < zero) then
   call sonic(us(1),us(2),us(3),uLm,aLm,pLm,vL,aL,SLm,SL)
   Osher = Osher - half*( two*sgnUpC * (physical_flux(us)-physical_flux(uLma)) )
  endif

 !*** The third Path: u-c    FUmC
     SRm = uRm-aRm
      SR = vR-aR  !*** wave speeds
  if (SRm*SR < zero) then
   call sonic(us(1),us(2),us(3),vR,aR,pR,uRm,aRm,SR,SRm)
   Osher = Osher - half*( two*sgnUmC * (physical_flux(uRma)-physical_flux(us)) )
  endif

 end function Osher
!-----------------------------------------------------------------------------


!*****************************************************************************
!* -- Roe's Flux Function ---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!* 
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function Roe(uL,uR)
 implicit none
 real :: uL(3), uR(3) !  Input (conservative variables rho*[1, v, E])
 real :: Roe(3)       ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, quarter, half, one, two, four
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: aL, aR, HL, HR               ! Speeds of sound.
 real :: RT,rho,v,H,a                 ! Roe-averages
 real :: drho,du,dP,dV(3)
 real :: ws(3),Da, R(3,3)
 integer :: j, k

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0
       two = 2.0
      four = 4.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(3) + pR ) / rhoR

!First compute the Roe Averages **************************
    RT = sqrt(rhoR/rhoL);
   rho = RT*rhoL
     v = (vL+RT*vR)/(one+RT)
     H = (HL+RT*HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*v*v) )

!Differences in primitive variables.
   drho = rhoR - rhoL
     du =   vR - vL
     dP =   pR - pL

!Wave strength (Characteristic Variables).
   dV(1) =  half*(dP-rho*a*du)/(a*a)
   dV(2) = -( dP/(a*a) - drho )
   dV(3) =  half*(dP+rho*a*du)/(a*a)

!Absolute values of the wave speeds (Eigenvalues)
   ws(1) = abs(v-a)
   ws(2) = abs(v  )
   ws(3) = abs(v+a)

!Modified wave speeds for nonlinear fields (to remove expansion shocks).
!There are various ways to implement an entropy fix. This is just one
!example.
   Da = max(zero, four*((vR-aR)-(vL-aL)) )
   if (ws(1) < half*Da) ws(1) = ws(1)*ws(1)/Da + quarter*Da
   Da = max(zero, four*((vR+aR)-(vL+aL)) )
   if (ws(3) < half*Da) ws(3) = ws(3)*ws(3)/Da + quarter*Da

!Right eigenvectors
   R(1,1) = one
   R(2,1) = v - a
   R(3,1) = H - v*a

   R(1,2) = one
   R(2,2) = v
   R(3,2) = half*v*v

   R(1,3) = one
   R(2,3) = v + a
   R(3,3) = H + v*a

!Compute the average flux.
   Roe = half*( physical_flux(uL) + physical_flux(uR) )

!Add the matrix dissipation term to complete the Roe flux.
  do j = 1, 3
   do k = 1, 3
    Roe(j) = Roe(j) - half*ws(k)*dV(k)*R(j,k) 
   end do
  end do

 end function Roe
!-----------------------------------------------------------------------------

!*****************************************************************************
!* --- Rusanov's Flux Function ---
!* 
!* V. V. Rusanov, Calculation of Interaction of Non-Steady Shock Waves with
!* Obstacles, J. Comput. Math. Phys. USSR, 1, pp. 267-279, 1961.
!*
!* This requires the function, physical_flux(u), included in this file.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function Rusanov(uL,uR)
 implicit none
 real :: uL(3), uR(3)   !  Input (conservative variables rho*[1, v, E])
 real :: Rusanov(3)     ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: half, one
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: HL, HR
 real :: RT,rho,v,H,a                 ! Roe-averages
 real :: smax

!Constants.
     gamma = 1.4
      half = 0.5
       one = 1.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      HR = ( uR(3) + pR ) / rhoR

!First compute the Roe Averages **************************
    RT = sqrt(rhoR/rhoL);
   rho = RT*rhoL
     v = (vL+RT*vR)/(one+RT)
     H = (HL+RT*HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*v*v) )

    smax = abs(v)+a

!Rusanov Flux (This is very similar to LaxFriedrichs, isn't it?)
 Rusanov = half*(physical_flux(uR) + physical_flux(uL) - smax*(uR-uL))

 end function Rusanov
!-----------------------------------------------------------------------------


!*****************************************************************************
!* --- HLL Flux Function ---
!*
!* A. Harten, P. D. Lax, and B. van Leer,On Upstream Differencing and 
!* Godunov-Type Schemes for Hyperbolic Conservation Laws, SIAM Review,
!* 25(1), pp. 35-61, 1983.
!*
!* With wave speeds evaluated by Einfeldt's method:
!* B. Einfeldt, On Godunov-Type Methods for Gas Dynamics, SIAM Journal of
!* Numerical Analysis, 25(2), pp. 294-318, 1988.
!* 
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function HLL(uL,uR)
 implicit none
 real :: uL(3), uR(3)     !  Input (conservative variables rho*[1, v, E])
 real :: HLL(3)           ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, quarter, half, one
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: aL, aR, HL, HR               ! Speeds of sound.
 real :: RT,rho,u,H,a
 real :: SRp,SLm, uma,upa

!Constants.
     gamma = 1.4
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(3) + pR ) / rhoR

!Evaluate the two wave speeds: Einfeldt.
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
     u = (vL+RT*vR)/(one+RT)
     H = (HL+RT*HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*u*u) )

   uma = u - a
   upa = u + a
   SLm = min(vL-aL, uma, zero)
   SRp = max(vR+aR, upa, zero)

!Compute the HLL flux.
   HLL = (SRp*physical_flux(uL)-SLm*physical_flux(uR) + SLm*SRp*(uR-uL))/(SRp-SLm)

 end function HLL
!-----------------------------------------------------------------------------


!*****************************************************************************
!* --- HLLL Flux Function ---
!*
!* T. Linde, A Practical, Geeral-Purpose, Two-State HLL Riemann Solver for
!* Hyperbolic COnservation Laws, International Journal for Numerical Methods
!* in Fluids, 40, pp. 391-402, 2002.
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function HLLL(uL,uR)
 implicit none
 real :: uL(3), uR(3) !  Input (conservative variables rho*[1, v, E])
 real :: HLLL(3)      ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma, gm1                   ! Ratio of specific heat.
 real :: zero, quarter, half, one
!Local variables
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: aL, aR                       ! Speeds of sound.
 real :: WsL(3),WsR(3),DWs(3),DF(3),FL(3),FR(3),DU(3)
 real :: lm_p,lm_m, alpha,V,Vp,Vm, DWsDF,DWsDU, Pdiag(3),DUDF,DU2,DF2
 real :: pm,um,rm, Hugoniot,sL,sR,s_flux_R,s_flux_L

!Constants.
     gamma = 1.4
       gm1 = 1.4 - 1.0
      zero = 0.0
    quarter = 0.25
      half = 0.5
       one = 1.0

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)


 FL = physical_flux(uL)
 FR = physical_flux(uR)
 
 WsL =(/ PL/((gamma-one)*rhoL)*(gamma-log(PL/rhoL**gamma))-HALF*vL*vL, vL, -one /)
 WsR =(/ PR/((gamma-one)*rhoR)*(gamma-log(PR/rhoR**gamma))-HALF*vR*vR, vR, -one /)
     DWs = WsR - WsL 
      DF =  FR -  FL
      DU =  UR -  UL
   DWsDF = DWs(1)*DF(1) + DWs(2)*DF(2) + DWs(3)*DF(3)
   DWsDU = DWs(1)*DU(1) + DWs(2)*DU(2) + DWs(3)*DU(3)
       V = DWsDF/(DWsDU+1.0e-25)

!Entropy
           sL = -rhoL*log(PL/rhoL**gamma)
           sR = -rhoR*log(PR/rhoR**gamma)
     s_flux_L = vL*sL
     s_flux_R = vR*sR
     Hugoniot = V*(sR-sL)-(s_flux_R-s_flux_L)

   if (Hugoniot < 1.0e-02) then
      alpha = zero
   else
         rm = HALF*(rhoL+rhoR)
         pm = HALF*(PL+PR)
         um = HALF*(vL+vR)
      Pdiag = (/ quarter*um**4+gamma*pm*pm/gm1/gm1/rm/rm,  um*um+pm/gm1/rm,  one /)
       DUDF = DU(1)*Pdiag(1)*DF(1) + DU(2)*Pdiag(2)*DF(2) + DU(3)*Pdiag(3)*DF(3)
       DU2  = DU(1)*Pdiag(1)*DU(1) + DU(2)*Pdiag(2)*DU(2) + DU(3)*Pdiag(3)*DU(3)
       DF2  = DF(1)*Pdiag(1)*DF(1) + DF(2)*Pdiag(2)*DF(2) + DF(3)*Pdiag(3)*DF(3)
      alpha = DUDF*DUDF/(DU2*DF2 + 1.0e-25);
   endif

     lm_m = min(V, vL-aL, vR-aR, zero)
     lm_p = max(V, vL+aL, vR+aR, zero)
       Vm = min(V, ZERO)
       Vp = max(V, ZERO)

   HLLL = ( (lm_p*FL-lm_m*FR) + ((one-alpha)*lm_m*lm_p+alpha*(lm_m*Vp+lm_p*Vm)) *DU ) /(lm_p-lm_m)

 end function HLLL
!-----------------------------------------------------------------------------

!*****************************************************************************
!* --- AUFS Flux Function ---
!*
!* M. Sun and K. Takayama, An artificially upstream flux vector splitting
!* scheme for the Euler equations.
!* Journal of Computational Physics, 189, pp. 305-329, 2003.
!*
!* Katate Masatsuka, February 2010. http://www.cfdbooks.com
!*****************************************************************************
 function AUFS(uL,uR)
 implicit none
 real :: uL(3), uR(3)!  Input (conservative variables rho*[1, v, E])
 real :: AUFS(3)     ! Output (numerical flux across L and R states)
!Local constants
 real :: gamma                        ! Ratio of specific heat.
 real :: zero, quarter, half, one, two
!Local variables
 real :: gm1, vstar
 real :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real :: HL, HR                       ! Total enthalpy
 real :: aL, aR, astar                ! Speeds of sound
 real :: ca                           ! Intermediate speed of sound
 real :: s1, s2, M                    ! Artificial wave speeds, weight
 real :: F1(3), F2(3), del_U(3)       ! Split fluxes, and dissipation

!Constants.
     gamma = 1.4
      zero = 0.0
   quarter = 0.25
      half = 0.5
       one = 1.0
       two = 2.0
       gm1 = gamma-1

!Primitive and other variables.
!  Left state
    rhoL = uL(1)
      vL = uL(2)/uL(1)
      pL = (gamma-one)*( uL(3) - half*rhoL*vL*vL )
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = uR(1)
      vR = uR(2)/uR(1)
      pR = (gamma-one)*( uR(3) - half*rhoR*vR*vR )
      aR = sqrt(gamma*pR/rhoR)
      HR = ( uR(3) + pR ) / rhoR

!Compute F1: The Steger-Warming flux with an intermediate speed of sound.
        ca = half*(aL+aR)
  del_U(1) = half*(pL-pR)/ca
  del_U(2) = half*(vL*pL-vR*pR)/ca
  del_U(3) = half*(ca/gm1*(pL-pR)+half*(pL*vL*vL-pR*vR*vR)/ca)

  F1(1) =                          del_U(1)
  F1(2) = half*(    pL + pR    ) + del_U(2)
  F1(3) = half*( vL*pL + vR*pR ) + del_U(3)

!Isentropic wave estimates
     vstar = half*(vL+vR) + (aL-aR)/gm1
     astar = half*(aL+aR) + quarter*gm1*(vL-vR)

!Compute F2: Fully one-sided flux (due to the special wave structure).
        s1 = half*(vL+vR)
    if(s1 > zero)then
          s2 = min(zero, vL-aL, vstar-astar)
           M = s1/(s1-s2)
       F2(1) = uL(1)*(vL-s2)
       F2(2) = uL(2)*(vL-s2) + pL
       F2(3) = uL(3)*(vL-s2) + pL*vL
   else
          s2 = max(zero,vstar+astar,vR+aR)
           M = s1/(s1-s2)
       F2(1) = uR(1)*(vR-s2)
       F2(2) = uR(2)*(vR-s2) + pR
       F2(3) = uR(3)*(vR-s2) + pR*vR
   end if

!Final form of the numerical flux (weighted average of F1 and F2).
  AUFS = (one-M)*F1 + M*F2

 end function AUFS
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
! Below are functions and subroutines used by the numerical fluxes above.
!-----------------------------------------------------------------------------

!*****************************************************************************
!* --- 1D physical Euler flux --- 
!*
!* This is called in LaxFriedrichs, Richtmyer, MacCormack, StegerWarming,
!* Godunov, Osher, Roe, HLL, HLLL.
!*
!* The vector f(U) in u_t + f(u)_x = 0, as a function of 
!* the conservative variables, u = [density, density*velocity, total energy].
!
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
 function physical_flux(u)
 real :: u(3)             !  Input (conservative variables [rho, rho*v, rho*E])
 real :: physical_flux(3) ! Output (physical flux of the Euler equations)
!Local variables
 real :: density, velocity, pressure, enthalpy, gamma

!Define and compute some quantities.
     gamma = 1.4
   density = u(1)
  velocity = u(2)/u(1)
  pressure = (gamma-1.0)*( u(3) - 0.5*density*velocity*velocity )
  enthalpy = u(3) + pressure

!Evaluate the physical flux (mass, momentum, and energy fluxes).
  physical_flux(1) =           density * velocity
  physical_flux(2) = (density*velocity)* velocity + pressure
  physical_flux(3) =          enthalpy * velocity

 end function physical_flux
!*****************************************************************************


!****************************************************************************
!* --- Mass flux used for Godunov's Flux Function ---
!*
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
  function massflux(r,c,pQ,pm)
  real :: r,c,pQ,pm !  Input
  real :: massflux  ! Output
!Local variables
  real :: gamma,gam1,gam2,one,half,eps
    one = 1.0
   half = 0.5
  gamma = 1.4
    eps = 1.0e-15

   gam1=half*(gamma+one)/gamma; gam2=half*(gamma-one)/gamma;
   if (pm/pQ >= one-eps) then ! eps to avoid zero-division
      massflux=(r*c)*sqrt( one+gam1*(pm/pQ-one) )
   else
      massflux=(r*c)*gam2*(one-pm/pQ)/( ONE-(pm/pQ)**(gam2) )
   endif
  end function massflux
!*****************************************************************************

!****************************************************************************
!* --- Solutions at Sonic points --- used in Godunov's and Osher's fluxes.
!
!* Katate Masatsuka, February 2009. http://www.cfdbooks.com
!*****************************************************************************
  subroutine sonic(US1,US2,US3, u1,c1,P1,u2,c2,a1,a2)
  real :: u1,c1,P1,u2,c2,a1,a2 ! Input
  real :: US1,US2,US3          ! Output
!Local variables
  real              :: us,cs,Ps,rs,R1,R2    
  real :: gamma, one, two, half
  gamma = 1.4
    one = 1.0
    two = 2.0
   half = 0.5

    R1 =  a2/(a2-a1)
    R2 = -a1/(a2-a1)
    us = R1*u1+R2*u2
    cs = R1*c1+R2*c2
    Ps = (cs/c1)**(two*gamma/(gamma-1))*P1
    rs = gamma*Ps/(cs*cs)

   US1 = rs
   US2 = rs*us
   US3 = Ps/(gamma-one)+half*rs*us*us

  end subroutine sonic
!*****************************************************************************
