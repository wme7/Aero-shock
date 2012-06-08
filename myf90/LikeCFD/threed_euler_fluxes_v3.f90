!********************************************************************************
!* A collection of three-dimensional Euler numerical fluxes, Version 3 (2012),
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!*------------------------
!* List of Flux Functions:
!*
!*   inviscid_roe          : Roe flux (using face normal and tangent vectors)
!*   inviscid_roe_n        : Roe flux (using face normal vector only)
!*   inviscid_rotated_rhll : Rotated-RHLL flux (one of the robust fluxes for carbuncle)
!*
!*------------------------
!*
!* These F90 routines were written and made available for an educational purpose.
!* Detailed descripstion of each numerical flux can be found in the original paper,
!* or in popular textbooks, or in the (will-be-available) second volume of
!* "I do like CFD". 
!*
!* Note: These subroutines have been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*       One way to make it efficient is to reduce the number of local variables,
!*       by re-using temporary variables as many times as possible.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!*
!* Version 01
!* 05-06-11: Minor bugs fixed.
!*
!* Version 02
!* 04-06-12: Major bugs fixed for Rotated-RHLL flux.
!*           - qL and qR were not computed with the correct normals
!*           - Tangent vector (2nd direction) is now computed more carefully.
!*           (Bugs reported by Dr. Yoshitaka Nakashima. Thank you!)
!*
!* Version 03
!* 04-16-12: Equal sign has been added in the if-statement for tangent vector to
!*           avoid a possible failure (reported by Dr. Yoshitaka Nakashima. Thank you!)
!*           Roe flux without tangent vectors (inviscid_roe_n) added.
!*           Many more comments added.
!*
!* This file may be updated in future. (to add more flux functions.)
!********************************************************************************


!********************************************************************************
!* -- 3D Roe's Flux Function with an entropy fix --
!*
!* This subroutine computes the Roe flux for the Euler equations
!* in the direction, njk=[nx,ny,nz].
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* Conservative form of the Euler equations:
!*
!*     dU/dt + dF/dx + dG/dy + dH/dz = 0
!*
!* This subroutine computes the numerical flux for the flux in the direction,
!* njk=[nx,ny,nz]:
!*
!*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
!*                               | rho*qn*u + p*nx |
!*                               | rho*qn*v + p*ny |
!*                               | rho*qn*w + p*nz |
!*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
!*
!* The Roe flux is implemented in the following form:
!*
!*   Numerical flux = 1/2 [ Fn(UR) + Fn(UL) - |An|dU ], 
!*
!*  where
!*
!*    An = dFn/dU,  |An| = R|Lambda|L, dU = UR - UL.
!*
!* The dissipation term, |An|dU, is actually computed as
!*
!*     sum_{k=1,5} |lambda_k| * (LdU)_k * r_k,
!*
!* where lambda_k is the k-th eigenvalue, (LdU)_k is the k-th wave strength,
!* and r_k is the k-th right-eigenvector.
!*
!* ------------------------------------------------------------------------------
!*  Input: primL(1:5) =  Left state (rhoL, uL, vL, wR, pL)
!*         primR(1:5) = Right state (rhoR, uR, vR, wR, pR)
!*           njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right.
!*
!*           njk
!*  Face normal ^   o Right data point
!*              |  .
!*              | .
!*              |. 
!*       -------x-------- Face
!*             .                 Left and right states are
!*            .                   1. Values at data points for 1st-order accuracy
!*           .                    2. Extrapolated values at the face midpoint 'x'
!*          o Left data point        for 2nd/higher-order accuracy.
!*
!*
!* Output:  flux(1:5) = The Roe flux with an entropy fix
!* ------------------------------------------------------------------------------
!*
!* Note: This subroutine has been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*       One way to make it efficient is to reduce the number of local variables,
!*       by re-using temporary variables as many times as possible.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!*
!* Katate Masatsuka, March 2011. http://www.cfdbooks.com
!********************************************************************************
 subroutine inviscid_roe(primL, primR, njk,  num_flux)

 implicit none
 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

!Input
 real(p2), intent( in) :: primL(5), primR(5) ! Input: primitive variables
 real(p2), intent( in) :: njk(3)             ! Input: face normal vector

!Output
 real(p2), intent(out) :: num_flux(5)        ! Output: numerical flux

!Some constants
 real(p2) ::  zero = 0.0_p2
 real(p2) ::   one = 1.0_p2
 real(p2) ::   two = 2.0_p2
 real(p2) ::  half = 0.5_p2
 real(p2) :: fifth = 0.2_p2

!Local variables
 real(p2) :: nx, ny, nz                   ! Normal vector
 real(p2) :: mx, my, mz                   ! Orthogonal tangent vector
 real(p2) :: lx, ly, lz                   ! Another orthogonal tangent vector
 real(p2) :: abs_n_cross_l                ! Magnitude of n x l
 real(p2) :: uL, uR, vL, vR, wL, wR       ! Velocity components.
 real(p2) :: rhoL, rhoR, pL, pR           ! Primitive variables.
 real(p2) :: qnL, qnR, qmL, qmR, qlL, qlR ! Normal and tangent velocities
 real(p2) :: aL, aR, HL, HR               ! Speed of sound, Total enthalpy
 real(p2) :: RT,rho,u,v,w,H,a,qn, ql, qm  ! Roe-averages
 real(p2) :: drho,dqn,dql,dqm,dp,LdU(5)   ! Wave strengths
 real(p2) :: ws(5), R(5,5)                ! Wave speeds and right-eigenvectors
 real(p2) :: dws(5)                       ! Width of a parabolic fit for entropy fix
 real(p2) :: fL(5), fR(5), diss(5)        ! Fluxes ad dissipation term
 real(p2) :: gamma = 1.4_p2               ! Ratio of specific heats
 real(p2) :: temp, tempx, tempy, tempz    ! Temoprary variables

! Face normal vector (unit vector)

  nx = njk(1)
  ny = njk(2)
  nz = njk(3)

!Define face tangent vectors, m and l, such that the three vectors, n, m, and l
!are mutually orthogonal.

! l = (lx.ly,lz)
!
! We first carefully choose one tangent vector: l = (lx,ly,lz).
!  E.g., if n = (1,0,0), then, (0,-nz,ny) = (0,0,0), which is a problem.
!        To avoid such a situation, we choose the 2D tangential vector
!        having the maximum length.

     tempx = ny*ny + nz*nz
     tempy = nz*nz + nx*nx
     tempz = nx*nx + ny*ny

     if     ( tempx >= tempy .and. tempx >= tempz ) then
       lx =  zero
       ly = -nz
       lz =  ny
     elseif ( tempy >= tempx .and. tempy >= tempz ) then
       lx = -nz
       ly =  zero
       lz =  nx
     elseif ( tempz >= tempx .and. tempz >= tempy ) then
       lx = -ny
       ly =  nx
       lz =  zero
     else
      ! Impossible to happen
      write(*,*) "subroutine inviscid_roe: Impossible to happen. Please report the problem."
      stop
     endif

!     Make it the unit vector.
      temp = sqrt( lx*lx + ly*ly + lz*lz )
       lx = lx/temp
       ly = ly/temp
       lz = lz/temp

! m = (mx,my,mz)
!
! The other one, m = (mx,my,mz), is chosen as a vector orthogonal to both n and l
! defined by the vector product: m = n x l / |n x l|

  mx = ny*lz - nz*ly
  my = nz*lx - nx*lz
  mz = nx*ly - ny*lx

  abs_n_cross_l = sqrt(mx**2 + my**2 + mz**2)
  mx = mx / abs_n_cross_l
  my = my / abs_n_cross_l
  mz = mz / abs_n_cross_l

!(Do you like such ambiguous tangent vectors? Actually, the Roe flux can
! be implemented without any tangent vector. See "I do like CFD, VOL.1",
! or the subroutine "inviscid_roe_n" or "inviscid_rotated_rhll" for details.
! In fact, the resulting flux is independent of the choice of these tangent vectors
! as it should be.)

!Primitive and other variables.

!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
      wL = primL(4)
     qnL = uL*nx + vL*ny + wL*nz
     qlL = uL*lx + vL*ly + wL*lz
     qmL = uL*mx + vL*my + wL*mz
      pL = primL(5)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
      wR = primR(4)
     qnR = uR*nx + vR*ny + wR*nz
     qlR = uR*lx + vR*ly + wR*lz
     qmR = uR*mx + vR*my + wR*mz
      pR = primR(5)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

!First compute the Roe-averaged quantities

!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL                                        !Roe-averaged density
     u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
     v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
     w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
     H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
     a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
    qn = u*nx + v*ny + w*nz                             !Roe-averaged face-normal velocity
    ql = u*lx + v*ly + w*lz                             !Roe-averaged face-tangent velocity
    qm = u*mx + v*my + w*mz                             !Roe-averaged face-tangent velocity

!Wave Strengths

   drho = rhoR - rhoL !Density difference
     dp =   pR - pL   !Pressure difference
    dqn =  qnR - qnL  !Normal velocity difference
    dql =  qlR - qlL  !Tangent velocity difference in l
    dqm =  qmR - qmL  !Tangent velocity difference in m

  LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            !Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
  LdU(4) = rho*dql                     !Shear wave strength
  LdU(5) = rho*dqm                     !Shear wave strength

!Absolute values of the wave speeds

  ws(1) = abs(qn-a) !Left-moving acoustic wave speed
  ws(2) = abs(qn)   !Entropy wave speed
  ws(3) = abs(qn+a) !Right-moving acoustic wave speed
  ws(4) = abs(qn)   !Shear wave speed
  ws(5) = abs(qn)   !Shear wave speed

!Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
!NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!Right Eigenvectors

! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx
  R(3,1) = v - a*ny
  R(4,1) = w - a*nz
  R(5,1) = H - a*qn

! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)

! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx
  R(3,3) = v + a*ny
  R(4,3) = w + a*nz
  R(5,3) = H + a*qn

! Shear wave
  R(1,4) = zero
  R(2,4) = lx
  R(3,4) = ly
  R(4,4) = lz
  R(5,4) = ql

! Shear wave
  R(1,5) = zero
  R(2,5) = mx
  R(3,5) = my
  R(4,5) = mz
  R(5,5) = qm

!Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]

 diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) + ws(3)*LdU(3)*R(:,3) &
         + ws(4)*LdU(4)*R(:,4) + ws(5)*LdU(5)*R(:,5)

!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]

  num_flux = half * (fL + fR - diss)

!Normal max wave speed in the normal direction.
!  wsn = abs(qn) + a

 end subroutine inviscid_roe
!--------------------------------------------------------------------------------


!********************************************************************************
!* -- 3D Roe's Flux Function with an entropy fix and without tangent vectors --
!*
!* NOTE: This version does not use any tangent vector.
!*       See "I do like CFD, VOL.1" about how tangent vectors are eliminated.
!*
!* This subroutine computes the Roe flux for the Euler equations
!* in the direction, njk=[nx,ny,nz].
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* Conservative form of the Euler equations:
!*
!*     dU/dt + dF/dx + dG/dy + dH/dz = 0
!*
!* This subroutine computes the numerical flux for the flux in the direction,
!* njk=[nx,ny,nz]:
!*
!*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
!*                               | rho*qn*u + p*nx |
!*                               | rho*qn*v + p*ny |
!*                               | rho*qn*w + p*nz |
!*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
!*
!* The Roe flux is implemented in the following form:
!*
!*   Numerical flux = 1/2 [ Fn(UR) + Fn(UL) - |An|dU ], 
!*
!*  where
!*
!*    An = dFn/dU,  |An| = R|Lambda|L, dU = UR - UL.
!*
!* The dissipation term, |An|dU, is actually computed as
!*
!*     sum_{k=1,5} |lambda_k| * (LdU)_k * r_k,
!*
!* where lambda_k is the k-th eigenvalue, (LdU)_k is the k-th wave strength,
!* and r_k is the k-th right-eigenvector.
!*
!* ------------------------------------------------------------------------------
!*  Input: primL(1:5) =  Left state (rhoL, uL, vL, wR, pL)
!*         primR(1:5) = Right state (rhoR, uR, vR, wR, pR)
!*           njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right.
!*
!*           njk
!*  Face normal ^   o Right data point
!*              |  .
!*              | .
!*              |. 
!*       -------x-------- Face
!*             .                 Left and right states are
!*            .                   1. Values at data points for 1st-order accuracy
!*           .                    2. Extrapolated values at the face midpoint 'x'
!*          o Left data point        for 2nd/higher-order accuracy.
!*
!*
!* Output:  flux(1:5) = The Roe flux with an entropy fix
!* ------------------------------------------------------------------------------
!*
!* Note: This subroutine has been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*       One way to make it efficient is to reduce the number of local variables,
!*       by re-using temporary variables as many times as possible.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!*
!* Katate Masatsuka, April 2012. http://www.cfdbooks.com
!********************************************************************************
 subroutine inviscid_roe_n(primL, primR, njk,  num_flux)

 implicit none
 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

!Input
 real(p2), intent( in) :: primL(5), primR(5) ! Input: primitive variables
 real(p2), intent( in) :: njk(3)             ! Input: face normal vector

!Output
 real(p2), intent(out) :: num_flux(5)        ! Output: numerical flux

!Some constants
 real(p2) ::  zero = 0.0_p2
 real(p2) ::   one = 1.0_p2
 real(p2) ::   two = 2.0_p2
 real(p2) ::  half = 0.5_p2
 real(p2) :: fifth = 0.2_p2

!Local variables
 real(p2) :: nx, ny, nz                   ! Normal vector
 real(p2) :: uL, uR, vL, vR, wL, wR       ! Velocity components.
 real(p2) :: rhoL, rhoR, pL, pR           ! Primitive variables.
 real(p2) :: qnL, qnR                     ! Normal velocities
 real(p2) :: aL, aR, HL, HR               ! Speed of sound, Total enthalpy
 real(p2) :: RT,rho,u,v,w,H,a,qn          ! Roe-averages
 real(p2) :: drho,dqn,dp,LdU(4)           ! Wave strengths
 real(p2) :: du, dv, dw                   ! Velocity differences
 real(p2) :: ws(4), R(5,4)                ! Wave speeds and right-eigenvectors
 real(p2) :: dws(4)                       ! Width of a parabolic fit for entropy fix
 real(p2) :: fL(5), fR(5), diss(5)        ! Fluxes ad dissipation term
 real(p2) :: gamma = 1.4_p2               ! Ratio of specific heats

! Face normal vector (unit vector)

  nx = njk(1)
  ny = njk(2)
  nz = njk(3)

!Primitive and other variables.

!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
      wL = primL(4)
     qnL = uL*nx + vL*ny + wL*nz
      pL = primL(5)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
      wR = primR(4)
     qnR = uR*nx + vR*ny + wR*nz
      pR = primR(5)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

!First compute the Roe-averaged quantities

!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL                                        !Roe-averaged density
     u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
     v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
     w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
     H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
     a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound
    qn = u*nx + v*ny + w*nz                             !Roe-averaged face-normal velocity

!Wave Strengths

   drho = rhoR - rhoL !Density difference
     dp =   pR - pL   !Pressure difference
    dqn =  qnR - qnL  !Normal velocity difference

  LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            !Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
  LdU(4) = rho                         !Shear wave strength (not really, just a factor)

!Absolute values of the wave Speeds

  ws(1) = abs(qn-a) !Left-moving acoustic wave
  ws(2) = abs(qn)   !Entropy wave
  ws(3) = abs(qn+a) !Right-moving acoustic wave
  ws(4) = abs(qn)   !Shear waves

!Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
!NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!Right Eigenvectors
!Note: Two shear wave components are combined into one, so that tangent vectors
!      are not required. And that's why there are only 4 vectors here.
!      See "I do like CFD, VOL.1" about how tangent vectors are eliminated.

! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx
  R(3,1) = v - a*ny
  R(4,1) = w - a*nz
  R(5,1) = H - a*qn

! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)

! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx
  R(3,3) = v + a*ny
  R(4,3) = w + a*nz
  R(5,3) = H + a*qn

! Two shear wave components combined into one (wave strength incorporated).
  du = uR - uL
  dv = vR - vL
  dw = wR - wL
  R(1,4) = zero
  R(2,4) = du - dqn*nx
  R(3,4) = dv - dqn*ny
  R(4,4) = dw - dqn*nz
  R(5,4) = u*du + v*dv + w*dw - qn*dqn

!Dissipation Term: |An|(UR-UL) = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ]

 diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) &
         + ws(3)*LdU(3)*R(:,3) + ws(4)*LdU(4)*R(:,4)

!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

! This is the numerical flux: Roe flux = 1/2 *[  Fn(UL)+Fn(UR) - |An|(UR-UL) ]

  num_flux = half * (fL + fR - diss)

!Normal max wave speed in the normal direction.
!  wsn = abs(qn) + a

 end subroutine inviscid_roe_n
!--------------------------------------------------------------------------------


!********************************************************************************
!* -- 3D Rotated-Roe-HLL (Rotated-RHLL) Flux Function ---
!*
!* This subroutine computes the Rotated-RHLL flux for the Euler equations
!* in the direction, njk=[nx,ny,nz].
!*
!* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
!* Resolving, Rotated-Hybrid Riemann Solvers,
!* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
!*
!* Very robust for nonlinear instability (carbuncle).
!* Recommended for high-speed flows involving strong shocks.
!* It is also capable of resolving boundary layers.
!*
!*
!* Conservative form of the Euler equations:
!*
!*     dU/dt + dF/dx + dG/dy + dH/dz = 0
!*
!* This subroutine computes the numerical flux for the flux in the direction,
!* njk=[nx,ny,nz]:
!*
!*     Fn = F*nx + G*ny + H*nz = | rho*qn          |
!*                               | rho*qn*u + p*nx |
!*                               | rho*qn*v + p*ny |
!*                               | rho*qn*w + p*nz |
!*                               | rho*qn*H        |    (qn = u*nx + v*ny + w*nz)
!*
!* The Rotated-RHLL flux is defined as
!*
!*   Numerical flux = alpha1*HLL(n1) + alpha2*Roe(n2),  if |dq| > eps
!*                    Roe(n)                         ,  if |dq| < eps
!*
!* where n1 is taken as the velocity difference vector (normal to shock or
!* tangent to shear waves), n2 is a vector perpendicular to n1,
!* alpha1 = n*n1, alpha2 = n*n2. That is, we decompose the normal vector, n,
!* in the two directions, n1 and n2, and apply the HLL in n1 and the Roe in n2.
!* However, the resulting flux can be simplified and can be made to look like
!* a single Roe-type flux. The additional cost over the Roe flux is reported
!* to be as small as 14% of the Roe flux.
!*
!* Note: HLL is introduced only near discontinuities and only by a fraction, alpha1.
!* Note: For full convergence to machine zero for staedy state computations,
!*       the vectors, n1 and n2, may have to be frozen. For time-accurate 
!        calculation, it is not necessary.
!* Note: Here, the Roe flux is computed without tangent vectors.
!*
!* ------------------------------------------------------------------------------
!*  Input: primL(1:5) =  Left state (rhoL, uL, vL, wR, pL)
!*         primR(1:5) = Right state (rhoR, uR, vR, wR, pR)
!*           njk(1:3) = unit face normal vector (nx, ny, nz), pointing from Left to Right
!*
!*           njk
!*  Face normal ^   o Right data point
!*              |  .
!*              | .
!*              |. 
!*       -------x-------- Face
!*             .                 Left and right states are
!*            .                   1. Values at data points for 1st-order accuracy
!*           .                    2. Extrapolated values at the face midpoint 'x'
!*          o Left data point        for 2nd/higher-order accuracy.
!*
!* Output:  flux(1:5) = the Roe flux
!* ------------------------------------------------------------------------------
!*
!* Note: This subroutine has been prepared for an educational purpose.
!*       It is not at all efficient. Think about how you can optimize it.
!*
!* Note: Please let me know if you find bugs. I'll greatly appreciate it and
!*       fix the bugs.
!* 
!* Katate Masatsuka, March 2011. http://www.cfdbooks.com
!********************************************************************************
 subroutine inviscid_rotated_rhll(primL, primR, njk,  num_flux)

 implicit none
 integer , parameter :: p2 = selected_real_kind(15) ! Double precision

!Input
 real(p2), intent( in) :: primL(5), primR(5) ! Input: primitive variables
 real(p2), intent( in) :: njk(3)             ! Input: face normal vector

!Output
 real(p2), intent(out) :: num_flux(5)        ! Output: numerical flux

!Some constants
 real(p2) ::  zero = 0.0_p2
 real(p2) ::   one = 1.0_p2
 real(p2) ::   two = 2.0_p2
 real(p2) ::  half = 0.5_p2
 real(p2) :: fifth = 0.2_p2

!Local variables
 real(p2) :: nx, ny, nz                     ! Face normal vector
 real(p2) :: uL, uR, vL, vR, wL, wR         ! Velocity components.
 real(p2) :: rhoL, rhoR, pL, pR             ! Primitive variables.
 real(p2) :: qnL, qnR                       ! Normal velocity
 real(p2) :: aL, aR, HL, HR                 ! Speed of sound, Total enthalpy
 real(p2) :: RT,rho,u,v,w,H,a,qn            ! Roe-averages
 real(p2) :: drho,dqn,dp,LdU(4)             ! Wave strengths
 real(p2) :: du, dv, dw                     ! Velocity conponent differences
 real(p2) :: eig(4)                         ! Eigenvalues
 real(p2) :: ws(4), R(5,4)                  ! Absolute Wave speeds and right-eigenvectors
 real(p2) :: dws(4)                         ! Width of a parabolic fit for entropy fix
 real(p2) :: fL(5), fR(5), diss(5)          ! Fluxes ad dissipation term

 real(p2) :: gamma = 1.4_p2                 ! Ratio of specific heats

 real(p2) :: SRp,SLm                        ! Wave speeds for the HLL part
 real(p2) :: nx1, ny1, nz1                  ! Vector along which HLL is applied
 real(p2) :: nx2, ny2, nz2                  ! Vector along which Roe is applied
 real(p2) :: alpha1, alpha2                 ! Projections of the new normals
 real(p2) :: abs_dq                         ! Magnitude of the velocity difference
 real(p2) :: temp, tempx, tempy, tempz      ! Temporary variables

! Face normal vector (unit vector)

  nx = njk(1)
  ny = njk(2)
  nz = njk(3)

!Primitive and other variables.

!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
      wL = primL(4)
     qnL = uL*nx + vL*ny + wL*nz
      pL = primL(5)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL+wL*wL)

!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
      wR = primR(4)
     qnR = uR*nx + vR*ny + wR*nz
      pR = primR(5)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR+wR*wR)

!Compute the physical flux: fL = Fn(UL) and fR = Fn(UR)

  fL(1) = rhoL*qnL
  fL(2) = rhoL*qnL * uL + pL*nx
  fL(3) = rhoL*qnL * vL + pL*ny
  fL(4) = rhoL*qnL * wL + pL*nz
  fL(5) = rhoL*qnL * HL

  fR(1) = rhoR*qnR
  fR(2) = rhoR*qnR * uR + pR*nx
  fR(3) = rhoR*qnR * vR + pR*ny
  fR(4) = rhoR*qnR * wR + pR*nz
  fR(5) = rhoR*qnR * HR

!--------------------------------------------------------------------------------
!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
! Note: n1 and n2 may need to be frozen at some point during 
!       a steady calculation to fully make it converge. For time-accurate 
!       calculation, it is not necessary.
! Note: For a boundary face, you may want to set (nx2,ny2,nz2)=(nx,ny,nz),
!       (nx1,ny1)=tangent vector, i.e., use the Roe flux.

    abs_dq = sqrt( (uR-uL)**2 + (vR-vL)**2 + (wR-wL)**2 )

!n1 = Velocity difference vector: normal to shock or tangent to shear
  if ( abs_dq > 1.0e-12_p2) then

       nx1 = (uR-uL)/abs_dq
       ny1 = (vR-vL)/abs_dq
       nz1 = (wR-wL)/abs_dq

!n1 = Face tangent vector if abs_dq is too small.
!     Note: There are infinitely many choices for the tangent vector.
!           The best choice may be discovered in future.
! Here, we choose a vector in a plane (essentially 2D vector).
! Note that we must be careful and make sure that the vector is not a zero vector.
  else

!  E.g., if n = (1,0,0), then, (0,-nz,ny) = (0,0,0). This is a problem.
!        To avoid such a situation, we choose the 2D tangential vector
!        having the maximum length.

     tempx = ny*ny + nz*nz
     tempy = nz*nz + nx*nx
     tempz = nx*nx + ny*ny

     if     ( tempx >= tempy .and. tempx >= tempz ) then
       nx1 =  zero
       ny1 = -nz
       nz1 =  ny
     elseif ( tempy >= tempx .and. tempy >= tempz ) then
       nx1 = -nz
       ny1 =  zero
       nz1 =  nx
     elseif ( tempz >= tempx .and. tempz >= tempy ) then
       nx1 = -ny
       ny1 =  nx
       nz1 =  zero
     else
      ! Impossible to happen
      write(*,*) "inviscid_rotated_rhll: Impossible to happen. Please report the problem."
      stop
     endif

!     Make it the unit vector.
      temp = sqrt( nx1*nx1 + ny1*ny1 + nz1*nz1 )
       nx1 = nx1/temp
       ny1 = ny1/temp
       nz1 = nz1/temp

  endif

    alpha1 = nx*nx1 + ny*ny1 + nz*nz1

! Make alpha1 always positive.
      temp = sign(one,alpha1)
       nx1 = temp * nx1
       ny1 = temp * ny1
       nz1 = temp * nz1
    alpha1 = temp * alpha1

!n2 = direction perpendicular to n1.
!     Note: There are infinitely many choices for this vector.
!           The best choice may be discovered in future.
! Here, we employ the formula (4.4) in the paper:
!     (nx2,ny2,nz2) = (n1xn)xn1 / |(n1xn)xn1|    ('x' is the vector product.)

!  (tempx,tempy,tempz) = n1xn
     tempx = ny1*nz - nz1*ny
     tempy = nz1*nx - nx1*nz
     tempz = nx1*ny - ny1*nx

!  (nx2,ny2,nz2) = (n1xn)xn1
     nx2 = tempy*nz1 - tempz*ny1
     ny2 = tempz*nx1 - tempx*nz1
     nz2 = tempx*ny1 - tempy*nx1

!  Make n2 the unit vector
     temp = sqrt( nx2*nx2 + ny2*ny2 + nz2*nz2 )
       nx2 = nx2/temp
       ny2 = ny2/temp
       nz2 = nz2/temp

    alpha2 = nx*nx2 + ny*ny2 + nz*nz2

!  Make alpha2 always positive.
      temp = sign(one,alpha2)
       nx2 = temp * nx2
       ny2 = temp * ny2
       nz2 = temp * nz2
    alpha2 = temp * alpha2

!--------------------------------------------------------------------------------
!Now we are going to compute the Roe flux with n2 as the normal with modified 
!wave speeds (5.12). NOTE: the Roe flux here is computed without tangent vectors.
!See "I do like CFD, VOL.1" for details: page 57, Equation (3.6.31).

!First compute the Roe-averaged quantities

!  NOTE: See http://www.cfdnotes.com/cfdnotes_roe_averaged_density.html for
!        the Roe-averaged density.

      RT = sqrt(rhoR/rhoL)
     rho = RT*rhoL                                        !Roe-averaged density.
       u = (uL + RT*uR)/(one + RT)                        !Roe-averaged x-velocity
       v = (vL + RT*vR)/(one + RT)                        !Roe-averaged y-velocity
       w = (wL + RT*wR)/(one + RT)                        !Roe-averaged z-velocity
       H = (HL + RT*HR)/(one + RT)                        !Roe-averaged total enthalpy
       a = sqrt( (gamma-one)*(H-half*(u*u + v*v + w*w)) ) !Roe-averaged speed of sound

!----------------------------------------------------
!Compute the wave speed estimates for the HLL part,
!following Einfeldt:
!
! B. Einfeldt, On Godunov-type methods for gas dynamics,
! SIAM Journal on Numerical Analysis 25 (2) (1988) 294–318.
!
! Note: HLL is actually applied to n1, but this is
!       all we need to incorporate HLL. See JCP2008 paper.

     qn  = u *nx1 + v *ny1 + w *nz1
     qnL = uL*nx1 + vL*ny1 + wL*nz1
     qnR = uR*nx1 + vR*ny1 + wR*nz1
     SLm = min( zero, qn - a, qnL - aL ) !Minimum wave speed estimate
     SRp = max( zero, qn + a, qnR + aR ) !Maximum wave speed estimate

! This is the only place where n1=(nx1,ny1,nz1) is used.
! n1=(nx1,ny1,nz1) is never used below.
!----------------------------------------------------

!Wave Strengths

     qn  = u *nx2 + v *ny2 + w *nz2
     qnL = uL*nx2 + vL*ny2 + wL*nz2
     qnR = uR*nx2 + vR*ny2 + wR*nz2

    drho = rhoR - rhoL  !Density difference
      dp =   pR - pL    !Pressure difference
     dqn =  qnR - qnL   !Normal velocity difference

  LdU(1) = (dp - rho*a*dqn )/(two*a*a) !Left-moving acoustic wave strength
  LdU(2) =  drho - dp/(a*a)            !Entropy wave strength
  LdU(3) = (dp + rho*a*dqn )/(two*a*a) !Right-moving acoustic wave strength
  LdU(4) = rho                         !Shear wave strength (not really, just a factor)

!Wave Speed (Eigenvalues)

  eig(1) = qn-a !Left-moving acoustic wave velocity
  eig(2) = qn   !Entropy wave velocity
  eig(3) = qn+a !Right-moving acoustic wave velocity
  eig(4) = qn   !Shear wave velocity

!Absolute values of the wave speeds (Eigenvalues)

   ws(1) = abs(qn-a) !Left-moving acoustic wave speed
   ws(2) = abs(qn)   !Entropy wave speed
   ws(3) = abs(qn+a) !Right-moving acoustic wave speed
   ws(4) = abs(qn)   !Shear wave speed

!Harten's Entropy Fix JCP(1983), 49, pp357-393: only for the nonlinear fields.
!NOTE: It avoids vanishing wave speeds by making a parabolic fit near ws = 0.

  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(3) = fifth
   if ( ws(3) < dws(3) ) ws(3) = half * ( ws(3)*ws(3)/dws(3)+dws(3) )

!Combine the wave speeds for Rotated-RHLL: Eq.(5.12) in the original JCP2008 paper.

      ws = alpha2*ws - (alpha1*two*SRp*SLm + alpha2*(SRp+SLm)*eig)/(SRp-SLm)

!Below, we compute the Roe dissipation term in the direction n2
!with the above modified wave speeds. HLL wave speeds act something like
!the entropy fix or eigenvalue limiting; they contribute only by the amount
!given by the fraction, alpha1 (less than or equal to 1.0). See JCP2008 paper.

!Right Eigenvectors:
!Note: Two shear wave components are combined into one, so that tangent vectors
!      are not required. And that's why there are only 4 vectors here.

! Left-moving acoustic wave
  R(1,1) = one    
  R(2,1) = u - a*nx2
  R(3,1) = v - a*ny2
  R(4,1) = w - a*nz2
  R(5,1) = H - a*qn

! Entropy wave
  R(1,2) = one
  R(2,2) = u
  R(3,2) = v 
  R(4,2) = w
  R(5,2) = half*(u*u + v*v + w*w)

! Right-moving acoustic wave
  R(1,3) = one
  R(2,3) = u + a*nx2
  R(3,3) = v + a*ny2
  R(4,3) = w + a*nz2
  R(5,3) = H + a*qn

! Two shear wave components combined into one (wave strength incorporated).
  du = uR - uL
  dv = vR - vL
  dw = wR - wL
  R(1,4) = zero
  R(2,4) = du - dqn*nx2
  R(3,4) = dv - dqn*ny2
  R(4,4) = dw - dqn*nz2
  R(5,4) = u*du + v*dv + w*dw - qn*dqn

!Dissipation Term: Roe dissipation with the modified wave speeds.
! |An|dU = R|Lambda|L*dU = sum_k of [ ws(k) * R(:,k) * L*dU(k) ], where n=n2.

 diss(:) = ws(1)*LdU(1)*R(:,1) + ws(2)*LdU(2)*R(:,2) &
         + ws(3)*LdU(3)*R(:,3) + ws(4)*LdU(4)*R(:,4)

!Compute the Rotated-RHLL flux. (It looks like the HLL flux with Roe dissipation.)

  num_flux = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss

!Normal max wave speed in the normal direction.
!  wsn = abs(qn) + a

 end subroutine inviscid_rotated_rhll
!--------------------------------------------------------------------------------
