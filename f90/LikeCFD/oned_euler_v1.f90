!********************************************************************************
!* One-dimensional Euler solver for a Sod's shock tube problem.
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* This is Version 1 (2010).
!* 
!* This F90 program was written and made available for download
!* for an educational purpose. Comments are welcome.
!*
!* This file may be updated in future.
!*
!* Katate Masatsuka, December 2010. http://www.cfdbooks.com
!********************************************************************************

!********************************************************************************
!* --- Main program for the 1D Euler shock-tube solver.
!*
!* This code solves the Sod's shock tube problem which is described
!* in Section 7.13.3 of "I do like CFD, VOL.1": Sod's problem 1, Figure 7.12.2.
!*
!* - t=0                               - t=tf
!* Density                             Density
!*   ****************|                 *********\
!*                   |                           \
!*                   |                            \
!*                   |                             ****|
!*                   |                                 |
!*                   |                                 ****|
!*                   ***************                       ***********
!*
!* Methods employed:
!*   - Roe's flux
!*   - Minmod limiter
!*   - Two-stage Runge-Kutta time stepping scheme
!*
!* Input ---------------------------------------------------
!*    ncells = # of cells on a grid.
!*        tf = Final time
!*       cfl = CFL number (<1)
!*      xmin = Left end of the domain
!*      xmax = Right end of the domain
!* 
!* Output --------------------------------------------------
!*  "solution.dat" = Data file containing for each cell,
!*                   cell-center coordinate, density, velocity, pressure, 
!*                   entropy, in the following format:
!*
!*       write(*,*) ncells
!*      do i=1,ncells
!*       write(*,*) xc(i), rho(i), u(i), p(i), entropy(i)
!*      end do
!* 
!*     Use the matlab program, oned_euler_v1.m, to plot the solutions.
!*
!*
!*  Note: Explore other shock tube problems by changing the initial condition
!*        and the final time (Other problems are described in Section 7.13.3
!*        of "I do like CFD, VOL.1").
!*
!*  Note: Other limiters may be explored (see CFD textbooks).
!*
!*  Note: Other flux functions may be explored.
!*        Various subroutines are available at cfdbooks.com: Osher, Van Leer, etc.
!*
!*  Note: Boundary condition need to be modified for problems with hard
!*        boundaries.
!*
!*
!* Katate Masatsuka, December 2010. http://www.cfdbooks.com
!*
!*
!* 12-29-10: Some compiler warnings fixed.
!*
!********************************************************************************
 program oned_euler

 implicit none
!Numeric parameters: [Note: Change 13 -> 5 to make everything single precision.]
   integer , parameter :: p2 = selected_real_kind(13) !Double Precision
   real(p2), parameter ::  zero = 0.0_p2
   real(p2), parameter ::   one = 1.0_p2
   real(p2), parameter ::  half = 0.5_p2
   real(p2), parameter :: gamma = 1.4_p2  !Ratio of specific heats for air

!Derived data type: these data are stored in each cell.
 type cell_data
  real(p2) :: xc     ! Cell-center coordinate
  real(p2) :: u(3)   ! Conservative variables = [rho, rho*u, rho*E]
  real(p2) :: u0(3)  ! Conservative variables at the previous time step
  real(p2) :: w(3)   ! Primitive variables = [rho, u, p]
  real(p2) :: dw(3)  ! Slope (difference) of primitive variables
  real(p2) :: res(3) ! Residual = f_{j+1/2) - f_{j-1/2)
 end type cell_data

!Local variables
 type(cell_data), allocatable :: cell(:)    !Array of cell-data
 real(p2)                     :: xmin, xmax !Left and right ends of the domain
 real(p2)                     :: dx         !Cell spacing (uniform grid)
 real(p2)                     :: t, tf      !Current time and final time
 real(p2)                     :: cfl, dt    !CFL number and global time step
 integer                      :: ncells     !Total number of cells
 integer                      :: nsteps     !Number of time steps
 integer                      :: itime      !Index for time stepping
 integer                      :: istage     !Index for Runge-Kutta stages
 integer                      :: i, j

!Local variables used for computing numerical fluxes.
 real(p2), dimension(3) :: dwl, dwr !Slopes between j and j-1, j and j+1
 real(p2), dimension(3) :: wL, wR   !Extrapolated states at a face
 real(p2), dimension(3) :: flux     !Numerical flux

!--------------------------------------------------------------------------------
! 0. Input parameters and initial condition.

!Parameters
  ncells = 80     ! Number of cells
      tf = 1.7_p2 ! Finial time
     cfl = 0.8_p2 ! CFL number
    xmin =-5.0_p2 ! Left boundary coordinate
    xmax = 5.0_p2 ! Right boundary coordinate

! Allocate the cell array: 2 ghost cells, 0 on the left and ncells+1 on the right.
!  E.g., data in cell j is accessed by cell(j)%xc, cell(j)%u, cell(j)%w, etc.
  allocate(cell(0:ncells+1))

! Cell spacing (grid is uniform)
  dx = (xmax-xmin)/real(ncells)

! The initial condition for Sod's shock tube problem (I Do Like CFD, VOL.1, page 199).
! [Note: Change these data (and tf) to solve different problems.]

   do i = 0, ncells+1 !Include ghost cells at i=0 and ncells+1.

      if (i <= ncells/2) then
       cell(i)%w(1) = 1.0_p2   !Density  on the left
       cell(i)%w(2) = 0.0_p2   !Velocity on the left
       cell(i)%w(3) = 1.0_p2   !Pressure on the left
      else
       cell(i)%w(1) = 0.125_p2 !Density  on the right
       cell(i)%w(2) = 0.0_p2   !Velocity on the right
       cell(i)%w(3) = 0.1_p2   !Pressure on the right
      endif

      cell(i)%u = w2u(cell(i)%w,gamma) !Compute the conservative variables
      cell(i)%xc=xmin+real(i-1)*dx     !Cell center coordinate

   end do

!--------------------------------------------------------------------------------
! Time stepping loop to reach t = tf 
!--------------------------------------------------------------------------------

      t = zero !Initialize the current time.
 nsteps = 0    !Initialize the number of time steps.
 time_step : do itime = 1 , 50000 !50000 is large enough to reach tf=1.7.

   if (t==tf) exit                    !Finish if the final time is reached.
   dt = timestep(cfl,dx,gamma,ncells) !Compute the global time step.
   if (t+dt > tf) dt =  tf - t        !Adjust dt to finish exactly at t=tf.
        t = t + dt                    !Update the current time.
   nsteps = nsteps + 1                !Count the number of time steps.

!---------------------------------------------------
! Runge-Kutta Stages
!
! Two-stage Runge-Kutta scheme:
!  1. u^*     = u^n - dt/dx*Res(u^n)
!  2. u^{n+1} = 1/2*u^n + 1/2*[u^*- dt/dx*Res(u^*)]
!---------------------------------------------------
  rk_stages : do istage = 1, 2

!(1) Residual computation: compute cell(:)%res(1:3).

! Compute the slopes (as difference) at every cell.
! NB: for uniform meshes, difference (du) can be used in place of gradient (du/dx).

   reconstruction : do j = 1, ncells

    dwl = cell(j  )%w-cell(j-1)%w !Simple central-difference between j   and j-1.
    dwr = cell(j+1)%w-cell(j  )%w !Simple central-difference between j+1 and j.

!    Apply a slope limiter.
!    (minmod: zero if opposite sign, otherwise the one of smaller magnitude.)
     do i = 1, 3
      cell(j)%dw(i) = minmod(dwl(i),dwr(i))
     end do

   end do reconstruction

! Initialize the residuals.
   initialize_res : do j = 1, ncells
    cell(j)%res = zero
   end do initialize_res

! Compute the residuals: residual_j = flux_{j+1/2} - flux_{j-1/2}.
! Here, compute the flux at j+1/2 and add it to the left cell and subtract
! from the right cell. Only the internal faces are considered; the left
! and right most faces are considered later.

!     j+1/2
!   | wL|   |
!   |  /|wR |
!   | / |\  |
!   |/  | \ |
!   |   |  \|
!   |   |   |
!     j  j+1
!
   flux_comp : do j = 1, ncells-1
      wL = cell(j  )%w + half*cell(j  )%dw !State extrapolated to j+1/2 from j
      wR = cell(j+1)%w - half*cell(j+1)%dw !State extrapolated to j+1/2 from j+1
    flux = roe_flux(wL,wR,gamma)           !Numerical flux at j+1/2
    cell(j  )%res = cell(j  )%res + flux   !Add it to the left cell.
    cell(j+1)%res = cell(j+1)%res - flux   !Subtract from the right cell.
   end do flux_comp

! Add boundary fluxes: left end and right end.
! For the problem considered here, it suffices to simply copy the state
!  from inside the domain to the ghost cell (no gradient condition).

!  Left most face: left face of cell i=1.
   wR = cell(2)%w - half*cell(2)%dw  !State extrapolated to j-1/2 from j=1
   wL = wR                           !The same state
   flux = roe_flux(wL,wR,gamma)      !Use Roe flux to compute the flux.
   cell(1)%res = cell(1)%res - flux  !Subtract the flux: -flux_{j-1/2}.

!  Right most face: right face of cell i=ncells.
   wL = cell(ncells)%w + half*cell(ncells)%dw !State extrapolated to ncells+1/2 from j=ncells
   wR = wL                                    !The same state
   flux = roe_flux(wL,wR,gamma)               !Use Roe flux to compute the flux.
   cell(ncells)%res = cell(ncells)%res + flux !Add the flux: +flux_{j+1/2}.

!(2) Solution update

  if (istage==1) then
!  1st Stage of Runge-Kutta: save u^n as u0(:); u^* is stored at u(:).
   stage01_update : do j = 1, ncells
    cell(j)%u0 = cell(j)%u            !Save the solution at n for 2nd stage.
    cell(j)%u  = cell(j)%u - (dt/dx)*cell(j)%res
    cell(j)%w  = u2w(cell(j)%u,gamma) !Update primitive variables
   end do stage01_update

  else
!  2nd Stage of Runge-Kutta:
   stage02_update : do j = 1, ncells
    cell(j)%u = cell(j)%u - (dt/dx)*cell(j)%res
    cell(j)%u = half*(cell(j)%u0 + cell(j)%u )
    cell(j)%w = u2w(cell(j)%u,gamma)  !Update primitive variables
   end do stage02_update

  endif

! Copy the solutions to the ghost cells.
! In this program, the ghost cell values are used only in the reconstruction.
  cell(0)%w        = cell(1)%w
  cell(ncells+1)%w = cell(ncells)%w

  end do rk_stages
!---------------------------------------------------
! End of Runge-Kutta Stages
!---------------------------------------------------

 end do time_step
!--------------------------------------------------------------------------------
! End of time stepping
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! Write out a data file.

  write(*,*)
  write(*,*) "final time t(sec) = ",t," by ",nsteps," time steps"
  call output(gamma, ncells)

 stop
!********************************************************************************
! End of program
!********************************************************************************



!Below is a list of subroutines and functions used in the main program.

 contains

!*******************************************************************************
!* Compute the global time step dt restricted by a given CFL number.
!*
!* ------------------------------------------------------------------------------
!*  Input: CFL number
!* Output: Global time step, dt.
!* ------------------------------------------------------------------------------
!*
!*******************************************************************************
 function timestep(cfl,dx,gamma,ncells) result(dt)

 implicit none
 real(p2), intent(in) :: cfl, dx, gamma !Input
 integer , intent(in) :: ncells         !Input
 real(p2)             ::  dt            !Output
!Local variables
 real(p2) :: u, c, max_speed
 integer  :: i

  max_speed = -one

  do i = 1, ncells
   u = cell(i)%w(2)                          !Velocity
   c = sqrt(gamma*cell(i)%w(3)/cell(i)%w(1)) !Speed of sound
   max_speed = max( max_speed, abs(u)+c )
  end do

   dt = cfl*dx/max_speed !CFL condition: dt = CFL*dx/max_wavespeed, CFL <= 1.

 end function timestep

!****************************************************************************
!* Minmod limiter
!* --------------------------------------------------------------------------
!*  Input: two real values, a and b
!* Output: minmod of a and b.
!* --------------------------------------------------------------------------
!* 
!****************************************************************************
 function minmod(a,b)

 implicit none
 real(p2), intent(in) :: a, b   !Input
 real(p2)             :: minmod !Output

!Local parameter
 real(p2), parameter  :: zero = 0.0_p2

  if (a*b <= zero) then
    minmod = zero               !a>0 and b<0; or a<0 and b>0
  elseif (abs(a)<abs(b)) then
    minmod = a                  !|a| < |b|
  elseif (abs(b)<abs(a)) then
    minmod = b                  !|a| > |b|
  else
    minmod = a                  !Here, a=b, so just take a or b.
  endif

 end function minmod
!--------------------------------------------------------------------------------


!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  w = primitive variables (rho, u, p, 0, 0)
!* Output:  u = conservative variables (rho, rho*u, rho*E, 0, 0)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
 function w2u(w,gamma) result(u)
 implicit none
 real(p2), intent(in) :: w(3), gamma !Input
 real(p2)             :: u(3)        !output

  u(1) = w(1)
  u(2) = w(1)*w(2)
  u(3) = w(3)/(gamma-one)+half*w(1)*w(2)*w(2)

 end function w2u
!--------------------------------------------------------------------------------

!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  u = conservative variables (rho, rho*u, rho*E, 0, 0)
!* Output:  w = primitive variables (rho, u, p, 0, 0)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
 function u2w(u,gamma) result(w)
 implicit none
 real(p2), intent(in) :: u(3), gamma !Input
 real(p2)             :: w(3)        !output

  w(1) = u(1)
  w(2) = u(2)/u(1)
  w(3) = (gamma-one)*( u(3) - half*w(1)*w(2)*w(2) )

 end function u2w
!--------------------------------------------------------------------------------


!********************************************************************************
!* -- Roe's Flux Function without entropy fix---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* ------------------------------------------------------------------------------
!*  Input:   wL(1:3) =  left state (rhoL, uL, pL)
!*           wR(1:3) = right state (rhoR, uR, pR)
!*
!* Output:  flux(1:3) = numerical flux for the Euler equations (the Roe flux)
!* ------------------------------------------------------------------------------
!* 
!* Katate Masatsuka, December 2010. http://www.cfdbooks.com
!********************************************************************************
 function roe_flux(wL,wR,gamma) result(flux)

 implicit none
 real(p2), intent(in) :: wL(3), wR(3), gamma !  Input (conservative variables rho*[1, v, E])
 real(p2)             :: flux(3)             ! Output (numerical flux across L and R states)

!Local parameters
 real(p2), parameter ::    zero = 0.0_p2
 real(p2), parameter ::     one = 1.0_p2
 real(p2), parameter ::    four = 4.0_p2
 real(p2), parameter ::    half = 0.5_p2
 real(p2), parameter :: quarter = 0.25_p2
!Local variables
 real(p2) :: uL(3), uR(3)
 real(p2) :: rhoL, rhoR, vL, vR, pL, pR   ! Primitive variables.
 real(p2) :: aL, aR, HL, HR               ! Speeds of sound.
 real(p2) :: RT,rho,v,H,a                 ! Roe-averages
 real(p2) :: drho,du,dP,dV(3)
 real(p2) :: ws(3),Da, R(3,3)
 integer :: j, k

    uL = w2u(wL,gamma)
    uR = w2u(wR,gamma)

!Primitive and other variables.
!  Left state
    rhoL = wL(1)
      vL = wL(2)
      pL = wL(3)
      aL = sqrt(gamma*pL/rhoL)
      HL = ( uL(3) + pL ) / rhoL
!  Right state
    rhoR = wR(1)
      vR = wR(2)
      pR = wR(3)
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

!Modified wave speeds for nonlinear fields (the so-called entropy fix, which
!is often implemented to remove non-physical expansion shocks).
!There are various ways to implement the entropy fix. This is just one
!example. Try turn this off. The solution may be more accurate.
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
   flux = half*( euler_physical_flux(wL) + euler_physical_flux(wR) )

!Add the matrix dissipation term to complete the Roe flux.
  do j = 1, 3
   do k = 1, 3
    flux(j) = flux(j) - half*ws(k)*dV(k)*R(j,k) 
   end do
  end do

 end function roe_flux
!--------------------------------------------------------------------------------

!********************************************************************************
!* Physical flux of the Euler equations (inviscid part only).
!*
!* ------------------------------------------------------------------------------
!*  Input:     w(1:5) = primitive variables (rho, u, p, 0, 0)
!*
!* Output:  flux(1:5) = physical inviscid flux (rho, rho*u, rho*H, 0, 0)
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 function euler_physical_flux(w) result(flux)

 implicit none
 real(p2), dimension(3), intent(in) :: w    !Input
 real(p2), dimension(3)             :: flux !Output

!Local parameters
 real(p2), parameter :: half = 0.5_p2
!Local variables
 real(p2) :: rho, u, p
 real(p2) ::  a2

  rho = w(1)
    u = w(2)
    p = w(3)

   a2 = gamma*p/rho

  flux(1) = rho*u
  flux(2) = rho*u*u + p
  flux(3) = rho*u*( a2/(gamma-one) + half*u*u ) ! H = a2/(gamma-one) + half*u*u

 end function euler_physical_flux
!--------------------------------------------------------------------------------

!********************************************************************************
!* Write output data file
!*
!* ------------------------------------------------------------------------------
!* Output:  Data file "solution.dat" containing for each cell the following:
!*          cell-center coordinate, density, velocity, pressure, entropy
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine output(gamma,ncells)

 implicit none
 real(p2), intent(in) :: gamma  !Input
 integer , intent(in) :: ncells !Input

!Local parameters
 real(p2), parameter :: one = 1.0_p2
!Local variables
 real(p2) :: entropy
 integer  :: i, os

  open(unit=1, file = "solution.dat", status="unknown", iostat=os)

  do i = 1, ncells
   entropy = log(cell(i)%w(3)*cell(i)%w(1)**(-gamma))/(gamma-one) 
   write(1,'(5es25.15)') cell(i)%xc, cell(i)%w(1), cell(i)%w(2), cell(i)%w(3), entropy
  end do 

  close(1)

 end subroutine output
!--------------------------------------------------------------------------------

 end program oned_euler
!********************************************************************************









