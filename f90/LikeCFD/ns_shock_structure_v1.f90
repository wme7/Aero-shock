!********************************************************************************
!* Program for computing the steady shock structure based on the Navier-Stokes
!* equations. (for use in the code verification of Navier-Stokes codes)
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* This is Version 1 (2010).
!* 
!* These F90 routines were written and made available for download
!* for an educational purpose. Also, this is intended to provide all CFD
!* students and researchers with a tool for their code verification.
!*
!* This file may be updated in future.
!*
!*
!* 12-29-10: Minor bugs fixed.
!*
!*
!* Katate Masatsuka, September 2010. http://www.cfdbooks.com
!********************************************************************************

!********************************************************************************
!* --- Main program for the NS shock structure.
!*
!*   uL(>1) **********     Shock width depends on the Renolds number.
!*                    *
!*                     *
!*                      *********** uR
!*
!* This code solves the 1D Navier-Stokes equations across a steady shock.
!* The steady shock is defined in Section 7.11.8 of "I do like CFD, VOL.1".
!*
!* Input ---------------------------------------------------
!*    nnodes = # of nodes on a grid.
!*     gamma = Ratio of specific heats (1.4 for air)
!*       Pra = Prandtl number
!*     M_inf = Upstream Mach number
!*    Re_inf = Upstream Reynolds number
!*     T_inf = Upstream temperature (K)
!*      xmin = Left end of the domain
!*      xmax = Right end of the domain
!*       eps = Perturbation to the initial velocity.
!* 
!* Output --------------------------------------------------
!*  "ns_shock_structure.dat" = a data file containing grid coordinates and
!*                             the exact solution values in the following format:
!*
!*       write(*,*) nnodes
!*      do i=1,nnodes
!*       write(*,*) x(i), rho(i), u(i), p(i), tau(i), q(i)
!*      end do
!* 
!*     Read this file in your Navier-Stokes code and start running the code.
!*     (Use the exact solution as an initial solution for your calculation;
!*     fix the pressure somewhere in the middle of the shock as an internal
!*     boundary condition.) Can your Navier-Stokes code compute the shock 
!*     structure accurately?
!*
!*     Use the matlab program, ns_shock_structure_plot_v1.m, to plot the 
!*     exact solution.
!*
!*         x = node coordinate
!*       rho = density
!*         u = velocity
!*         p = pressure
!*       tau = viscous stress
!*         q = heat flux
!*
!*  Note: the solutions are nondimensionalized values (compressible scale),
!*          rho=rho/rho_inf, rho=u/a_inf, rho=p/(rho_inf*a_inf^2), 
!*            T=T/T_inf    ,  mu=mu/mu_inf.
!*
!*  Note: the shock location depends on 'eps'. If eps=0, nothting will happen and
!*        the computation of the viscous shock fails. The shock appear sooner
!*        for larger value of eps. 
!*
!*  Note: The exact solutions are point values. For high-order (>2nd) cell-centered
!*        codes (that compute cell-averages, you need to modify this code to output
!*        cell-averages.)
!*
!* Katate Masatsuka, September 2010. http://www.cfdbooks.com
!********************************************************************************
 program ns_shock_structure_main

 implicit none

!Numeric parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter :: zero=0.0_dp, one=1.0_dp, two=2.0_dp, three=3.0_dp
 real(dp), parameter :: four=4.0_dp
 real(dp), parameter :: half=0.5_dp

!Local Variables
 integer  :: nnodes ! Number of nodes in the grid
 real(dp) ::  gamma ! Ratio of specific heats
 real(dp) ::    Pra ! Prandtl number
 real(dp) ::  M_inf ! Upstream Mach number
 real(dp) :: Re_inf ! Upstream Reynolds number
 real(dp) ::  T_inf ! Upstream temperature (K)
 real(dp) ::   xmin ! Left end of the domain
 real(dp) ::   xmax ! Right end of the domain
 real(dp) :: eps    ! Perturbation to the initial velocity.

 integer  :: i
 real(dp) :: h                !Grid spacing of a uniform grid
 real(dp) :: ML               !Mach number of the left state (=M_inf)
 real(dp) :: rhoL, uL, pL, TL ! Left state of the shock
 real(dp) :: rhoR, uR, pR, TR !Right state of the shock
 real(dp) :: C(3)             !Rankine-Hugoniot constants
 real(dp) :: u0, T0           !Initial solutions for the ODE system
 real(dp) :: ui, Ti           !Computed solutions at node i
 real(dp) :: x0, xi           !Initial location, final location to solve the ODE
                              !X-coordinates, primitive variables
 real(dp), dimension(:), allocatable :: x, rho, u, p, T, tau, q

!--------------------------------------------------------------------------------
! 0. Input parameters

  nnodes =  81         ! Number of nodes in the grid
   gamma =  1.4_dp     ! Ratio of specific heats (1.4 for air)
     Pra =  three/four ! Prandtl number
   M_inf =  3.5_dp     ! Upstream Mach number
  Re_inf =  25.0_dp    ! Upstream Reynolds number
   T_inf =  400.0_dp   ! Upstream temperature (K)
    xmin = -1.0_dp     ! Left end of the domain
    xmax =  1.0_dp     ! Right end of the domain
     eps = 1.0e-10_dp  ! Perturbation to the initial velocity.

 allocate(x(nnodes), rho(nnodes), u(nnodes), p(nnodes))
 allocate(T(nnodes), tau(nnodes), q(nnodes))
!--------------------------------------------------------------------------------
! 1. Generate a grid

  write(*,*) "Generating a grid ......."

   h = (xmax-xmin)/real(nnodes-1,dp) !Domain length divided by ncells=nnodes-1.

! Uniform grid
  do i = 1, nnodes
   x(i) = xmin + h*real(i-1,dp)
  end do

  write(*,*) "Finished generating a grid."

!--------------------------------------------------------------------------------
! 2. Set up the left and right states of the shock (Rankine-Hugoniot).
!    (cf. Section 7.11.8 of "I do like CFD, VOL.1".)

  write(*,*) "Setting up the states and constants ......."

    ML = M_inf

! Left state
  rhoL = one
    uL = ML
    pL = one/gamma
    TL = gamma*pL/rhoL

! Right state
  rhoR = (gamma+one)*ML*ML/( (gamma-one)*ML*ML + two )
    uR = ML/rhoR
    pR = (one+two*gamma/(gamma+one)*(ML*ML-one))/gamma
    TR = gamma*pR/rhoR

!Rankine-Hugoniot constants
 write(*,*) "Conservation check"
 write(*,*) " Left state:"
   C(1) = rhoL*uL                                       !Mass flux
   C(2) = rhoL*uL*uL + pL                               !Momentum flux
   C(3) = uL*( gamma*pL/(gamma-one) + half*rhoL*uL*uL ) !Energy flux
 write(*,'(3es25.15)') C(1), C(2), C(3)
 write(*,*) " Right state:"
 write(*,'(3es25.15)') C(1), C(2), C(3)
   C(1) = rhoR*uR                                       !Mass flux
   C(2) = rhoR*uR*uR + pR                               !Momentum flux
   C(3) = uR*( gamma*pR/(gamma-one) + half*rhoR*uR*uR ) !Energy flux

  write(*,*) "Finished setting up the states and constants ......."

!--------------------------------------------------------------------------------
! 3. Set initial values for the ODE: the right end of the domain (the right end
!    of the grid).

   rho(nnodes) = rhoR
     u(nnodes) =   uR
     p(nnodes) =   pR
     T(nnodes) =   TR
   tau(nnodes) = zero
     q(nnodes) = zero

!--------------------------------------------------------------------------------
! 4. Loop over nodes: compute the exact solution at each node
!    by integrating the ODE for u and T from the previous node.

  write(*,*) "Solving the ODE from one node to the next......."
  node_loop: do i = nnodes-1, 1, -1

! Previous node location and the solution there.
   x0 = x(i+1)
   u0 = u(i+1)
   T0 = T(i+1)
! The current node location where we compute the solution.
   xi = x(i)

! Perturbation to the initial velocity. This is required to generate
! the viscous shock. Without this, nothing happens and you'll never
! reach the left state.
  if (i == nnodes-1) u0 = u0*(one+eps)

! Integrate ODE from x0 to xi.
   call ode_solver_rk4(x0,u0,T0,  xi,ui,Ti)
   write(*,'(i5,3es12.3,3x,a4,es10.3)') i, xi, ui, Ti, " TL=",TL

! Store the computed solution at the current node (x=xi).
    rho(i) = C(1)/ui
      u(i) = ui
      p(i) = rho(i)*Ti/gamma
      T(i) = Ti
    tau(i) = ( C(1)*ui + (C(1)/ui)*(Ti/gamma) - C(2) )                     
      q(i) = ( C(1)*ui*ui*half-C(1)/(gamma-one)*(Ti/gamma)+C(3)-C(2)*ui )

  end do node_loop
  write(*,*) "Finished solving the ODE......."

!--------------------------------------------------------------------------------
! 5. Write a data file: "ns_shock_structure.dat"

  write(*,*) "Writing a data file......."
  write(*,*) " filename = ns_shock_structure.dat"
  open(unit=10, file ="ns_shock_structure.dat", status="unknown")
    write(10,*) nnodes
   do i = 1, nnodes
    write(10,'(7es28.15)') x(i),rho(i),u(i),p(i),tau(i),q(i),T(i)
   end do
  write(*,*) "Finish writing a data file......."
  close(10)

!--------------------------------------------------------------------------------
!
  write(*,*) "Exact solutions computed."
  write(*,*) "Use the matlab file ns_shock_structure.m to make plots."

 stop
!********************************************************************************

 contains

!********************************************************************************
!* This subroutine integrates the ODE from x=x0 to x=xi and find the velocity and
!* temperature at x=xi: ui and Ti.
!*
!* Input ---------------------------------------------------
!*   x0 =    initial position
!*   u0 =    velocity at x=x0
!*   T0 = temperature at x=x0
!*
!* Output --------------------------------------------------
!*   xi =      final position
!*   ui =    velocity at x=xi
!*   Ti = temperature at x=xi
!*
!* Note: the ODE is integrated from x=x0 to x=xi in 1000 steps by the classical
!*       4th-order Runge-Kutta method.
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* Katate Masatsuka, September 2010. http://www.cfdbooks.com
!********************************************************************************
 subroutine ode_solver_rk4(x0,u0,T0,  xi,ui,Ti)
 implicit none
!Parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter :: two=2.0_dp
 real(dp), parameter :: third=1.0_dp/3.0_dp, sixth=1.0_dp/6.0_dp, half=0.5_dp
 real(dp), parameter :: thousand=1000.0_dp
!Input and Output
 real(dp), intent( in)  :: x0,u0,T0 !Initial position, velocity, and temperature
 real(dp), intent( in)  :: xi       !Final position
 real(dp), intent(out)  ::    ui,Ti !velocity and temperature at xi
!Local variables
 real(dp)               :: dz     !Step size for ODE integration
 real(dp)               :: z0, zi !Initial and final locations
 real(dp)               ::  z     !Current location
 real(dp), dimension(2) :: V      !Solution vector = (velocity,temperature)
 real(dp), dimension(2) :: K1, K2, K3 !Vector used in Runge-Kutta
 logical                :: finish

! Increment for ODE integration: 1000 steps between x0 and xi.

     dz = (x0-xi)/thousand

! Integrate ODE, dV/dx=RHS, from x0 to xi by the classical RK4.
!      xi <--- x0
!  -----.-------.-------->x
!
! We integrate actually dV/dz=-RHS where z=-x.
!      z0 ---> zi
!  -----.-------.-------->z
!
! 1. Initial values and final location, zi, in the reversed coordinate.

      z0 = -x0
    V(1) =  u0
    V(2) =  T0
      zi = -xi

! 2. Stepping from -x0 to -xi. Solve the ODE with x=-x.

      z = z0
      finish = .false.

    z0_to_zi: do

!    To finish up at z=zi precisely
     if (z + dz > zi) then
      dz = zi - z
      finish = .true.
     endif

       z = z + dz
      K1 = V + half*dz*( rhs( V) )
      K2 = V + half*dz*( rhs(K1) )
      K3 = V +      dz*( rhs(K2) )
       V = (K1 + two*K2 + K3 - V)*third + dz*( rhs(K3) )*sixth

     if (finish) exit

    end do z0_to_zi

! Solution at z=zi, i.e., x=xi.

   ui =  V(1)
   Ti =  V(2)

 return
 end subroutine ode_solver_rk4
!*******************************************************************************

!*******************************************************************************
!* Right hand side of the original ODE: dV/dx = RHS (not with respect to z)
!*
!* This function is used in the subroutine, ode_solver_rk4().
!*
!* Katate Masatsuka, September 2010. http://www.cfdbooks.com
!*******************************************************************************
 function rhs(V) result(right_hand_side)
 implicit none
 integer , parameter  :: sp = kind(1.0)
 integer , parameter  :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter  :: one=1.0_dp,four=4.0_dp
 real(dp), parameter  :: half=0.5_dp,third=1.0_dp/3.0_dp

!Input and output
 real(dp), dimension(2), intent(in) :: V               !input
 real(dp), dimension(2)             :: right_hand_side !Output
!Local variables
 real(dp) :: Cu, CT     !Viscous stress and heat flux
 real(dp) ::  u,  T, mu !velocity, temperature, viscosity

    u = V(1)
    T = V(2)
   mu = viscosity(T)
                                    ! Nondimensionalization factor(compressible scale)
   Cu =  four*third*mu              * M_inf/Re_inf
   CT = -gamma*mu/(Pra*(gamma-one)) * M_inf/Re_inf/gamma

  right_hand_side(1) =-( C(1)*u + C(1)*T/gamma/u - C(2) )                    /Cu
  right_hand_side(2) =-( C(1)*u*u*half-C(1)/(gamma-one)*T/gamma+C(3)-C(2)*u )/CT

 end function rhs
!*******************************************************************************

!*******************************************************************************
!* Viscousity: Sutherland's law
!* NB: You may replace this by any other law according to the formula you use
!*     in your Navier-Stokes code (e.g., a power law).
!*******************************************************************************
 function viscosity(T) result(mu)
 implicit none
 integer , parameter  :: sp = kind(1.0)
 integer , parameter  :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter  :: zero=0.0_dp,one=1.0_dp,three_half=3.0_dp/2.0_dp

!Input and output
 real(dp), intent(in) :: T  ! input
 real(dp)             :: mu !output
!Parameter for Sutherland's law:
 real(dp), parameter  :: SC = 110.5_dp

 if (T < zero) then
  write(*,*) " Negative temperature entered viscosity(): T = ", T
  write(*,*) " Stop"
  stop
 endif

! Sutherland's law (nondimensional form)
  mu = (one+SC/T_inf)/(T+SC/T_inf)*T**three_half

 end function viscosity
!*******************************************************************************

 end program ns_shock_structure_main
