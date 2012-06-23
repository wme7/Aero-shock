!********************************************************************************
!* Subroutine for computing the Blasius solution for a flow over a flat plate.
!* (often used for code verification of Navier-Stokes codes)
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
!* 12-25-10: Minor bugs fixed.
!*
!* Katate Masatsuka, March 2010. http://www.cfdbooks.com
!********************************************************************************

!********************************************************************************
!* --- Driver code for Blasius boundary layer solution ---
!*
!*                  Exact solutions are computed along this vertical line.
!*
!*                                                   |
!*  ^ y                                              |
!*  |     U_inf --->                                 .(xp,yp)
!*  |                                                |
!*  ----->x         ------------------------------------------------------------
!*                x=0                  flat plate
!*
!* This driver code calls the subroutine that computes the Blasius solution
!* at a given position (xp,yp) of the flow described above, and generate a set of
!* solution along a vertical line at x=xp where xp is fixed.
!*
!* Input ---------------------------------------------------
!*        Np = # of points along the vertical line.
!*   (xp,yp) = points along a line: xp is fixed but yp is a set of values.
!*        nu = kinematic viscosity
!*     U_inf = Freestream velocity (x-direction only)
!* 
!* Output --------------------------------------------------
!*  blasius_solutions.dat: a data file containing 
!*                              etap,u,v,f,f1,f2, xp,yp,nu,U_inf
!*                         along the vertical line, i.e., Np rows.
!* 
!* Katate Masatsuka, March 2010. http://www.cfdbooks.com
!********************************************************************************
 program blasius_main
 implicit none
!numeric parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter :: one=1.0_dp, ten=10.0_dp
!local variables
 integer  :: j, ny
 real(dp) :: u,v,f,f1,f2,etap, xp,yp
!flow parameters
 real(dp) :: U_inf = one        ! Freesteam velocity (x-velocity only)
 real(dp) ::    nu = 1.0e-02_dp ! kinematic viscosity
 real(dp) ::    dy = 1.0e-01_dp ! Spacing of the data points in the line.

!-------------------------------------------------------------------------------
! Compute the exact solution along the vertical line at x=xp.

  open(unit=10, file ="blasius_solutions.dat", status="unknown")
  write(*,*) "Computing exact solutions ......."
  xp = ten ! x=xp: fixed, defining the location of the vertical line.
  ny = 50  ! # of points along the line.

!--------------------------------------------------------------
!Loop over points along the line: from y=0 to y=ymax=dy*(ny-1).
  do j = 1, ny

! y-coordinate of the next point in the line.
   yp = dy*real(j-1)
! Compute the exact solution at (xp,yp): integrate the ODE.
   call blasius(etap,u,v,f,f1,f2, xp,yp,nu,U_inf)
! Write the exact solution to the file.
   write(10,'(10es25.15)') xp,yp, etap,u,v,f,f1,f2, nu,U_inf

  end do
!--------------------------------------------------------------

  write(*,*) " Rex = ", U_inf*xp/nu
  write(*,*) "Exact solutions computed."
  write(*,*) "Data written in the file: blasius_solutions.dat"
  write(*,*) "Use the matlab file blasius_solutions.m to make plots."

 close(10)
 stop
!********************************************************************************

 contains

!********************************************************************************
!* This computes the solution set (u,v,f,f1,f2) at a given position (xp,yp),
!* using the one-step integration method described in the book,
!* "I do like CFD, VOL.1" (See Section 7.13.9).
!* 
!*  U_inf --->
!*
!*  ^ y
!*  |                                                .(xp,yp)
!*  |
!*  ----->x         ------------------------------------------------------------
!*                x=0                  flat plate
!*
!* Input ---------------------------------------------------
!*   (xp,yp) = position where the exact solution is sought.
!*        nu = kinematic viscosity
!*     U_inf = Freestream velocity (x-direction only)
!*
!* Output --------------------------------------------------
!*   (u,v) = velocity vector
!*    etap = self-similar variable, eta
!*       f = function f in the streamfunction: Eq.(7.13.46)
!*      f1 = derivative of f
!*      f2 = second derivative of f
!*
!* NB: The location (xp,yp) must be far enough from the leading edge of the
!*     plate, i.e., x >> 0. Otherwise the solution will not be self-similar,
!*     and the Blasius solution does not apply.
!*
!* NB: This subroutine requires the function, rhs(Fv), for computing the right
!*     hand side of the ordinary differential equation.
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* Katate Masatsuka, March 2010. http://www.cfdbooks.com
!********************************************************************************
 subroutine blasius(etap,u,v,f,f1,f2, xp,yp,nu,U_inf)
 implicit none
!Parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter :: zero=0.0_dp, two=2.0_dp
 real(dp), parameter :: third=1.0_dp/3.0_dp, sixth=1.0_dp/6.0_dp, half=0.5_dp
!Input and Output
 real(dp), intent( in)  :: xp, yp    ! position at which the solution is computed.
 real(dp), intent( in)  :: nu, U_inf ! viscosity and freestrem velocity
 real(dp), intent(out)  :: u, v, f, f1, f2, etap ! velocity, f-functions, etap
!Local variables
 real(dp)               :: Rex, f2_0
 real(dp)               :: eta, deta
 real(dp), dimension(3) :: Fv, K1, K2, K3
 logical                :: finish

  f2_0 = 0.3320573362151946_dp ! Pre-computed initial value: Eq.(7.13.62)
   Rex = U_inf*xp/nu           ! Reynolds number based on x: Eq.(7.13.47)
  etap = yp/xp*sqrt(Rex)       ! Variable eta              : Eq.(7.13.47)

!Increment for ODE integration: default = 1.0e-04
   deta = 1.0e-04_dp

!Integrate ODE up to eta = etap by the classical RK4.
! 1. Initial values.
      eta = zero
    Fv(1) = zero
    Fv(2) = zero
    Fv(3) = f2_0
   finish = .false.
! 2. Stepping to etap! (NB: no need to integrate if etap==0.)
   if (etap > zero) then
    do
     if (eta + deta > etap) then
      deta = etap - eta
      finish = .true.
     endif
     eta = eta + deta
      K1 = Fv + half*deta*rhs(Fv)
      K2 = Fv + half*deta*rhs(K1)
      K3 = Fv +      deta*rhs(K2)
      Fv = (K1 + two*K2 + K3 - Fv)*third + deta*rhs(K3)*sixth
     if (finish) exit
    end do
   endif

!Solution at eta=etap, i.e., (x,y)=(xp,yp).
   f = Fv(1)
  f1 = Fv(2)
  f2 = Fv(3)
   u = f1
   v = half/sqrt(Rex)*(etap*f1-f)

 return
 end subroutine blasius
!*******************************************************************************
!* Right hand side of the ODE: Equation (7.13.53) in "I do like CFD, VOL.1"
!*
!* This function is used in the subroutine, blasius().
!*
!* Katate Masatsuka, March 2010. http://www.cfdbooks.com
!*******************************************************************************
 function rhs(Fv) result(G)
 implicit none
 integer , parameter  :: sp = kind(1.0)
 integer , parameter  :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter  :: half=0.5_dp
 real(dp), intent(in) :: Fv(3)
 real(dp) :: G(3)
 real(dp) :: f, f1, f2
    f  = Fv(1) 
    f1 = Fv(2)
    f2 = Fv(3)
  G(1) = f1
  G(2) = f2
  G(3) = -half*f*f2
 end function rhs
!*******************************************************************************

 end program blasius_main
