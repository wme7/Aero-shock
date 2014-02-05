!********************************************************************************
!* Subroutines for generating grids and computing exact solutions for the famous
!* (often used for code verification of Euler codes) Ringleb's flow [F. Ringleb,
!*  Exakte Losungen der Differentialgleichungen einer adiabatischen Gasstromung,
!*  Z.A.M.M., 20(4), pp.185-198, 1940].
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* This is Version 1 (2009).
!* 
!* These F90 routines were written and made available for download
!* for an educational purpose. Also, this is intended to provide all CFD
!* students and researchers with a tool for their code verification.
!*
!* This file may be updated in future. (to improve the boundary point 
!*                                      distribution, extend the domain, etc.)
!*
!* Katate Masatsuka, March 2009. http://www.cfdbooks.com
!********************************************************************************

!********************************************************************************
!* --- Driver code for Ringleb's Flow ---
!* 
!* 1. Generate boundary points for a specified domain in Ringleb's flow.
!*    : The data is written in the file `ringleb_boundary_nodes.dat'.
!*      It can be used as an input for other grid generation programs.
!*
!* 2. Generate structured grids over the domain.
!*    : The grid data is written as a Tecplot file:
!*         ringleb_grid_soln_tec.dat     for quadrilateral grid,
!*         ringleb_grid_soln_tri_tec.dat for triangular grid.       
!*
!* 3. Compute the exact solution on the grid points.
!*    : To demonstrate the ability of the subroutine 'exact_solution_Ringleb()'
!*      to compute the exact solution at a given point location in a physical
!*      domain.
!*
!* Input ---------------------------------------------------
!*   psi_min = right boundary streamline.
!*   psi_max =  left boundary streamline.
!*   Vmin    =  upper boundary (constant speed inflow).
!*       nX  = # of cells in the horizontal direction.
!*       nY  = # of cells in the vertical direction.
!* 
!* Output --------------------------------------------------
!*  ringleb_boundary_nodes.dat: boundary point data file
!*  ringleb_grid_soln_tec.dat : Tecplot grid data file for quadrilateral grid
!*                              (with exact solutions)
!*  ringleb_grid_soln_tri_tec.dat: Tecplot grid data file for triangular grid
!*                                 (with exact solutions)
!*
!* 12-29-2010: Bugs fixed.
!*
!* Katate Masatsuka, March 2009. http://www.cfdbooks.com
!********************************************************************************
 program Ringleb
 implicit none
!Parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter :: one=1.0_dp
!Local
 integer :: i, j

!Variables and parameters (Grid generation part)
 real(dp), dimension(:,:), allocatable :: x,y, V,psi,rho,p,vx,vy,theta
 real(dp)                              :: psi_min, psi_max, V_min
 integer                               :: nx,ny

!Variables and parameters (Exact solution computation part)
 real(dp), dimension(:,:), allocatable :: Vp,psip,rhop,pp,vxp,vyp,thetap
 real(dp) :: xp,yp, error(7), error_inf(7)

!Input parameters
! This defines the domain.
  psi_min = 0.69_dp ! suggested value = 0.69_dp
  psi_max = 1.20_dp ! suggested value = 1.20_dp
    V_min = 0.50_dp ! suggested value = 0.50_dp
! This defines the grid size.
       nx = 30      ! the number of nodes in psi-direction
       ny = 60      ! the number of nodes in V-direction

!Allocate arrays
! Array for grid generation
  allocate(     x(nx+1, ny+1) )
  allocate(     y(nx+1, ny+1) )
  allocate(     V(nx+1, ny+1) )
  allocate(   psi(nx+1, ny+1) )
  allocate(   rho(nx+1, ny+1) )
  allocate(     p(nx+1, ny+1) )
  allocate(    vx(nx+1, ny+1) )
  allocate(    vy(nx+1, ny+1) )
  allocate( theta(nx+1, ny+1) )
! Arrays for exact solution computation
  allocate(    Vp(nx+1, ny+1) )
  allocate(  psip(nx+1, ny+1) )
  allocate(  rhop(nx+1, ny+1) )
  allocate(    pp(nx+1, ny+1) )
  allocate(   vxp(nx+1, ny+1) )
  allocate(   vyp(nx+1, ny+1) )
  allocate(thetap(nx+1, ny+1) )

! 1. Grid generation:
!    Generate boundary points, structured quadrilateral and triangular grids.
     call generate_grid_Ringleb(x,y,V,psi,rho,p,vx,vy,theta, &
                                psi_min,psi_max,V_min,nx,ny  )

! 2. Exact solution computation on the grid generated above: 
!    Compute the exact solution at each grid point
!    from the coordinates (xp, yp).
   write(*,*) "Computing exact solutions and errors ......."
   error     =        0.0D0
   error_inf = -1000000.0D0

!    Check if the exact solutions (those with p) match those used to generate
!    the grid.
  do j = 1, ny+1
   do i = 1, nx+1

           xp = x(i,j)
           yp = y(i,j)
!        Compute exact solutions at (xp,yp).
         call exact_solution_Ringleb(Vp(i,j),rhop(i,j),pp(i,j),psip(i,j), &
                             thetap(i,j),vxp(i,j),vyp(i,j),  xp,yp)

	 error(1) = error(1) + abs(     Vp(i,j) -     V(i,j) )
	 error(2) = error(2) + abs( thetap(i,j) - theta(i,j) )
	 error(3) = error(3) + abs(   psip(i,j) -   psi(i,j) )
	 error(4) = error(4) + abs(   rhop(i,j) -   rho(i,j) )
	 error(5) = error(5) + abs(     pp(i,j) -     p(i,j) )
	 error(6) = error(6) + abs(    vxp(i,j) -    vx(i,j) )
	 error(7) = error(7) + abs(    vyp(i,j) -    vy(i,j) )

     error_inf(1) = max( error_inf(1), abs(     Vp(i,j) -     V(i,j) ) )
     error_inf(2) = max( error_inf(2), abs( thetap(i,j) - theta(i,j) ) )
     error_inf(3) = max( error_inf(3), abs(   psip(i,j) -   psi(i,j) ) )
     error_inf(4) = max( error_inf(4), abs(   rhop(i,j) -   rho(i,j) ) )
     error_inf(5) = max( error_inf(5), abs(     pp(i,j) -     p(i,j) ) )
     error_inf(6) = max( error_inf(6), abs(    vxp(i,j) -    vx(i,j) ) )
     error_inf(7) = max( error_inf(7), abs(    vyp(i,j) -    vy(i,j) ) )

   end do
  end do

!   Errors between exact solutions computed for given grid points and those
!   from which the grid was generated. All errors should be zero.
    error = error / real((nx+1)*(ny+1))
    write(*,'(A14,7ES9.2)') "   L1 Errors: ", (error(i)    , i=1,7)
    write(*,'(A14,7ES9.2)') "L_inf Errors: ", (error_inf(i), i=1,7)

 stop
!--------------------------------------------------------------------------------
!* End of Main Program
!--------------------------------------------------------------------------------

 contains

!********************************************************************************
!* This generates boundary points and a computational (structured quadrilateral
!* and triangular) grid for a domain of Ringleb's flow, based on the algorithm
!* described in the book,
!*    "I do like CFD, VOL.1" by Katate Masatuka (http://www.cfdbooks.com).
!*
!* Input ------------------------------------------------------------
!*       nx: the number of grid points in x (in psi, to be precise)
!*       ny: the number of grid points in y (in  V , to be precise)
!*  psi_min:  left boundary coordinate (suggested value = 0.69)
!*  psi_max: right boundary coordinate (suggested value = 1.20)
!*    V_min:   top boundary coordinate (suggested value = 0.50)
!*
!* (The physical domain is defined by the last three parameters: psi_min, 
!*  psi_max, V_min. See Figure 7.11.12 in "I do like CFD, VOL.1".)
!*
!* Output -----------------------------------------------------------
!*          x(i,j), y(i,j)    : grid points (i=1,nx+1, j=1,ny+1)
!*          V(i,j), psi(i,j)  : flow speed and stream function
!*        rho(i,j),   p(i,j)  : density and pressure
!*         vx(i,j),  vy(i,j)  : velocity components
!*                 theta(i,j) : flow angle
!*       ringleb_boundary.dat : boundary point data (ordered counter-clockwise)
!*  ringleb_grid_soln_tec.dat : Tecplot file for a quadrilateral grid
!*                              (grid and exact solution)
!*  ringleb_grid_soln_tri_tec.dat : Tecplot file for triangular grid
!*                                  (grid and exact solution)
!*
!*
!* See "I do like CFD, VOL.1" for details.
!*
!* Katate Masatsuka, March 2009. http://www.cfdbooks.com
!********************************************************************************
 subroutine generate_grid_Ringleb( x, y, V, psi, rho, p, vx, vy, theta, &
                                        psi_min, psi_max, V_min, nx, ny)
 implicit none
!Parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter :: fifth=0.2_dp, half=0.5_dp
!Input and Output
 real(dp), intent(in) :: psi_min, psi_max, V_min ! Input
 integer , intent(in) :: nx, ny                  ! Input
 real(dp), dimension(nx+1, ny+1), intent(out) :: x,y,V,psi         ! Output
 real(dp), dimension(nx+1, ny+1), intent(out) :: rho,p,vx,vy,theta ! Output
!Local variables
 integer  :: i, j
 real(dp) :: V_max, b
!******************************************************************
! 1. Generate a structured grid: (psi,V) -> (x,y)
!    (See "I do like CFD, VOL.1", p.190.)
!******************************************************************
 write(*,*) "Generating boundary points ......."

  !*** Right boundary: i = nx+1
     V_max = one / psi_min
    do j = 1, ny+1
     psi(nx+1,j) = psi_min
     call V_from_theta(psi_min, V_min, V_max, j, ny, V(nx+1,j),theta(nx+1,j))
!       V(nx+1,j) = V_from_theta(psi_min, V_min, V_max, j, ny,theta(nx+1,j))
     call xy_from_psi_V(x(nx+1,j), y(nx+1,j),  psi_min, V(nx+1,j))
    end do
    write(*,*) "Right boundary done."

  !*** Left boundary: i = 1
     V_max = one / psi_max
    do j = 1, ny+1
        psi(1,j) = psi_max
      call V_from_theta(psi_max, V_min, V_max, j, ny,  V(1,j),theta(1,j))
!          V(1,j) = V_from_theta(psi_max, V_min, V_max, j, ny,theta(1,j))
     call xy_from_psi_V(x(1, j), y(1, j),  psi_max, V(1,j))
    end do
    write(*,*) "left boundary done."

  !*** Top boundary: j = ny+1
    do i = 1, nx+1
     psi(i,ny+1) = psi_max + (psi_min-psi_max)/real(nx)*real(i-1)
       V(i,ny+1) = V_min
     call xy_from_psi_V(x(i, ny+1), y(i, ny+1),  psi(i,ny+1), V_min)
    end do
    write(*,*) "Top boundary done."

  !*** Bottom boundary: j = 1
    do i = 2, nx
        psi(i,1) = psi_max + (psi_min-psi_max)/real(nx)*real(i-1)
           V_max = one / psi(i,1)
          V(i,1) = V_max
     call xy_from_psi_V(x(i,1), y(i,1),  psi(i,1), V_max)
    end do
    write(*,*) "Bottom boundary done."

! Write out a boundary node data: (x,y) in counter-clockwise order.
  write(*,*) "Writing boundary node data ......."
  open(unit=1, file ="ringleb_boundary_nodes.dat", status="unknown")

   bottom : do i = 1, nx   ! left to right
    write(1,*) x(i, 1), y(i, 1)
   end do bottom

   right : do j = 1, ny    ! bottom to top
    write(1,*) x(nx+1, j), y(nx+1, j)
   end do right

   top : do i = nx, 1, -1  ! right to left
    write(1,*) x(i,ny+1), y(i,ny+1)
   end do top

   left : do j = ny, 1, -1 ! top to bottom
    write(1,*) x(1,j), y(1,j)
   end do left

  close(1)

!Generate interior nodes to construct a structured mesh.
   write(*,*) "Generating interior points ......."
   do i = 2, nx
     do j = 2, ny
       psi(i,j) = psi_max + (psi_min-psi_max)/real(nx)*real(i-1)
          V_max = one / psi(i,j)
       call V_from_theta(psi(i,j), V_min, V_max, j, ny,  V(i,j),theta(i,j))
!         V(i,j) = V_from_theta(psi(i,j), V_min, V_max, j, ny,theta(i,j))
       call xy_from_psi_V(x(i,j), y(i,j),  psi(i,j), V(i,j))
     end do
   end do

   write(*,*) "Interior Points done."

!******************************************************************
! 2. Compute exact solutions from the grid generation data, i.e.
!    from the known values of psi and V at every grid point.
!    ( Iterations are not required here; the grid was generated
!      from the solution. See "I do like CFD, VOL.1", p.191. )
!******************************************************************
 write(*,*) "Computing exact solutions at grid points..."
  do j = 1, ny+1
   do i = 1, nx+1

   if ( abs(psi(i,j)*V(i,j)-one) < 1.0e-15_dp ) then
    theta(i,j) = half*3.141592653589793238_dp ! Do not rely on `asin' here.
   else
    theta(i,j) = asin(psi(i,j)*V(i,j)) ! Eq.(7.11.65)
   endif

       vx(i,j) = -V(i,j)*cos(theta(i,j))
       vy(i,j) = -V(i,j)*sin(theta(i,j))
             b = sqrt(one-fifth*V(i,j)*V(i,j))
      rho(i,j) = b**5 ! Eq.(7.11.68)
        p(i,j) = b**7
   end do
  end do

!******************************************************************
! 3. Write a Tecplot file of a structured grid which includes the
!    exact solution values.
!******************************************************************
 write(*,*) "Writing a Tecplot file (quadrilateral grid)..."
 open(unit=10, file ="ringleb_grid_soln_tec.dat", status="unknown")
 write(10,*) "TITLE = GRID" 
 write(10,*) "VARIABLES = x, y, vx, vy, rho, p, V, theta, Mach"
 write(10,*) "ZONE  N=",(nx+1)*(ny+1), ",E=" , nx*ny, &
             " ,ET=QUADRILATERAL, F=FEPOINT"

! nodes: coordinates and exact solutions.
  do j = 1, ny+1
   do i = 1, nx+1
    write(10,'(9ES20.10)') x(i,j),y(i,j),vx(i,j),vy(i,j),rho(i,j),p(i,j), &
        V(i,j),theta(i,j), V(i,j)/sqrt(abs(1.4_dp*p(i,j)/rho(i,j)))*sqrt(1.4_dp)
   end do
  end do

! elements
  do j = 0, ny-1
    do i = 0, nx-1
      write(10,*) (nx+1)*j+i+1        , (nx+1)*j+(i+1)+1, &
                  (nx+1)*(j+1)+(i+1)+1, (nx+1)*(j+1)+i+1
    end do
  end do

 close(10)

!******************************************************************
! 4. Write a Tecplot file of a triangular grid which includes the
!    exact solution values.
!******************************************************************
 write(*,*) "Writing a Tecplot file (triangular grid)..."
 open(unit=11, file ="ringleb_grid_soln_tri_tec.dat", status="unknown")
 write(11,*) "TITLE = GRID" 
 write(11,*) "VARIABLES = x, y, vx, vy, rho, p, V, theta, Mach"
 write(11,*) "ZONE  N=",(nx+1)*(ny+1), ",E=" , 2*nx*ny, &
             " ,ET=TRIANGLE, F=FEPOINT"

! nodes: coordinates and exact solutions.
  do j = 1, ny+1
   do i = 1, nx+1
    write(11,'(9ES20.10)') x(i,j),y(i,j),vx(i,j),vy(i,j),rho(i,j),p(i,j), &
        V(i,j),theta(i,j), V(i,j)/sqrt(abs(1.4_dp*p(i,j)/rho(i,j)))*sqrt(1.4_dp)
   end do
  end do

! elements
  do j = 0, ny-1
    do i = 0, nx-1
      write(11,*) (nx+1)*j+i+1, (nx+1)*j+(i+1)+1    , (nx+1)*(j+1)+(i+1)+1
      write(11,*) (nx+1)*j+i+1, (nx+1)*(j+1)+(i+1)+1, (nx+1)*(j+1)+i+1
    end do
  end do

 close(11)

 end subroutine generate_grid_Ringleb
!--------------------------------------------------------------------------------



!********************************************************************************
!* This computes exact solutions at a given location (xp,yp).
!*
!* This is a useful subroutine when you need exact solutions at a point
!* on a grid which is not generated from exact solutions ( other than those
!* generated by the subtourine generate_grid_soln_Ringleb() ).
!*
!* See "I do like CFD, VOL.1" for details.
!*
!* Katate Masatsuka, March 2009. http://www.cfdbooks.com
!********************************************************************************
 subroutine exact_solution_Ringleb(V,rho,p,psi,theta,vx,vy,  xp,yp)
 implicit none
!Parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter :: zero=0.0_dp, fifth=0.2_dp, half=0.5_dp
 real(dp), parameter ::  one=1.0_dp, three=3.0_dp, five=5.0_dp
!Input and Output

 real(dp), intent( in) :: xp,yp                    ! Input
 real(dp), intent(out) :: V,rho,p,psi,theta,vx,vy  ! Output
!Local variables
 real(dp) :: b, L, Vpsi

     V = zero
   rho = zero
     p = zero
   psi = zero
 theta = zero
    vx = zero
    vy = zero

!Start computing the exact solution by finding V from (xp,yp)
      V = V_from_xy_fpitr(xp, yp, 0.0_dp, 2.0_dp)

!Once V is obtained, everything else can be computed.
      b = sqrt(one - fifth*V*V)
    rho = b**five
      p = b**7
      L = one/b + one/(three*b**three) + one/(five*b**five) &
        - half*log( (one+b) / (one-b) )
    psi = sqrt( half/V**2-(xp-half*L)*rho )

   Vpsi = sqrt( half-V*V*(xp-half*L)*rho )

 if ( abs( Vpsi - one ) < 1.0e-15_dp  ) then
  theta = half*3.141592653589793238_dp ! Do not rely on `asin' here.
 else
  theta = asin( psi*V )
 endif

     vx = -V*cos(theta)
     vy = -V*sin(theta)

 end subroutine exact_solution_Ringleb


!********************************************************************************
!* This computes (x,y) from psi and V.
!*
!* See "I do like CFD, VOL.1" for details.
!*
!* Katate Masatsuka, March 2009. http://www.cfdbooks.com
!********************************************************************************
 subroutine xy_from_psi_V(x,y, psi,V)
 implicit none
!Parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter :: zero=0.0_dp, fifth=0.2_dp, half=0.5_dp
 real(dp), parameter :: one=1.0_dp, three=3.0_dp, five=5.0_dp
!Input and Output
 real(dp), intent( in) :: psi, V ! Input
 real(dp), intent(out) :: x, y   ! Output
!Local variables
 real(dp) :: b, rho, L

   b  = sqrt((one - fifth*V*V))
  rho = b**five
   L  =   one/b + one/(three*b**three) + one/(five*b**five) &
        - half*log((one+b)/(one-b))
   x  = (one/rho)*(half/V/V- psi*psi) + half*L ! Eq.(7.11.66)
   if ( abs(psi*V - one) < 1.0e-15_dp ) then
    y = zero   ! to avoid sqrt(negative) below.
   else
    y  = psi/(rho*V)*sqrt(one-psi*psi*V*V) ! Eq.(7.11.67)
   endif

 end subroutine xy_from_psi_V
!--------------------------------------------------------------------------------


!********************************************************************************
!* Adjust V through theta along the streamline psi.
!*
!* This is to generate a smooth node distribution on a boundary.
!* See "I do like CFD, VOL.1" by Katate Masatuka (http://www.cfdbooks.com)
!* for details. It computes and returns theta also which is based on a uniform
!* spacing in theta..
!*
!* Katate Masatsuka, March 2009. http://www.cfdbooks.com
!********************************************************************************
 subroutine V_from_theta(psi, V_min, V_max, j, ny, V, theta)
 implicit none
!Parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
!Input and output
 real(dp), intent(in)  :: psi, V_min, V_max   ! Input
 integer , intent(in)  :: j, ny               ! Input
 real(dp), intent(out) :: V, theta            ! Output
!Local variables
 real(dp) :: th_min, th_max
 
           th_min = asin(psi*V_min)
           th_max = asin(psi*V_max)
            theta = th_max + (th_min-th_max)/real(ny)*real(j-1)
                V = sin(theta) / psi

 end subroutine V_from_theta
!--------------------------------------------------------------------------------

!********************************************************************************
!* This solves F(V) = (x-L/2)^2 + y^2 - (1/2/rho/V^2)^2 = 0 for V by the 
!* fixed-point iteration.
!*
!* See "I do like CFD, VOL.1" for details.
!*
!* Katate Masatsuka, March 2009. http://www.cfdbooks.com
!********************************************************************************
 function V_from_xy_fpitr(x, y, Vmin, Vmax)
 implicit none
!Parameters
 integer , parameter :: sp = kind(1.0)
 integer , parameter :: dp = selected_real_kind(2*precision(1.0_sp))
 real(dp), parameter :: fifth=0.2_dp, half=0.5_dp
 real(dp), parameter :: one=1.0_dp, two=2.0_dp, three=3.0_dp, five=5.0_dp
!Input and output
 real(dp), intent(in) :: x, y, Vmin, Vmax ! Iput
 real(dp)             :: V_from_xy_fpitr  ! Output
!Local variables
 integer  :: i, imax
 real(dp) :: b, rho,L, V, VP, tol

 imax = 500
  tol = 1.0e-15_dp

!Initial guess
 Vp = half*( Vmin + Vmax )

!Start iterations
 fp_itertions : do i=1, imax

   b  = sqrt(one - fifth*Vp*Vp)
  rho = b**five
   L  =   one/b + one/(three*b**three) + one/(five*b**five) &
        - half*log((one+b)/(one-b))
   V  = sqrt(one/sqrt( abs((x-half*L)**two + y**two) )/two/rho) ! Eq.(7.11.76)

  if (abs(V-Vp) < tol) exit
  if (i == imax) then
    write(*,*) "did not converge... Sorry.", i
    stop
  endif

   Vp = V

 end do fp_itertions

  V_from_xy_fpitr = V

 end function V_from_xy_fpitr
!--------------------------------------------------------------------------------

 end program Ringleb







