!********************************************************************************
!*
!* --------- OSSAN (Oh, Such a Simple 'Ansutorakucha' Navier-Stokes) -----------
!*
!*            This program is the inviscid version: OSSAN-Euler2D
!*
!*             specially set up for a shock-diffraction problem
!*
!*                                 Wall
!*                         --------------------
!*     Post-shock (inflow) |                  |
!*                         |->Shock           |            o: Corner node
!*                         |  Mach=5.09       |
!*                  .......o                  |Outflow
!*                    Wall |                  |
!*                         |                  |
!*                         |                  |
!*                         --------------------
!*                               Outflow
!*
!* - Node-centered finite-volume method for unstructured grids (quad/tri/mixed)
!* - Roe flux with an entropy fix and Rotated-RHLL flux
!* - Gradient reconstruction by unweighted least-squares method
!* - Van Albada slope limiter to the primitive variable gradients
!* - 2-Stage Runge-Kutta global time-stepping towards the final time
!* - All quantities are nondimensionalized; velocity and pressure are
!*   nondimensionalized based on the free stream speed of sound
!*   (see Section 4.8.3 in I do like CFD, VOL.1).
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!*
!* This is Version 0 (January 2012)
!*
!* ------------------------------------------------------------------------------
!* Files: There are 3 files.
!*
!*  (1) data_package.f90    : Contains 3 data modules
!*
!*                           1. constants      - Numerical values are defined
!*                           2. grid_data_type - Grid data type is defined
!*                           3. my_main_data   - Main grid data is declared
!*
!*  (2) euler_solver.f90    : Integrate the Euler equations to the final time
!*  (3) main.f90            : Main driver code
!*
!* ------------------------------------------------------------------------------
!* Notes:
!*
!*  The purpose of this code is to give a beginner an opportunity to learn how to
!*  write an unstructured CFD code. Hence, the focus here is on the simplicity.
!*  The code is not optimized for efficiency.
!*
!*  This code is set up specifically for a shock diffraction problem.
!*  It can be modified easily to solve other problems:
!*
!*   1. Define your own free stream Mach number, M_inf, at the beginning of the main.
!*   2. Modify the subroutine, initial_solution_shock_diffraction(), which is in this
!*      file, to set up an appropriate initial condition for your problem.
!*   3. Delete the special treatment at the corner in euler_solver_main.f90
!*      (Special treatment is done in two places in that file.)
!*
!*  If the code is not simple enough to understand, please send questions to Hiro
!*  at sunmasen(at)hotmail.com. I'll greatly appreciate it and revise the code.
!*
!*  If the code helps you understand how to write your own code that is more
!*  efficient and has more features, it'll have served its purpose.
!*
!* ------------------------------------------------------------------------------
!* Examples of additional features you might want to add.
!*
!*  1. Local time-stepping      (faster convergence for steady problems)
!*  2. Implicit time-stepping   (remove CFL restriction)
!*  3. More boundary conditions (periodic, symmetry, suction, etc.)
!*  4. Other reconstruction     (Van Leer's kappa scheme)
!*  5. Other limiters           (Venkat/Barth limiter,etc.)
!*  6. Other flux functions     (HLL, LDFSS, AUSM, etc.)
!*  7. Local-preconditioning    (low-Mach number accuracy and stiffness removal)
!*  8. CFL ramping              (increase CFL gradually for a stable start-up)
!*  9. More output              (convergence history, etc.)
!* 10. Parallelization          (large-scale problems)
!* 11. Grid adaptation          (h-refinement, steady or unsteady)
!* 12. Adjoint capability       (aerodynamic design, adaptation, etc.)
!* 13. Moving grid              (sliding mesh, overset grid, etc.)
!* 14. Multigrid                (grid-independent convergence)
!* ------------------------------------------------------------------------------
!*
!* Katate Masatsuka, January 2012. http://www.cfdbooks.com
!********************************************************************************

!********************************************************************************
!* Main program: Node-centered finite-volume Euler code
!*
!* This code computes an unsteady solution of the Euler equations.
!* It is set up to solve a shock diffraction problem.
!* So, it is not really (but almost) a general purpose code.
!*
!* Input -------------------------------------------------------
!*
!*   project.grid  = grid file containing boundary information
!*   project.bcmap = file that specifies boundary conditions
!*
!*   (Note: See the subroutine "read_grid", which is in this file,
!*          for the format of these files.)
!*
!*  Parameters to be specified inside the main program:
!*
!*          M_inf = Upstream Mach number
!*          gamma = Ratio of specific heats (1.4 for air)
!*            CFL = CFL number (global time step)
!*        t_final = Final time to stop the calculation
!*  time_step_max = Max iterations (just a big enough number)
!*  inviscid_flux = Inviscid flux selection (1 = Roe, 2 = Rotated-RHLL)
!*   limiter_type = Limiter selection (1=Van Albada limiter, 0=No limiter)
!*
!*
!* Output ------------------------------------------------------
!*
!*  "project_tecplot.dat" = Tecplot file containing the grid and the solution.
!*
!*
!*  NOTE: The variables are nondimensionalized values (compressible scale),
!*           rho=rho/rho_inf, u=u/a_inf, v=v/a_inf, rho=p/(rho_inf*a_inf^2)
!*
!*  NOTE: Many variables are passed to subroutines via USE statement.
!*        Each module contains specific data, and they are accessed by USE.
!*
!*
!* Katate Masatsuka, January 2012. http://www.cfdbooks.com
!********************************************************************************
 program ossan_euler2d

 use constants   , only : p2
 use my_main_data, only : M_inf, gamma, CFL, time_step_max, inviscid_flux, &
                          limiter_type, t_final
 use euler_solver, only : euler_solver_main

 implicit none

!Inout data files
 character(80) :: datafile_grid_in  !Grid file
 character(80) :: datafile_bcmap_in !Boundary condition file
!Output data file
 character(80) :: datafile_tec      !Tecplot file for viewing the result.

! Set file names

 datafile_grid_in  = 'project.grid'
 datafile_bcmap_in = 'project.bcmap'
 datafile_tec      = 'project_tecplot.dat'

!--------------------------------------------------------------------------------
! Input Parameters

         M_inf  = 0.0_p2      ! Freestream Mach number to be set in the subroutine
                              !    -> "initial_solution_shock_diffraction"
                              !    (Specify M_inf here for other problems.)
          gamma = 1.4_p2      ! Ratio of specific heats
            CFL = 0.95_p2     ! CFL number
        t_final = 0.18_p2     ! Final time to stop the calculation.
  time_step_max = 5000        ! Max time steps (just a big enough number)
  inviscid_flux = "rhll"      ! = Rotated-RHLL      , "roe"  = Roe flux
   limiter_type = "vanalbada" ! = Van Albada limiter, "none" = No limiter

!--------------------------------------------------------------------------------
! Solve the Euler equations and write the output datafile.
!
!    NOTE: The subroutine, "euler_solver_main", is contained in euler_solver.f90.
!          Other subroutines are contained in this file.

! (1) Read grid files
      call read_grid(datafile_grid_in, datafile_bcmap_in)

! (2) Construct grid data
      call construct_grid_data

! (3) Check the grid data (It is always good to check them before use!)
      call check_grid_data

! (4) Set initial solution for a shock diffraction problem
!     (Re-write or replace it by your own subroutine for other problems.)
      call initial_solution_shock_diffraction

! (5) Compute the solution (March in time to the final time)
      call euler_solver_main

! (6) Write out the tecplot data file (Solutions at nodes)
      call write_tecplot_file(datafile_tec)

!--------------------------------------------------------------------------------

  write(*,*) "Successfully completed. Stop."

 stop

 contains

!********************************************************************************
!* Read the grid and the exact solution.
!* ------------------------------------------------------------------------------
!*  Input: datafile_grid_in  = filename of the grid file
!*         datafile_bcmap_in = filename of the bc file
!*
!* Output: nnodes, ncells, node(:), elm(:), bound(:) = data used in the solver
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
!* 1. "datafile_grid_in" is assumed to have been written in the following format:
!*
!*   -----------------------------------------------------------------------
!*    write(*,*) nnodes, ntria, nquad !Numbers of nodes, triangles and quads
!*
!*   do i = 1, nnodes
!*    write(*,*) x(i), y(i) !(x,y) coordinates of each node
!*   end do
!*
!*   do i = 1, ntria        !Nodes of triangles ordered counterclockwise
!*    write(*,*) node_1(i), node_2(i), node_3(i)
!*   end do
!*
!*   do i = 1, nquad        !Nodes of quadrilaterals ordered counterclockwise
!*    write(*,*) node_1(i), node_2(i), node_3(i), node_4(i)
!*   end do
!* 
!*    write(*,*) nbound     !Number of boundary segments
!*
!*   do i = 1, nbound
!*    write(*,*) nbnodes(i) !Number of nodes on each segment
!*   end do
!*
!*   do i = 1, nbound
!*    do j = 1, nbnodes(i)
!*     write(*,*) bnode(j)  !Node number of each node j in segment i
!*    end do
!*   end do
!*   -----------------------------------------------------------------------
!*
!*   NOTE: Add the first node to the end if the segment is closed
!*         (e.g., airfoil) The number of nodes will be the actual number + 1
!*         in that case.
!*
!*   NOTE: Boundary nodes must be ordered such that the domain is on the left.
!*
!********************************************************************************
!*
!* 2. "datafile_bcmap_in" is assumed have been written in the following format:
!*
!*   -----------------------------------------------------------------------
!*    write(*,*) "Boundary Segment              Boundary Condition"
!*   do i = 1, nbound
!*    write(*,*) i, bc_name
!*   end do
!*   -----------------------------------------------------------------------
!*
!*   NOTE: bc_name is the name of the boundary condition.
!*         Only four BCs are available in this version:
!*
!*         1. "freestream"
!*             Roe flux with freestream condition on the right state.
!*
!*         2. "slip_wall"
!*             Solid wall condition. Mass flux through the boundary is set zero.
!*
!*         3. "outflow_supersonic"
!*             Just compute the boundary flux by the physical Euler flux
!*             (equivalent to the interior-extrapolation condition.)
!*
!*         4. "outflow_back_pressure"
!*             Fix the back pressure. This should work for subsonic flows in a
!*             large enough domain.
!*
!********************************************************************************
!* Data to be read and stored:
!*
!* 1. Some numbers
!*    nnodes        = Number of nodes
!*    ntria         = Number of triangular elements
!*    nquad         = Number of quadrilateral elements
!*    nelms         = Total number of elements (=ntria+nquad)
!*
!* 2. Element data:
!*    elm(1:nelms)%nvtx   =  Number of vertices of each element
!*    elm(1:nelms)%vtx(:) = Pointer to vertices of each element
!*
!* 3. Node data: nodes are stored in a 1D array
!*    node(1:nnodes)%x     = x-coordinate of the nodes
!*    node(1:nnodes)%y     = y-coordinate of the nodes
!*
!* 4. Boundary Data:
!*    nbound                   = Number of boundary segments
!*    bound(1:nbound)%nbnodes  = Number of nodes in each segment
!*    bound(1:nbound)%bnode(:) = List of node numbers for each segment
!*    bound(1:nbound)%bc_type  = Boundary condition name for each segment
!*    bound(1:nbound)%bc_type  = Boundary condition name for each segment
!*
!********************************************************************************
 subroutine read_grid(datafile_grid_in, datafile_bcmap_in)

 use constants   , only : p2
 use my_main_data, only : nnodes, node, ntria, nquad, nelms, elm, nbound, bound

 implicit none
 character(80), intent(in) :: datafile_grid_in, datafile_bcmap_in

!Local variables
 integer  :: i, j, os, dummy_int

!--------------------------------------------------------------------------------
! 1. Read grid file>: datafile_grid_in

  write(*,*) "Reading the grid file....", datafile_grid_in

!  Open the input file.
   open(unit=1, file=datafile_grid_in, status="unknown", iostat=os)

! READ: Get the size of the grid.
  read(1,*) nnodes, ntria, nquad
  nelms = ntria + nquad

!  Allocate node and element arrays.
   allocate(node(nnodes))
   allocate(elm(  nelms))

! READ: Read the nodal coordinates
  do i = 1, nnodes
   read(1,*) node(i)%x, node(i)%y
  end do

! Read element-connectivity information

! Triangles: assumed that the vertices are ordered counterclockwise
!
!         v3
!         /\
!        /  \
!       /    \
!      /      \
!     /        \
!    /__________\
!   v1           v2

! READ: read connectivity info for triangles
  if ( ntria > 0 ) then
   do i = 1, ntria
    elm(i)%nvtx = 3
    allocate(elm(i)%vtx(3))
    read(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3)
   end do
  endif

! Quads: assumed that the vertices are ordered counterclockwise
!
!        v4________v3
!         /        |
!        /         |
!       /          |
!      /           |
!     /            |
!    /_____________|
!   v1             v2

! READ: read connectivity info for quadrilaterals
  if ( nquad > 0 ) then
   do i = 1, nquad
    elm(ntria+i)%nvtx = 4
    allocate( elm(ntria+i)%vtx(4))
    read(1,*) elm(ntria+i)%vtx(1), elm(ntria+i)%vtx(2), &
              elm(ntria+i)%vtx(3), elm(ntria+i)%vtx(4)
   end do
  endif

!  Write out the grid data.

   write(*,*) " Total numbers:"
   write(*,*) "      nodes = ", nnodes
   write(*,*) "  triangles = ", ntria
   write(*,*) "      quads = ", nquad
   write(*,*)

! Read the boundary grid data

! READ: Number of boundary condition types
  read(1,*) nbound
  allocate(bound(nbound))

! READ: Number of Boundary nodes (including the starting one at the end if
! it is closed such as an airfoil.)
  do i = 1, nbound
   read(1,*) bound(i)%nbnodes
   allocate(bound(i)%bnode(bound(i)%nbnodes))
  end do

! READ: Read boundary nodes
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes
   read(1,*) bound(i)%bnode(j)
   end do
  end do

!  Print the boundary grid data.
   write(*,*) " Boundary nodes:"
   write(*,*) "    segments = ", nbound
    do i = 1, nbound
     write(*,'(a9,i3,2(a11,i5))') " boundary", i, "  bnodes = ", bound(i)%nbnodes, &
                                                  "  bfaces = ", bound(i)%nbnodes-1
    end do
   write(*,*)

  close(1)

! End of Read grid file>: datafile_grid_in
!--------------------------------------------------------------------------------

!--------------------------------------------------------------------------------
! 2. Read the boundary condition data file

   write(*,*) "Reading the boundary condition file....", datafile_bcmap_in

! Open the input file.
  open(unit=2, file=datafile_bcmap_in, status="unknown", iostat=os)

    read(2,*) 

! READ: Read the boundary condition type
  do i = 1, nbound
    read(2,*) dummy_int, bound(i)%bc_type
   end do

!  Print the data
    write(*,*) " Boundary conditions:"
   do i = 1, nbound
    write(*,'(a9,i3,a12,a35)') " boundary", i, "  bc_type = ", trim(bound(i)%bc_type)
   end do

    i = dummy_int !Never mind. Just to avoid a compilation warning.

    write(*,*)

  close(2)

! End of Read the boundary condition data file
!--------------------------------------------------------------------------------

 end subroutine read_grid

!********************************************************************************
!* Construct the grid data:
!*
!* The following data, needed for NCFV method, will be constructed based on the
!* data read from the grid file.
!*
!* 1. Element data:
!*    elm(:)%nnghbrs  = Number of element neighbors of each element
!*    elm(:)%nghbr(:) = List of element neighbors of each element
!*    elm(:)%x        = x-coordinate of the centroid
!*    elm(:)%y        = y-coordinate of the centroid
!*    elm(:)%vol      = Volume of the element
!*
!*
!* 2. Node data:
!*    node(:)%nnghbrs = Number of node neighbors of each node
!*    node(:)%nghbr(:)= List of node neighbors of each node
!*    node(:)%nelms   = Number of adjacent elements of each node
!*    node(:)%elm     = List of adjacent elements of each node
!*    node(:)%vol     = Volume of the dual volume around each node
!*
!* 3. Edge data:
!*    edge(:)%n1, n2  = End nodes of each edge (edge points n1 -> n2)
!*    edge(:)%e1, e2  = Left and right elements of each edge
!*    edge(:)%dav     = Unit directed area vector of each edge
!*    edge(:)%da      = Magnitude of the directed area vector for each edge
!*    edge(:)%ev      = Unit edge vector of each edge (vector n1 -> n2)
!*    edge(:)%e       = Magnitude of the edge vector for each edge
!*
!*
!* 4. Boudnary data
!*    bound(:)%bnx    = Outward normal at boundary nodes (x-component of unit vector)
!*    bound(:)%bny    = Outward normal at boundary nodes (y-component of unit vector)
!*    bound(:)%bn     = Magnitude of (bnx,bny)
!*    bound(:)%bfnx   = Outward normal at boundary nodes (x-component of unit vector)
!*    bound(:)%bfny   = Outward normal at boundary nodes (y-component of unit vector)
!*    bound(:)%bfn    = Magnitude of (bfnx,bfny)
!*    bound(:)%belm   = Element to which the boundary face belongs
!*
!********************************************************************************
 subroutine construct_grid_data

 use my_main_data, only : nnodes, node, nelms, elm, nedges, edge, nbound, bound
 use constants   , only : p2, zero, half, third, fourth
 use euler_solver, only : lsq01_2x2_matrix_nc

 implicit none

!Local variables
 integer  ::  i, j, k, ii, in, im, jelm, v1, v2, v3, v4
 real(p2) :: x1, x2, x3, x4, y1, y2, y3, y4, xm, ym, xc, yc
 real(p2) :: xj, yj, xm1, ym1, xm2, ym2
 logical  :: found
 integer  :: vL, vR, n1, n2, e1, e2
 integer  :: vt1, vt2, ielm

 integer  :: ave_nghbr, min_nghbr, max_nghbr, imin, imax

! Some initialization
 v2 = 0
 vL = 0
 im = 0
 jelm = 0

  write(*,*) "Constructing grid data...."

! Initializations
  do i = 1, nnodes
   node(i)%nelms = 0
  end do
   nedges = 0

!--------------------------------------------------------------------------------
! Loop over elements and construct the fololowing data.
!
! 1. Surrounding elements: node(:)%nelms, node(:)%elm(:)
!
!    Example: Node i is surrounded by the eleemnts, 23, 101, 13, 41.
!             node(i)%nelms = 4
!             node(i)%elm(1) = 23
!             node(i)%elm(2) = 13
!             node(i)%elm(3) = 41
!             node(i)%elm(4) = 101
!
!        o-------o-------------o
!       /        |   .         |
!      /    23   |      41     |
!     o----------o-------------o
!      \        i \            |
!       \   101    \     13    |
!        \          \          | 
!         o----------o---------o
!
! 2. Element quantities  : elm(:)%x,elm(:)%y,elm(:)%vol
!
!  o-----------o            
!   \          |            o
!    \    (x,y)|           / \
!     \   .    |          /   \
!      \       |         /  .  \    (x,y): centroid coordinates
!       \      |        / (x,y) \     vol: volume of element
!        o-----o       o---------o

  elements : do i = 1, nelms

   v1 = elm(i)%vtx(1)
   v2 = elm(i)%vtx(2)
   v3 = elm(i)%vtx(3)

   x1 = node(v1)%x
   x2 = node(v2)%x
   x3 = node(v3)%x

   y1 = node(v1)%y
   y2 = node(v2)%y
   y3 = node(v3)%y

! Distribute the element index to nodes.

   node(v1)%nelms = node(v1)%nelms + 1
   call my_alloc_int_ptr(node(v1)%elm, node(v1)%nelms)
   node(v1)%elm(node(v1)%nelms) = i

   node(v2)%nelms = node(v2)%nelms + 1
   call my_alloc_int_ptr(node(v2)%elm, node(v2)%nelms)
   node(v2)%elm(node(v2)%nelms) = i

   node(v3)%nelms = node(v3)%nelms + 1
   call my_alloc_int_ptr(node(v3)%elm, node(v3)%nelms)
   node(v3)%elm(node(v3)%nelms) = i

! Compute the cell center and cell volume.
   tri_or_quad : if (elm(i)%nvtx==3) then

!   Triangle centroid and volume
    elm(i)%x   = third*(x1+x2+x3)
    elm(i)%y   = third*(y1+y2+y3)
    elm(i)%vol = tri_area(x1,x2,x3,y1,y2,y3)

   elseif (elm(i)%nvtx==4) then

!   OK, this is a quad. Get the 4th vertex.
    v4 = elm(i)%vtx(4)
    x4 = node(v4)%x
    y4 = node(v4)%y
!   Centroid: median dual
!   (Note: There is an alternative. See Appendix B in Nishikawa AIAA2010-5093.)
    xm1 = half*(x1+x2)
    ym1 = half*(y1+y2)
    xm2 = half*(x3+x4)
    ym2 = half*(y3+y4)
    elm(i)%x   = half*(xm1+xm2)
    elm(i)%y   = half*(ym1+ym2)
!   Volume is computed as a sum of two triangles: 1-2-3 and 1-3-4.
    elm(i)%vol = tri_area(x1,x2,x3,y1,y2,y3) + tri_area(x1,x3,x4,y1,y3,y4)

     xc = elm(i)%x
     yc = elm(i)%y
    if (tri_area(x1,x2,xc,y1,y2,yc)<zero) then
     write(*,*) " Centroid outside the quad element 12c: i=",i
     write(*,'(a10,2es10.2)') "  (x1,y1)=",x1,y1
     write(*,'(a10,2es10.2)') "  (x2,y2)=",x2,y2
     write(*,'(a10,2es10.2)') "  (x3,y3)=",x3,y3
     write(*,'(a10,2es10.2)') "  (x4,y4)=",x4,y4
     write(*,'(a10,2es10.2)') "  (xc,yc)=",xc,yc
     stop
    endif

    if (tri_area(x2,x3,xc,y2,y3,yc)<zero) then
     write(*,*) " Centroid outside the quad element 23c: i=",i
     write(*,'(a10,2es10.2)') "  (x1,y1)=",x1,y1
     write(*,'(a10,2es10.2)') "  (x2,y2)=",x2,y2
     write(*,'(a10,2es10.2)') "  (x3,y3)=",x3,y3
     write(*,'(a10,2es10.2)') "  (x4,y4)=",x4,y4
     write(*,'(a10,2es10.2)') "  (xc,yc)=",xc,yc
     stop
    endif

    if (tri_area(x3,x4,xc,y3,y4,yc)<zero) then
     write(*,*) " Centroid outside the quad element 34c: i=",i
     write(*,'(a10,2es10.2)') "  (x1,y1)=",x1,y1
     write(*,'(a10,2es10.2)') "  (x2,y2)=",x2,y2
     write(*,'(a10,2es10.2)') "  (x3,y3)=",x3,y3
     write(*,'(a10,2es10.2)') "  (x4,y4)=",x4,y4
     write(*,'(a10,2es10.2)') "  (xc,yc)=",xc,yc
     stop
    endif

    if (tri_area(x4,x1,xc,y4,y1,yc)<zero) then
     write(*,*) " Centroid outside the quad element 41c: i=",i
     write(*,'(a10,2es10.2)') "  (x1,y1)=",x1,y1
     write(*,'(a10,2es10.2)') "  (x2,y2)=",x2,y2
     write(*,'(a10,2es10.2)') "  (x3,y3)=",x3,y3
     write(*,'(a10,2es10.2)') "  (x4,y4)=",x4,y4
     write(*,'(a10,2es10.2)') "  (xc,yc)=",xc,yc
     stop
    endif

!  Distribution of element number to the 4th node of the quadrilateral
   node(v4)%nelms = node(v4)%nelms + 1
   call my_alloc_int_ptr(node(v4)%elm, node(v4)%nelms)
   node(v4)%elm(node(v4)%nelms) = i

   endif tri_or_quad

  end do elements

! Median dual volume

  do i = 1, nnodes
   node(i)%vol = zero
  end do

  elementsv : do i = 1, nelms

   v1 = elm(i)%vtx(1)
   v2 = elm(i)%vtx(2)
   v3 = elm(i)%vtx(3)

   tri_or_quadv : if (elm(i)%nvtx==3) then
!   Dual volume is exactly 1/3 of the volume of the triangle.
    node(v1)%vol = node(v1)%vol + third*elm(i)%vol
    node(v2)%vol = node(v2)%vol + third*elm(i)%vol
    node(v3)%vol = node(v3)%vol + third*elm(i)%vol

   elseif (elm(i)%nvtx==4) then
    v4 = elm(i)%vtx(4)

    x1 = node(v1)%x
    x2 = node(v2)%x
    x3 = node(v3)%x
    x4 = node(v4)%x
    xc = elm(i)%x

    y1 = node(v1)%y
    y2 = node(v2)%y
    y3 = node(v3)%y
    y4 = node(v4)%y
    yc = elm(i)%y

! - Vertex 1
     xj = node(v1)%x
     yj = node(v1)%y
    xm1 = half*(xj+x2)
    ym1 = half*(yj+y2)
    xm2 = half*(xj+x4)
    ym2 = half*(yj+y4)

!   Median volume is computed as a sum of two triangles.
    node(v1)%vol = node(v1)%vol + & 
                   tri_area(xj,xm1,xc,yj,ym1,yc) + tri_area(xj,xc,xm2,yj,yc,ym2)

! - Vertex 2
     xj = node(v2)%x
     yj = node(v2)%y
    xm1 = half*(xj+x3)
    ym1 = half*(yj+y3)
    xm2 = half*(xj+x1)
    ym2 = half*(yj+y1)

!   Median volume is computed as a sum of two triangles.
    node(v2)%vol = node(v2)%vol + &
                   tri_area(xj,xm1,xc,yj,ym1,yc) + tri_area(xj,xc,xm2,yj,yc,ym2)

! - Vertex 3
     xj = node(v3)%x
     yj = node(v3)%y
    xm1 = half*(xj+x4)
    ym1 = half*(yj+y4)
    xm2 = half*(xj+x2)
    ym2 = half*(yj+y2)

!   Median volume is computed as a sum of two triangles.
    node(v3)%vol = node(v3)%vol + &
                   tri_area(xj,xm1,xc,yj,ym1,yc) + tri_area(xj,xc,xm2,yj,yc,ym2)

! - Vertex 4
     xj = node(v4)%x
     yj = node(v4)%y
    xm1 = half*(xj+x1)
    ym1 = half*(yj+y1)
    xm2 = half*(xj+x3)
    ym2 = half*(yj+y3)

!   Median volume is computed as a sum of two triangles.
    node(v4)%vol = node(v4)%vol + &
                   tri_area(xj,xm1,xc,yj,ym1,yc) + tri_area(xj,xc,xm2,yj,yc,ym2)
 
   endif tri_or_quadv

  end do elementsv

!--------------------------------------------------------------------------------
! Loop over elements 2
!
!  Allocate elm(:)%nghbr(:) : elm(:)%nnghrs, elm(:)%nghr(:)
!  Construct element nghbr data: elm(:)%nghbr(:)
!  Order of neighbor elements [e1,e2,e3,..] are closely related to
!  the order of vertices [v1,v2,v3,..] (see below).
!
!          o------o
!          |      |                
!        v4|  e1  |v3                     v3
!    o-----o------o------o      o---------o------------o
!    |     |      |      |       .      .   .        .
!    | e2  |      |  e4  |        . e2 .     . e1  .
!    o-----o------o------o         .  .       .  .
!       v1 |     .v2              v1 o---------o v2   
!          | e3 .                     .   e3  .
!          |   .                        .    .
!          |  .                           . .
!          | .                             o
!          o
!

! Allocate the neighbor array

  do i = 1, nelms

!  3 neighbors for triangle
   if (elm(i)%nvtx==3) then

    elm(i)%nnghbrs = 3
    allocate(elm(i)%nghbr(3))

!  4 neighbors for quadrilateral
   elseif (elm(i)%nvtx==4) then

    elm(i)%nnghbrs = 4
    allocate(elm(i)%nghbr(4))

   endif

  end do

! Begin constructing the element-neighbor data

  elements2 : do i = 1, nelms

   elm_vertex : do k = 1, elm(i)%nvtx

!   Get the face of the element i:
!
!             vL      vR
!              o------o
!             /       |
!            /        |
!           o---------o
!
    if (k  < elm(i)%nvtx) vL = elm(i)%vtx(k+1)
    if (k == elm(i)%nvtx) vL = elm(i)%vtx(1)     
    vR = elm(i)%vtx(k)

!   Loop over the surrounding elements of the node vR,
!   and find the element neighbor from them.
    found = .false.
    elms_around_vR : do j = 1, node(vR)%nelms
    jelm = node(vR)%elm(j)

     edge_matching : do ii = 1, elm(jelm)%nvtx
                   v1 = elm(jelm)%vtx(ii)
      if (ii  > 1) v2 = elm(jelm)%vtx(ii-1)
      if (ii == 1) v2 = elm(jelm)%vtx(elm(jelm)%nvtx)

      if (v1==vR .and. v2==vL) then
       found = .true.
       im = ii+1
       if (im > elm(jelm)%nvtx) im = im - elm(jelm)%nvtx
       exit edge_matching
      endif
     end do edge_matching

     if (found) exit elms_around_vR

    end do elms_around_vR

     in = k + 2
     if (in > elm(i)%nvtx) in = in - elm(i)%nvtx

    if (found) then
     elm(   i)%nghbr(in) = jelm
     elm(jelm)%nghbr(im) = i
    else
     elm(   i)%nghbr(in) = 0
    endif

   end do elm_vertex

  end do elements2

!--------------------------------------------------------------------------------
! Edge-data for node-centered (edge-based) scheme.
!
! Loop over elements 3
! Construct edge data: edge(:)%n1, n2, e1, e2.
! Edge points from node n1 to node n2.
!
!      n2
!       o------------o
!     .  \         .
!    .    \   e2  .
!   .  e1  \    .
!  .        \ .         Directed area is positive: n1 -> n2
! o----------o         e1: left element
!             n1       e2: right element (e2 > e1 or e2 = 0)

! First count the number of edges.
!
! NOTE: Count edges only if the neighbor element number is
!       greater than the current element (i) to avoid double
!       count. Zero element number indicates that it is outside
!       the domain (boundary face).

  elements0 : do i = 1, nelms

   v1 = elm(i)%vtx(1)
   v2 = elm(i)%vtx(2)
   v3 = elm(i)%vtx(3)

   tri_quad0 : if (elm(i)%nvtx==3) then

    if ( elm(i)%nghbr(3) > i  .or. elm(i)%nghbr(3)==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(1) > i .or. elm(i)%nghbr(1)==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(2) > i .or. elm(i)%nghbr(2)==0 ) then
     nedges = nedges + 1
    endif

   elseif (elm(i)%nvtx==4) then

    v4 = elm(i)%vtx(4)

    if ( elm(i)%nghbr(3) > i .or. elm(i)%nghbr(3) ==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(4) > i .or. elm(i)%nghbr(4) ==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(1) > i .or. elm(i)%nghbr(1) ==0 ) then
     nedges = nedges + 1
    endif

    if ( elm(i)%nghbr(2) > i .or. elm(i)%nghbr(2) ==0 ) then
     nedges = nedges + 1
    endif

   endif tri_quad0

  end do elements0

! Allocate the edge array.
  allocate(edge(nedges))
  nedges = 0
  edge(:)%e1 = 0
  edge(:)%e2 = 0

! Construct the edge data:
!  two end nodes (n1, n2), and left and right elements (e1, e2)

  elements3 : do i = 1, nelms

   v1 = elm(i)%vtx(1)
   v2 = elm(i)%vtx(2)
   v3 = elm(i)%vtx(3)

! Triangular element
   tri_quad2 : if (elm(i)%nvtx==3) then

    if ( elm(i)%nghbr(3) > i  .or. elm(i)%nghbr(3)==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v1
     edge(nedges)%n2 = v2
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(3)
    endif

    if ( elm(i)%nghbr(1) > i .or. elm(i)%nghbr(1)==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v2
     edge(nedges)%n2 = v3
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(1)
    endif

    if ( elm(i)%nghbr(2) > i .or. elm(i)%nghbr(2)==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v3
     edge(nedges)%n2 = v1
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(2)
    endif

!  Quadrilateral element
   elseif (elm(i)%nvtx==4) then

    v4 = elm(i)%vtx(4)

    if ( elm(i)%nghbr(3) > i .or. elm(i)%nghbr(3) ==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v1
     edge(nedges)%n2 = v2
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(3)
    endif

    if ( elm(i)%nghbr(4) > i .or. elm(i)%nghbr(4) ==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v2
     edge(nedges)%n2 = v3
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(4)
    endif

    if ( elm(i)%nghbr(1) > i .or. elm(i)%nghbr(1) ==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v3
     edge(nedges)%n2 = v4
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(1)
    endif

    if ( elm(i)%nghbr(2) > i .or. elm(i)%nghbr(2) ==0 ) then
     nedges = nedges + 1
     edge(nedges)%n1 = v4
     edge(nedges)%n2 = v1
     edge(nedges)%e1 = i
     edge(nedges)%e2 = elm(i)%nghbr(2)
    endif

   endif tri_quad2

  end do elements3

! Loop over edges
! Construct edge vector and directed area vector.
!
! Edge vector is a simple vector pointing froom n1 to n2.
! For each edge, add the directed area vector (dav) from
! the left and right elements.
!
!              n2
!   o-----------o-----------o
!   |     dav   |  dav      |
!   |       ^   |   ^       |
!   |       |   |   |       |
!   |   c - - - m - - -c    |
!   |           |           |
!   |           |           |    m: edge midpoint
!   |           |           |    c: element centroid
!   o-----------o-----------o
!                n1
!
  edges : do i = 1, nedges

   n1 = edge(i)%n1
   n2 = edge(i)%n2
   e1 = edge(i)%e1
   e2 = edge(i)%e2
   xm = half*( node(n1)%x + node(n2)%x )
   ym = half*( node(n1)%y + node(n2)%y )

   edge(i)%dav = zero

! Contribution from the left element
  if (e1 > 0) then
   xc = elm(e1)%x
   yc = elm(e1)%y
   edge(i)%dav(1) = -(ym-yc)
   edge(i)%dav(2) =   xm-xc
  endif

! Contribution from the right element
  if (e2 > 0) then
   xc = elm(e2)%x
   yc = elm(e2)%y
   edge(i)%dav(1) = edge(i)%dav(1) -(yc-ym)
   edge(i)%dav(2) = edge(i)%dav(2) + xc-xm
  endif

  if (e1 < 0 .and. e2 < 0) then
   write(*,*) "!!!!! e1 and e2 are both negative... No way..."
  endif

! Magnitude and unit vector
   edge(i)%da  = sqrt( edge(i)%dav(1)**2 + edge(i)%dav(2)**2 )
   edge(i)%dav = edge(i)%dav / edge(i)%da

! Edge vector

  edge(i)%ev(1) = node(n2)%x - node(n1)%x
  edge(i)%ev(2) = node(n2)%y - node(n1)%y
  edge(i)%e     = sqrt( edge(i)%ev(1)**2 + edge(i)%ev(2)**2 )
  edge(i)%ev    = edge(i)%ev / edge(i)%e

  end do edges

!--------------------------------------------------------------------------------
! Construct node neighbor data:
!  pointers to the neighbor nodes(o)
!
!        o     o
!         \   / 
!          \ /
!     o-----*-----o
!          /|
!         / |
!        /  o        *: node in interest
!       o            o: neighbors (edge-connected nghbrs)
!

  do i = 1, nnodes
   node(i)%nnghbrs = 0
  end do

! Loop over edges and distribute the node numbers:

  edges4 : do i = 1, nedges

   n1 = edge(i)%n1
   n2 = edge(i)%n2

! (1) Add node1 to the neighbor list of n2
   node(n1)%nnghbrs = node(n1)%nnghbrs + 1
   call my_alloc_int_ptr(node(n1)%nghbr, node(n1)%nnghbrs)
   node(n1)%nghbr(node(n1)%nnghbrs) = n2

! (2) Add node2 to the neighbor list of n1
   node(n2)%nnghbrs = node(n2)%nnghbrs + 1
   call my_alloc_int_ptr(node(n2)%nghbr, node(n2)%nnghbrs)
   node(n2)%nghbr(node(n2)%nnghbrs) = n1

  end do edges4

!--------------------------------------------------------------------------------
! Boundary normal at nodes constructed by accumulating the contribution
! from each boundary face normal. This vector will be used to enforce
! the tangency condition, for example.
!
!
!        Interior domain      /
!                            o
!                  .        /
!                  .       /
! --o-------o-------------o
!           j   |  .  |   j+1
!               v  .  v
!
!        Left half added to the node j, and
!       right half added to the node j+1.
!

! Allocate and initialize the normal vector arrays
  do i = 1, nbound

   allocate(bound(i)%bnx(bound(i)%nbnodes))
   allocate(bound(i)%bny(bound(i)%nbnodes))
   allocate(bound(i)%bn( bound(i)%nbnodes))

   do j = 1, bound(i)%nbnodes
    bound(i)%bnx(j) = zero
    bound(i)%bny(j) = zero
    bound(i)%bn( j) = zero
   end do

  end do

! Compute the outward normals
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes-1

    x1 = node(bound(i)%bnode(j  ))%x
    y1 = node(bound(i)%bnode(j  ))%y

    x2 = node(bound(i)%bnode(j+1))%x
    y2 = node(bound(i)%bnode(j+1))%y

    bound(i)%bnx(j) = bound(i)%bnx(j) + half*( -(y1-y2) )
    bound(i)%bny(j) = bound(i)%bny(j) + half*(   x1-x2  )

    bound(i)%bnx(j+1) = bound(i)%bnx(j+1) + half*( -(y1-y2) )
    bound(i)%bny(j+1) = bound(i)%bny(j+1) + half*(   x1-x2  )

   end do
  end do

! Compute the magnitude and turn (bnx,bny) into a unit vector
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes

    bound(i)%bn(j)  = sqrt( bound(i)%bnx(j)**2 + bound(i)%bny(j)**2 )
    bound(i)%bnx(j) = bound(i)%bnx(j) / bound(i)%bn(j)
    bound(i)%bny(j) = bound(i)%bny(j) / bound(i)%bn(j)

   end do
  end do

!--------------------------------------------------------------------------------
! Boundary face data
!
!      |     Domain      |
!      |                 |
!      o--o--o--o--o--o--o  <- Boundary segment
!   j= 1  2  3  4  5  6  7
!
!   In the above case, nbnodes = 7, nbfaces = 6
!

  do i = 1, nbound
   bound(i)%nbfaces = bound(i)%nbnodes-1
   allocate(bound(i)%bfnx(    bound(i)%nbfaces   ))
   allocate(bound(i)%bfny(    bound(i)%nbfaces   ))
   allocate(bound(i)%bfn(     bound(i)%nbfaces   ))
   allocate(bound(i)%belm(    bound(i)%nbfaces   ))
  end do

! Boundary face vector: outward normal
  do i = 1, nbound
   do j = 1, bound(i)%nbfaces

    x1 = node(bound(i)%bnode(j  ))%x
    y1 = node(bound(i)%bnode(j  ))%y
    x2 = node(bound(i)%bnode(j+1))%x
    y2 = node(bound(i)%bnode(j+1))%y

    bound(i)%bfn(j)  =  sqrt( (x1-x2)**2 + (y1-y2)**2 )
    bound(i)%bfnx(j) = -(y1-y2) / bound(i)%bfn(j)
    bound(i)%bfny(j) =  (x1-x2) / bound(i)%bfn(j)

   end do
  end do

! Boundary normal vector at nodes: outward normal
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes-1

    x1 = node(bound(i)%bnode(j  ))%x
    y1 = node(bound(i)%bnode(j  ))%y
    x2 = node(bound(i)%bnode(j+1))%x
    y2 = node(bound(i)%bnode(j+1))%y

    bound(i)%bfn(j)  =  sqrt( (x1-x2)**2 + (y1-y2)**2 )
    bound(i)%bfnx(j) = -(y1-y2) / bound(i)%bfn(j)
    bound(i)%bfny(j) =  (x1-x2) / bound(i)%bfn(j)

   end do
  end do

! Find element adjacent to the face: belm
!
!  NOTE: This is useful to figure out what element
!        each boundary face belongs to. Boundary flux needs
!        special weighting depending on the element.
!
!      |_________|_________|________|
!      |         |         |        | 
!      |         |         |        | 
!      |_________|_________|________|
!      |         |         |        |     <- Grid (e.g., quads)
!      |         | elmb(j) |        |
!   ---o---------o---------o--------o---  <- Boundary segment
!                 j-th face
!
! elmb(j) is the element number of the element having the j-th boundary face.
!

  do i = 1, nbound
   do j = 1, bound(i)%nbfaces

!   bface is defined by the nodes v1 and v2.
    v1 = bound(i)%bnode(j  )
    v2 = bound(i)%bnode(j+1)

    found = .false.
!   Find the element having the bface from the elements
!   around the node v1.
    do k = 1, node(v1)%nelms
     ielm = node(v1)%elm(k)
     do ii = 1, elm(ielm)%nvtx
      in = ii
      im = ii+1
      if (im > elm(ielm)%nvtx) im = im - elm(ielm)%nvtx !return to 1
      vt1 = elm(ielm)%vtx(in)
      vt2 = elm(ielm)%vtx(im)
       if (vt1 == v1 .and. vt2 == v2) then
        found = .true.
        exit
       endif
     end do
     if (found) exit
    end do

    if (found) then
     bound(i)%belm(j) = ielm
    else
     write(*,*) " Boundary-adjacent element not found. Error..."
     stop
    endif

   end do
  end do

!--------------------------------------------------------------------------------
! Construct least-squares matrix for node-centered schemes.
!
!        o     o
!         \   / 
!          \ /
!     o-----*-----o
!          /|
!         / |
!        /  o        *: node in interest
!       o            o: neighbors (edge-connected nghbrs)
!

! Check the number of neighbor nodes (must have at least 2 neighbors)
  write(*,*) " --- Node neighbor data:"

  ave_nghbr = node(1)%nnghbrs
  min_nghbr = node(1)%nnghbrs
  max_nghbr = node(1)%nnghbrs
       imin = 1
       imax = 1
   if (node(1)%nnghbrs==2) then
    write(*,*) "--- 2 neighbors for the node = ", 1
   endif

  do i = 2, nnodes
   ave_nghbr = ave_nghbr + node(i)%nnghbrs
   if (node(i)%nnghbrs < min_nghbr) imin = i
   if (node(i)%nnghbrs > max_nghbr) imax = i
   min_nghbr = min(min_nghbr, node(i)%nnghbrs)
   max_nghbr = max(max_nghbr, node(i)%nnghbrs)
   if (node(i)%nnghbrs==2) then
    write(*,*) "--- 2 neighbors for the node = ", i
   endif
  end do

  write(*,*) "      ave_nghbr = ", ave_nghbr/nnodes
  write(*,*) "      min_nghbr = ", min_nghbr, " at node ", imin
  write(*,*) "      max_nghbr = ", max_nghbr, " at node ", imax
  write(*,*)

! Now, compute the inverse of the LSQ matrix at each node.

  do i = 1, nnodes
   call lsq01_2x2_matrix_nc(i)
  end do

 return

 end subroutine construct_grid_data

!********************************************************************************
!* Compute the area of the triangle defined by the nodes, 1, 2, 3.
!*
!*              3 (x3,y3)
!*              o 
!*             / \ 
!*            /   \
!* (x1,y1) 1 o-----o 2 (x2,y2)
!*
!* Nodes must be ordered counterclockwise (otherwise it gives negative area)
!*
!********************************************************************************
 function tri_area(x1,x2,x3,y1,y2,y3) result(area)
 use constants, only : p2, half
 implicit none
 real(p2), intent(in) :: x1,x2,x3,y1,y2,y3
 real(p2) :: area

  area = half*( x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2) )

 end function tri_area


!********************************************************************************
!* Check the grid data.
!*
!* 1. Directed area must sum up to zero around every node.
!* 2. Directed area must sum up to zero over the entire grid.
!* 3. Global sum of the boundary normal vectors must vanish.
!* 4. Global sum of the boundary face normal vectors must vanish.
!* 5. Check element volumes which must be positive.
!* 6. Check dual volumes which must be positive.
!* 7. Global sum of the dual volumes must be equal to the sum of element volumes.
!* 8. Linear LSQ gradients must be exact for linear functions.
!*
!* Add more tests you can think of.
!*
!********************************************************************************
 subroutine check_grid_data

 use my_main_data  , only : nnodes, node,  nelms,   elm, nedges, edge, &
                            nbound, bound
 use constants     , only : p2, zero, half, third, fourth
 use euler_solver  , only : lsq01_2x2_gradients_nc

 implicit none
!Local variables
 integer  :: i, j, k, n1, n2, ierr
 real(p2), dimension(nnodes,2) :: sum_dav_i
 real(p2), dimension(2) :: sum_dav, sum_bn
 real(p2), dimension(2) :: sum_bfn
 real(p2)               :: sum_volc, sum_vol
 real(p2) :: c0_linear, cx_linear, cy_linear

  write(*,*) "Checking grid data...."

!--------------------------------------------------------------------------------
! Directed area sum check
!--------------------------------------------------------------------------------

! Compute the sum of the directed area for each node.

   sum_dav_i = zero
  do i = 1, nedges
   n1 = edge(i)%n1
   n2 = edge(i)%n2
   sum_dav_i(n1,:) = sum_dav_i(n1,:) + edge(i)%dav(:)*edge(i)%da
   sum_dav_i(n2,:) = sum_dav_i(n2,:) - edge(i)%dav(:)*edge(i)%da
  end do

! Compute also the sum of the boundary normal vector (at nodes).

  sum_bn = 0
  do i = 1, nbound
   do j = 1, bound(i)%nbnodes
     k = bound(i)%bnode(j)
    sum_dav_i(k,1) = sum_dav_i(k,1) + bound(i)%bnx(j)*bound(i)%bn(j)
    sum_dav_i(k,2) = sum_dav_i(k,2) + bound(i)%bny(j)*bound(i)%bn(j)
    sum_bn(1) = sum_bn(1) + bound(i)%bnx(j)*bound(i)%bn(j)
    sum_bn(2) = sum_bn(2) + bound(i)%bny(j)*bound(i)%bn(j)
   end do
  end do

! Global sum of boundary normal vectors must vanish.

  if (sum_bn(1) > 1.0e-12_p2 .and. sum_bn(2) > 1.0e-12_p2) then
   write(*,*) "--- Global sum of the boundary normal vector:"
   write(*,'(a19,es10.3)') "    sum of bn_x = ", sum_bn(1)
   write(*,'(a19,es10.3)') "    sum of bn_y = ", sum_bn(2)
   write(*,*) "Error: boundary normal vectors do not sum to zero..."
   stop
  endif

! Sum of the directed area vectors must vanish at every node.

  do i = 1, nnodes
   if (abs(sum_dav_i(i,1))>1.0e-12_p2 .or. abs(sum_dav_i(i,2))>1.0e-12_p2) then
   write(*,'(a11,i5,a7,2es10.3,a9,2es10.3)') &
    " --- node=", i, " (x,y)=", node(i)%x, node(i)%y, " sum_dav=",sum_dav_i(i,:)
   endif
  end do

   write(*,*) "--- Max sum of directed area vector around a node:"
   write(*,*) "  max(sum_dav_i_x) = ", maxval(sum_dav_i(:,1))
   write(*,*) "  max(sum_dav_i_y) = ", maxval(sum_dav_i(:,2))

  if (maxval(abs(sum_dav_i(:,1)))>1.0e-12_p2 .or. &
      maxval(abs(sum_dav_i(:,2)))>1.0e-12_p2) then
   write(*,*) "--- Max sum of directed area vector around a node:"
   write(*,*) "  max(sum_dav_i_x) = ", maxval(sum_dav_i(:,1))
   write(*,*) "  max(sum_dav_i_y) = ", maxval(sum_dav_i(:,2))
   write(*,*) "Error: directed area vectors do not sum to zero..."
   stop
  endif

! Of course, the global sum of the directed area vector sum must vanish.
   sum_dav = zero
  do i = 1, nnodes
   sum_dav = sum_dav + sum_dav_i(i,:)
  end do

   write(*,*) "--- Global sum of the directed area vector:"
   write(*,'(a19,es10.3)') "    sum of dav_x = ", sum_dav(1)
   write(*,'(a19,es10.3)') "    sum of dav_y = ", sum_dav(2)

  if (sum_dav(1) > 1.0e-12_p2 .and. sum_dav(2) > 1.0e-12_p2) then
   write(*,*) "Error: directed area vectors do not sum globally to zero..."
   write(*,*) "--- Global sum of the directed area vector:"
   write(*,'(a19,es10.3)') "    sum of dav_x = ", sum_dav(1)
   write(*,'(a19,es10.3)') "    sum of dav_y = ", sum_dav(2)
   stop
  endif

!--------------------------------------------------------------------------------
! Global sum check for boundary face vector
!--------------------------------------------------------------------------------
  sum_bfn = 0
  do i = 1, nbound
   do j = 1, bound(i)%nbfaces
     sum_bfn(1) =  sum_bfn(1) + bound(i)%bfnx(j)*bound(i)%bfn(j)
     sum_bfn(2) =  sum_bfn(2) + bound(i)%bfny(j)*bound(i)%bfn(j)
   end do
  end do

   write(*,*) "--- Global sum of the boundary face vector:"
   write(*,'(a19,es10.3)') "    sum of bfn_x = ", sum_bfn(1)
   write(*,'(a19,es10.3)') "    sum of bfn_y = ", sum_bfn(2)

  if (sum_bn(1) > 1.0e-12_p2 .and. sum_bn(2) > 1.0e-12_p2) then
   write(*,*) "Error: boundary face normals do not sum globally to zero..."
   write(*,*) "--- Global sum of the boundary face normal vector:"
   write(*,'(a19,es10.3)') "    sum of bfn_x = ", sum_bn(1)
   write(*,'(a19,es10.3)') "    sum of bfn_y = ", sum_bn(2)
   stop
  endif

!--------------------------------------------------------------------------------
! Volume check
!--------------------------------------------------------------------------------
! (1)Check the element volume: make sure there are no zero or negative volumes
       ierr = 0
   sum_volc = zero
  do i = 1, nelms

   sum_volc = sum_volc + elm(i)%vol

   if (elm(i)%vol < zero) then
     write(*,*) "Negative volc=",elm(i)%vol, " elm=",i, " stop..."
     ierr = ierr + 1
   endif

   if (abs(elm(i)%vol) < 1.0e-14_p2) then
     write(*,*) "Vanishing volc=",elm(i)%vol, " elm=",i, " stop..."
     ierr = ierr + 1
   endif

  end do

!--------------------------------------------------------------------------------
! (2)Check the dual volume (volume around a node)
      ierr = 0
   sum_vol = zero
  do i = 1, nnodes

   sum_vol = sum_vol + node(i)%vol

   if (node(i)%vol < zero) then
     write(*,*) "Negative vol=",node(i)%vol, " node=",i, " stop..."
     ierr = ierr + 1
   endif

   if (abs(node(i)%vol) < 1.0e-14_p2) then
     write(*,*) "Vanishing vol=",node(i)%vol, " node=",i, " stop..."
     ierr = ierr + 1
   endif

  end do

  if (ierr > 0) stop

  if (abs(sum_vol-sum_volc) > 1.0e-11_p2) then
   write(*,*) "--- Global sum of volume: must be the same"
   write(*,'(a19,es10.3)') "    sum of volc = ", sum_volc
   write(*,'(a19,es10.3)') "    sum of vol  = ", sum_vol
   write(*,'(a22,es10.3)') " sum_vol-sum_volc  = ", sum_vol-sum_volc
   write(*,*) "Error: sum of dual volumes and cell volumes do not match..."
   stop
  endif

!--------------------------------------------------------------------------------
! Check the least-squares matrix for node-centered scheme.
! Check if the LSQ gradients are exact for a linear function.
!--------------------------------------------------------------------------------
  write(*,*) "Checking the least-squares matrix(nc,linear)..."

! Store a linear function
  c0_linear = 5.0_p2
  cx_linear = 3.7_p2
  cy_linear = 2.3_p2

  do i = 1, nnodes
   node(i)%w = cx_linear*node(i)%x + cy_linear*node(i)%y + c0_linear
  end do

!-----------------------------------------------------------------
! Check the 2x2 gradients for node-centered schemes.
!-----------------------------------------------------------------
! Compute gradients

  do i = 1, nnodes
   call lsq01_2x2_gradients_nc(i)
  end do

! Check the gradients
  do i = 1, nnodes

!  Check dw/dx
   if (abs( maxval(node(i)%gradw(:,1)) - cx_linear) > 1.0e-11_p2) then
    write(*,*) " i = ", i, "  nnghbrs=",node(i)%nnghbrs
    write(*,*) " Max error = ", maxval(node(i)%gradw(:,1)) - cx_linear
    write(*,'(a9,9es10.2)') " gradw_x=", node(i)%gradw(:,1)
    stop
   endif

!  Check dw/dy
   if (abs( maxval(node(i)%gradw(:,2)) - cy_linear) > 1.0e-11_p2) then
    write(*,*) " i = ", i, "  nnghbrs=",node(i)%nnghbrs
    write(*,*) " Max error = ", maxval(node(i)%gradw(:,2)) - cy_linear
    write(*,'(a9,9es10.2)') " gradw_y=", node(i)%gradw(:,2)
    stop
   endif

  end do

! Re-initialuze the solution.
  do i = 1, nnodes
   node(i)%w = zero
  end do

  write(*,*) "LSQ gradients passed the linear function test! Good."

!--------------------------------------------------------------------------------

  write(*,*)
  write(*,*) "Grid data look good!"

 end subroutine check_grid_data



!********************************************************************************
!* Initial solution for the shock diffraction problem:
!*
!* NOTE: So, this is NOT a general purpose subroutine.
!*       For other problems, specify M_inf in the main program, and
!*       modify this subroutine to set up an appropriate initial solution.
!
!  Shock Diffraction Problem:
!
!                             Wall
!                     --------------------
! Post-shock (inflow) |                  |
! (rho,u,v,p)_inf     |->Shock (M_shock) |            o: Corner node
!    M_inf            |                  |
!              .......o  Pre-shock       |Outflow
!                Wall |  (rho0,u0,v0,p0) |
!                     |                  |
!                     |                  |
!                     --------------------
!                           Outflow
!
!********************************************************************************
 subroutine initial_solution_shock_diffraction

 use constants   , only : p2, zero, one, two, half
 use my_main_data, only : nnodes, node, gamma, M_inf, rho_inf, u_inf, v_inf, p_inf
 use euler_solver, only : w2u

 implicit none

!Local variables
 integer  :: i
 real(p2) :: M_shock, u_shock, rho0, u0, v0, p0

  do i = 1, nnodes

! Pre-shock state: uniform state; no disturbance has reahced yet.

       rho0 = one
         u0 = zero
         v0 = zero
         p0 = one/gamma

! Incoming shock speed

    M_shock = 5.09_p2
    u_shock = M_shock * sqrt(gamma*p0/rho0)

! Post-shock state: These values will be used in the inflow boundary condition.
   rho_inf = rho0 * (gamma + one)*M_shock**2/( (gamma - one)*M_shock**2 + two )
     p_inf =   p0 * (   two*gamma*M_shock**2 - (gamma - one) )/(gamma + one)
     u_inf = (one - rho0/rho_inf)*u_shock
     M_inf = u_inf / sqrt(gamma*p_inf/rho_inf)
     v_inf = zero

! Set the initial solution: set the pre-shock state inside the domain.

   node(i)%w = (/ rho0, u0, v0, p0 /)
   node(i)%u = w2u(node(i)%w)

  end do

 end subroutine initial_solution_shock_diffraction


!********************************************************************************
!* Write a tecplot file: grid and solution
!********************************************************************************
 subroutine write_tecplot_file(datafile_tec)

 use constants   , only : one
 use my_main_data, only : nnodes, node, elm, nelms, gamma

 implicit none
 character(80), intent(in) :: datafile_tec
 integer :: i, k, os
!--------------------------------------------------------------------------------
 open(unit=1, file=datafile_tec, status="unknown", iostat=os)

 write(1,*) 'title = "grid"'

 write(1,'(a80)') 'variables = "x","y","rho","u","v","p","Mach"'

 write(1,*) 'zone n=',nnodes,' e =', nelms,' et=quadrilateral, f=fepoint'

!--------------------------------------------------------------------------------
! Nodal quantities: x, y, rho, u, v, p, Mach number

   do i = 1, nnodes
    write(1,*) node(i)%x, node(i)%y, (node(i)%w(k),k=1,4), &
               sqrt( (node(i)%w(2)**2+node(i)%w(3)**2)/(gamma*node(i)%w(4)/node(i)%w(1)) )
   end do

!--------------------------------------------------------------------------------
! Both quad and tria elements in quad format:

 do i = 1, nelms

  if (elm(i)%nvtx == 3) then

   write(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3), elm(i)%vtx(3)

  elseif (elm(i)%nvtx == 4) then

   write(1,*) elm(i)%vtx(1), elm(i)%vtx(2), elm(i)%vtx(3), elm(i)%vtx(4)

  else

   !Impossible
   write(*,*) " Error in elm%vtx data... Stop..: elm(i)%nvtx=",elm(i)%nvtx
   stop

  endif

 end do

!--------------------------------------------------------------------------------
 close(1)
 end subroutine write_tecplot_file
!********************************************************************************

!********************************************************************************
!* This subroutine is useful to expand or shrink integer arrays.
!*
!*  Array, x, will be allocated if the requested dimension is 1 (i.e., n=1)
!*  Array, x, will be expanded to the requested dimension, n, if (n > dim(x)).
!*  Array, x, will be shrinked to the requested dimension, n, if (n < dim(x)).
!*
!********************************************************************************
  subroutine my_alloc_int_ptr(x,n)
  implicit none
  integer, intent(in) :: n
  integer, dimension(:), pointer :: x
  integer, dimension(:), pointer :: temp
  integer :: i

  if (n <= 0) then
   write(*,*) "my_alloc_int_ptr received non-positive dimension. Stop."
   stop
  endif

! If not allocated, allocate and return
  if (.not.(associated(x))) then
   allocate(x(n))
   return
  endif

! If reallocation, create a pointer with a target of new dimension.
  allocate(temp(n))
   temp = 0

! (1) Expand the array dimension
  if ( n > size(x) ) then

   do i = 1, size(x)
    temp(i) = x(i)
   end do

! (2) Shrink the array dimension: the extra data, x(n+1:size(x)), discarded.
  else

   do i = 1, n
    temp(i) = x(i)
   end do

  endif

! Destroy the target of x
  deallocate(x)

! Re-assign the pointer
   x => temp

  return

  end subroutine my_alloc_int_ptr
!********************************************************************************

 end program ossan_euler2d




