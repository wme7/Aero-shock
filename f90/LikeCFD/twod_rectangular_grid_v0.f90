!********************************************************************************
!* This program generates 2D quad and triangular grids in a rectangular domain.
!*
!*
!*      Boundary information is set up for a shock-diffraction problem
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
!* 3-Step Generation:
!*
!* 1. Generate a temporary structured grid data for nodes: xs(i,j) and ys(i,j)
!* 2. Generate a 1D node array: x(1:nnodes), y(1:nnodes)
!* 3. Generate element connectivity data: tria(1:ntria,3), quad(1:nquad,4)
!*
!*
!* Inuput: 
!*        xmin, xmax = x-coordinates of the left and right ends
!*        ymin, ymax = y-coordinates of the bottom and top ends
!*                nx = number of nodes in x-direction
!*                ny = number of nodes in y-direction
!*
!*        (NOTE: All input parameters are defined inside the program.)
!*
!* Output:
!*        tria_grid_tecplot.dat = tecplot file of the triangular grid
!*        quad_grid_tecplot.dat = tecplot file of the quadrilateral grid
!*                    tria.grid = triangular-grid file for OSSAN solver
!*                    quad.grid = quadrilateral-grid file for OSSAN solver
!*                project.bcmap = file that contains boundary condition info
!*
!*        (NOTE: The grid file format is specific to the OSSAN solver.)
!*
!*        (NOTE: The grid file contains 5 boundary segments which is designed
!*               for the shock-diffraction problem. See write_grid_file() on the
!*               bottom of this file. Boudnary condition files project.bcmap
!*               specify the boundary condition to be applied to each segment.)
!*
!*
!*
!*        written by Dr. Katate Masatsuka (info[at]cfdbooks.com),
!*
!* the author of useful CFD books, "I do like CFD" (http://www.cfdbooks.com).
!*
!* This is Version 0 (December 2011).
!* This F90 code is written and made available for an educational purpose as well
!* as for generating grids for the OSSAN-Euler2D code.
!* This file may be updated in future.
!*
!* Katate Masatsuka, January 2012. http://www.cfdbooks.com
!********************************************************************************
 program twod_rectangular_grid

 implicit none

!Parameters
  integer, parameter :: sp = kind(1.0)
  integer, parameter :: p2 = selected_real_kind(2*precision(1.0_sp))
  real(p2) :: zero=0.0_p2, one=1.0_p2

!Input  - domain size and grid dimensions
 real(p2) :: xmin, xmax !Minimum x and Max x.
 real(p2) :: ymin, ymax !Minimum y and Max y
 integer  :: nx         !Number of nodes in x-direction
 integer  :: ny         !Number of nodes in y-direction

!Output - grid files
 character(80) :: datafile_tria_tec = "tria_grid_tecplot.dat"
 character(80) :: datafile_quad_tec = "quad_grid_tecplot.dat"
 character(80) :: datafile_tria = "tria.grid" !Triangular grid file for OSSAN solver
 character(80) :: datafile_quad = "quad.grid" !      Quad grid file for OSSAN solver
 character(80) :: datafile_bcmap = "project.bcmap" ! bc file for OSSAN solver

!Local variables
 real(p2), dimension(:,:), allocatable :: xs, ys !Structured grid data

 integer :: nnodes !Total number of nodes
 integer ::  ntria !Total number of triangles
 integer ::  nquad !Total number of quadrilaterals
 integer ::  inode !Local variables used in the 1D nodal array

 integer , dimension(:,:), allocatable :: tria !Triangle connectivity data
 integer , dimension(:,:), allocatable :: quad !Quad connectivity data
 real(p2), dimension(:)  , allocatable :: x, y !Nodal coordinates, 1D array

 real(p2) :: dx !Uniform grid spacing in x-direction = (xmax-xmin)/nx
 real(p2) :: dy !Uniform grid spacing in y-direction = (ymax-ymin)/ny
 integer  :: i, j, os

!--------------------------------------------------------------------------------
! 0. Define the grid size and allocate the structured grid data array.
!
! ymax --------------------------------------
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
!      .                                    .
! ymin --------------------------------------
!    xmin                                  xmax

!  Define the domain: here we define a unit square.

      xmin = zero
      xmax = one

      ymin = zero
      ymax = one

!  Define the grid size: the number of nodes in each direction.
!  NOTE: "ny" is better to be an odd number to place a node at the midpoint
!        on the left boundary, which will be a corner node in shock diffraction problem.
    nx = 401
    ny = 401

!  Allocate arrays.
   allocate(xs(nx,ny),ys(nx,ny))

!--------------------------------------------------------------------------------
! 1. Generate a structured 2D grid data, (i,j) data: go up in y-direction!
!
! j=5 o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |   On the left is an example:
!     |        |        |        |        |         nx = 5
! j=4 o--------o--------o--------o--------o         ny = 5
!     |        |        |        |        |
!     |        |        |        |        |
!     |        |        |        |        |
! j=3 o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!     |        |        |        |        |
! j=2 o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!     |        |        |        |        |
! j=1 o--------o--------o--------o--------o
!     i=1      i=2      i=3      i=4      i=5

     write(*,*) "Generating structured data..."

!  Compute the grid spacing in x-direction
    dx = (xmax-xmin)/real(nx-1)

!  Compute the grid spacing in y-direction
    dy = (ymax-ymin)/real(ny-1)

!  Generate nodes in the domain.

   do j = 1, ny  ! Go up in y-direction.
    do i = 1, nx ! Go to the right in x-direction.

     xs(i,j) = xmin + dx*real(i-1)
     ys(i,j) = ymin + dy*real(j-1)

    end do
   end do

!--------------------------------------------------------------------------------
! 2. Generate unstructured data: 1D array to store the node information.
!
!    - The so-called lexcographic ordering -
!
!   21       22       23       24       25
!     o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |   On the left is an example:
!   16|      17|      18|      19|      20|           nx = 5
!     o--------o--------o--------o--------o           ny = 5
!     |        |        |        |        |
!     |        |        |        |        |   nnodes = 5x5 = 25
!   11|      12|      13|      14|      15|
!     o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!    6|       7|       8|       9|      10|
!     o--------o--------o--------o--------o
!     |        |        |        |        |
!     |        |        |        |        |
!    1|       2|       3|       4|       5|
!     o--------o--------o--------o--------o
!
   write(*,*) "Generating 1D node array for unstructured grid data..."

!  Total number of nodes
   nnodes = nx*ny

!  Allocate the arrays
   allocate(x(nnodes),y(nnodes))

! Node data: the nodes are ordered in 1D array.

  do j = 1, ny   !Go up in y-direction.
   do i = 1, nx  !Go to the right in x-direction.

    inode = i + (j-1)*nx   !<- Node number in the lexcographic ordering
      x(inode) =   xs(i,j)
      y(inode) =   ys(i,j)

   end do
  end do

! Deallocate the structured data: xs and ys, which are not needed any more.
! - You guys helped me create the 1D array. Thanks!
  deallocate(xs, ys)

 write(*,*)
 write(*,*) " Nodes have been generated:"
 write(*,*) "       nx  =", nx
 write(*,*) "       ny  =", ny
 write(*,*) "    nx*ny  =", nx*ny
 write(*,*) "    nnodes =", nnodes
 write(*,*)
 write(*,*) " Now, generate elements..."
 write(*,*)

!--------------------------------------------------------------------------------
! 3. Generate unstructured element data:
!
!    We generate both quadrilateral and triangular grids.
!    Both grids are constructed in the unstructured (finite-element) data.
!

! Allocate arrays of triangular and quad connectivity data.

!   Number of quadrilaterals = (nx-1)(ny-1)

    allocate( quad((nx-1)*(ny-1)  ,4) )

!   Number of triangles = 2*(nx-1)*(ny-1)

    allocate(  tria(2*(nx-1)*(ny-1),3) )

! (1)Genearte a triangular grid

  write(*,*) "Generating triangular grid..."
  call generate_tria_grid
  write(*,*)
  write(*,*) " Number of triangles =", ntria
  write(*,*)
  write(*,*) "Writing a tecplot file for the triangular grid..."
  call write_tecplot_file(datafile_tria_tec)
  write(*,*) " --> File generated: ", datafile_tria_tec

  write(*,*) "Writing a grid file for the triangular grid..."
  call write_grid_file(datafile_tria)
  write(*,*) " --> File generated: ", datafile_tria

  write(*,*)

! (2)Generate a quadrilateral grid

  write(*,*) "Generating quad grid..."
  call generate_quad_grid
  write(*,*)
  write(*,*) " Number of quads =", nquad
  write(*,*)
  write(*,*) "Writing a tecplot file for the quadrilateral grid..."
  call write_tecplot_file(datafile_quad_tec)
  write(*,*) " --> File generated: ", datafile_quad_tec

  write(*,*) "Writing a grid file for the quadrilateral grid..."
  call write_grid_file(datafile_quad)
  write(*,*) " --> File generated: ", datafile_quad

! (3)Generate a mixed grid. (not implemented. I'll leave it to you! You can do it!)
!

! (4)Write a boundary condition file: to be read by OSSAN-Euler2D code
  write(*,*) "Generating bcmap file..."
  open(unit=1, file=datafile_bcmap, status="unknown", iostat=os)
  write(1,*) "Boundary Segment  Boundary Condition"
  write(1,*) "               1          freestream"
  write(1,*) "               2           slip_wall"
  write(1,*) "               3  outflow_supersonic"
  write(1,*) "               4  outflow_supersonic"
  write(1,*) "               5           slip_wall"
  close(1)

!--------------------------------------------------------------------------------

 write(*,*)
 write(*,*) "Successfully completed. Stop."

 stop

 contains



!********************************************************************************
! This subroutine generates triangles by constructing the connectivity data.
!********************************************************************************
 subroutine generate_tria_grid
 implicit none
!Local variables
 integer  :: i, j, inode, i1, i2, i3, i4

! No quads
 nquad = 0
  quad = 0

! Trianguler grid with right-up diagonals (i.e., / ).
!
!  inode+nx   inode+nx+1     i4      i3
!       o--------o           o--------o
!       |     .  |           |     .  |
!       |   .    |     or    |   .    |
!       | .      |           | .      |
!       o--------o           o--------o
!    inode    inode+1        i1      i2
!
! Triangle is defined by the counterclockwise ordering of nodes.

 ntria = 0

 do j = 1, ny-1
  do i = 1, nx-1

   inode = i + (j-1)*nx

!     Define the local numbers (see figure above)
      i1 = inode
      i2 = inode + 1
      i3 = inode + nx + 1
      i4 = inode + nx

           ntria = ntria + 1
    tria(ntria,1) = i1
    tria(ntria,2) = i2
    tria(ntria,3) = i3

           ntria = ntria + 1
    tria(ntria,1) = i1
    tria(ntria,2) = i3
    tria(ntria,3) = i4

  end do
 end do

 end subroutine generate_tria_grid
!********************************************************************************


!********************************************************************************
! This subroutine generates quads by constructing the connectivity data.
!********************************************************************************
 subroutine generate_quad_grid
 implicit none
!Local variables
 integer :: i, j, inode, i1, i2, i3, i4

! No triangles
  ntria = 0
   tria = 0

!
!  inode+nx   inode+nx+1     i4      i3
!       o--------o           o--------o
!       |        |           |        |
!       |        |     or    |        |
!       |        |           |        |
!       o--------o           o--------o
!     inode   inode+1        i1      i2
!
! Quad is defined by the counterclockwise ordering of nodes.

! Quadrilateral grid

 nquad = 0

 do j = 1, ny-1
  do i = 1, nx-1

   inode = i + (j-1)*nx
!     Define the local numbers (see figure above)
      i1 = inode
      i2 = inode + 1
      i3 = inode + nx + 1
      i4 = inode + nx

!  Order the quad counterclockwise:
          nquad = nquad + 1
   quad(nquad,1) = i1
   quad(nquad,2) = i2
   quad(nquad,3) = i3
   quad(nquad,4) = i4

  end do
 end do

 end subroutine generate_quad_grid
!********************************************************************************

!********************************************************************************
! This subroutine writes a tecplot file.
!********************************************************************************
 subroutine write_tecplot_file(datafile)
 implicit none
 character(80),            intent(in) :: datafile
 integer :: os
!--------------------------------------------------------------------------------
 open(unit=1, file=datafile, status="unknown", iostat=os)
 write(1,*) 'title = "grid"'
 write(1,*) 'variables = "x","y",'
 write(1,*) 'zone N=',nnodes,',E=',ntria+nquad,',ET=quadrilateral,F=FEPOINT'
!--------------------------------------------------------------------------------
 do i = 1, nnodes
  write(1,*) x(i),y(i)
 end do
!--------------------------------------------------------------------------------
!Triangles
 if (ntria > 0) then
  do i = 1, ntria
   write(1,*)  tria(i,1),  tria(i,2), tria (i,3),  tria(i,3) !The last one is a dummy.
  end do
 endif

!Quadrilaterals
 if (nquad > 0) then
  do i = 1, nquad
   write(1,*) quad(i,1), quad(i,2), quad(i,3), quad(i,4)
  end do
 endif
!--------------------------------------------------------------------------------
 close(1)
 end subroutine write_tecplot_file
!********************************************************************************

!********************************************************************************
! This subroutine writes a grid file to be read by a solver.
! NOTE: Unlike the tecplot file, this files contains boundary info.
!********************************************************************************
 subroutine write_grid_file(datafile)
 implicit none
 character(80),            intent(in) :: datafile
 integer :: os
!--------------------------------------------------------------------------------
 open(unit=1, file=datafile, status="unknown", iostat=os)

!--------------------------------------------------------------------------------
! Grid size: # of nodes, # of triangles, # of quadrilaterals
  write(1,*) nnodes, ntria, nquad

!--------------------------------------------------------------------------------
! Node data
  do i = 1, nnodes
   write(1,*) x(i), y(i)
  end do

!--------------------------------------------------------------------------------
! Triangle connectivity
  if (ntria > 0) then
   do i = 1, ntria
    write(1,*) tria(i,1), tria(i,2), tria(i,3) 
   end do
  endif

! Quad connectivity
  if (nquad > 0) then
   do i = 1, nquad
    write(1,*) quad(i,1), quad(i,2), quad(i,3), quad(i,4) 
   end do
  endif

! Boundary data:
! NOTE: These boundary data are specific to the shock diffraction problem.
!
!  Example: nx=ny=7
!
!   in = inflow
!    w = wall
!    e = outflow
!    o = interior nodes
!  inw = this node belongs to both inflow and wall boundaries.
!   we = this node belongs to both wall and outflow boundaries.
!
!   inw----w----w----w----w----w----we
!     |    |    |    |    |    |    |
!    in----o----o----o----o----o----e
!     |    |    |    |    |    |    |
!    in----o----o----o----o----o----e
!     |    |    |    |    |    |    |
!   inw----o----o----o----o----o----e
!     |    |    |    |    |    |    |
!     w----o----o----o----o----o----e
!     |    |    |    |    |    |    |
!     w----o----o----o----o----o----e
!     |    |    |    |    |    |    |
!    we----e----e----e----e----e----e
!

! Number of boundary segments
  write(1,*) 5

  write(1,*) (ny-1)/2+1  !Inflow
  write(1,*) (ny-1)/2+1  !Left Wall
  write(1,*)  nx         !Bottom Outflow
  write(1,*)  ny         !Right  Outflow
  write(1,*)  nx         !Top Wall

  write(1,*)

! Inflow boundary
  do j = ny, (ny-1)/2+1, -1
   i = 1
    write(1,*) i + (j-1)*nx
  end do

! Left wall boundary
  do j = (ny-1)/2+1, 1, -1
   i = 1
    write(1,*) i + (j-1)*nx
  end do

! Bottom outflow boundary
  do i = 1, nx
   j = 1
    write(1,*) i + (j-1)*nx
  end do

! Right outflow boundary
  do j = 1, ny
   i = nx
    write(1,*) i + (j-1)*nx
  end do

! Top wall boundary
  do i = nx, 1, -1
   j = ny
    write(1,*) i + (j-1)*nx
  end do

!--------------------------------------------------------------------------------
 close(1)
 end subroutine write_grid_file
!********************************************************************************


 end program twod_rectangular_grid
