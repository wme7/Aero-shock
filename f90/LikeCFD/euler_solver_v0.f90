!********************************************************************************
!* --------- OSSAN (Oh, Such a Simple 'Ansutorakucha' Navier-Stokes) -----------
!*
!*         This module belongs to the inviscid version: OSSAN-Euler2D
!*
!*
!* This module contains all subroutines needed to advance the solution in time.
!*
!* Contained subroutines/functions:
!* -----------------------------------------------------------------------------
!* euler_solver_main          : Main time-stepping routine
!* eliminate_normal_mass_flux : Eliminate any normal mass flux for slip condition
!* compute_residual_ncfv      : Compute the residual (RHS, spatial part)
!* update_solution            : Update the solution at all nodes
!* va_slope_limiter           : Van Albada slope limiter
!* roe                        : Roe flux with an entropy fix
!* rotated_rhll               : Rotated-RHLL flux (JCP2008) for strong shocks
!* compute_time_step          : Compute global/local time step
!* residual_norm              : Compute the residual norms (L1, L2, Linf)
!* lsq01_2x2_gradients_nc     : Gradient computation at a node
!* lsq01_2x2_matrix_nc        : Computes the inverse of LSQ matrix
!* w2u                        : Primitive to conservative variables
!* u2w                        : Conservative to primitive variables
!* -----------------------------------------------------------------------------
!*
!* Written by Katate Masatsuka (http://www.cfdbooks.com)
!********************************************************************************
 module euler_solver

 private

 public :: euler_solver_main
 public :: lsq01_2x2_gradients_nc
 public :: lsq01_2x2_matrix_nc
 public :: w2u

 contains

!********************************************************************************
!********************************************************************************
!********************************************************************************
!* Euler solver: Node-Centered Finite-Volume Method (Edge-Based)
!*
!* - Node-centered finite-volume method for unstructured grids(quad/tri/mixed)
!* - Roe flux with an entropy fix and Rotated-RHLL flux
!* - Reconstruction by unweighted least-squares method (2x2 system for gradients)
!* - Van Albada slope limiter to the primitive variable gradients
!* - 2-Stage Runge-Kutta time-stepping
!*
!********************************************************************************
 subroutine euler_solver_main

 use constants   , only : p2, zero, one, half
 use my_main_data, only : CFL, time_step_max, M_inf, nnodes, node, &
                          inviscid_flux, t_final, limiter_type
 implicit none

!Local variables
 real(p2), dimension(4,3)              :: res_norm !Residual norms(L1,L2,Linf)
 real(p2), dimension(:,:), allocatable :: u0       !Saved solution
 real(p2) :: dt, time    !Time step and actual time
 integer  :: i_time_step !Number of time steps
 integer  :: i

! Allocate the temporary solution array needed for the Runge-Kutta method.
  allocate(u0(nnodes,4))

! These parameters are set in main.f90. Here just print them on display.
  write(*,*)
  write(*,*) "Calling the Euler solver..."
  write(*,*)
  write(*,*) "                  M_inf = ", M_inf
  write(*,*) "                    CFL = ", CFL
  write(*,*) "             final time = ", t_final
  write(*,*) "          time_step_max = ", time_step_max
  write(*,*) "          inviscid_flux = ", trim(inviscid_flux)
  write(*,*) "           limiter_type = ", trim(limiter_type)
  write(*,*)

!--------------------------------------------------------------------------------
! First, make sure that normal mass flux is zero at all solid boundary nodes.
! NOTE: Necessary because initial solution may generate the normal component.
!--------------------------------------------------------------------------------

    call eliminate_normal_mass_flux

!--------------------------------------------------------------------------------
! Time-stepping toward the final time
!--------------------------------------------------------------------------------
  time = zero

  time_step : do i_time_step = 1, time_step_max

!------------------------------------------------------
! Two-stage Runge-Kutta scheme: u^n is saved as u0(:,:)
!  1. u^*     = u^n - (dt/vol)*Res(u^n)
!  2. u^{n+1} = 1/2*u^n + 1/2*[u^* - (dt/vol)*Res(u^*)]
!------------------------------------------------------

!-----------------------------
!- 1st Stage of Runge-Kutta:
!-----------------------------

!   Compute Res(u^n)
    call compute_residual_ncfv

!   Compute residual norms (divided residual)
    call residual_norm(res_norm)

!   Display the residual norm.
    if (i_time_step==1) write(*,'(a88)') "Density    X-momentum  Y-momentum   Energy"
    if (mod(i_time_step,10)==0) write(*,'(a2,es9.3,a9,i10,a12,4es12.2)') &
       "t=", time, "steps=", i_time_step, " L1(res)=",res_norm(:,1)

!   Stop if the final time is reached.
    if (time == t_final) exit time_step

!   Save off the previous solution (to be used in the 2nd stage).
    do i = 1, nnodes
     u0(i,:) = node(i)%u
    end do

!   Compute the time step (local and global)
    call compute_time_step(dt)
!   Adjust dt so as to finish exactly at the final time
    if (time + dt > t_final) dt = t_final - time

!   Update the solution
!   1st Stage => u^* = u^n - dt/dx*Res(u^n)
    call update_solution(one,dt,CFL)

!-----------------------------
!- 2nd Stage of Runge-Kutta:
!-----------------------------

!   Compute Res(u^*)
    call compute_residual_ncfv

!   Compute 1/2*(u^n + u^*)
    do i = 1, nnodes
     node(i)%u = half*( node(i)%u + u0(i,:) )
    end do

!   2nd Stage => u^{n+1} = 1/2*(u^n + u^*) - 1/2*dt/dx*Res(u^*)
    call update_solution(half,dt,CFL)

   time = time + dt

  end do time_step
!--------------------------------------------------------------------------------
! End of Time-stepping to the final time.
!--------------------------------------------------------------------------------

  if (i_time_step==time_step_max .and. (time /= t_final)) then
   write(*,*) " Final time not reached... Sorry..."
   write(*,*) "   Max time step reached... time_step_max=", time_step_max
   write(*,*) "   Increase time_step_max, and try again."
   write(*,*) " ... Stop"
   stop
  endif

  write(*,*) " Congratulations. Final time is reached."
  write(*,*) " Final time_step      =", i_time_step
  write(*,*) " Final time           = ", time
  write(*,*) " Final time requested = ", t_final
  write(*,*)
  write(*,*) "Finished the Euler solver... Bye!"

 end subroutine euler_solver_main
!--------------------------------------------------------------------------------
!********************************************************************************
!********************************************************************************



!********************************************************************************
!* Prepararion for Tangency condition (slip wall):
!*
!* Eliminate normal mass flux component at all solid-boundary nodes at the
!* beginning. The normal component will never be changed in the solver: the
!* residuals will be constrained to have zero normal component.
!*
!********************************************************************************
 subroutine eliminate_normal_mass_flux

 use constants   , only : p2, zero
 use my_main_data, only : node, nbound, bound

!Local variables
 integer  :: i, j, inode
 real(p2) :: normal_mass_flux, n12(2)

  bc_loop : do i = 1, nbound

   only_slip_wall : if (trim(bound(i)%bc_type) == "slip_wall") then

    write(*,*) " Eliminating the normal momentum on slip wall boundary ", i

    bnodes_slip_wall : do j = 1, bound(i)%nbnodes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS IS A SPECIAL TREATMENT FOR SHOCK DIFFRACTION PROBLEM.
!
! NOTE: This is a corner point between the inflow boundary and
!       the lower-left wall. Enforce zero y-momentum, which is
!       not ensured by the standard BCs.
!       This special treatment is necessary because the domain
!       is rectangular (the left boundary is a straight ine) and
!       the midpoint node on the left boundary is actually a corner.
!
!       Our computational domain:
!
!                 ---------------
!          Inflow |             |
!                 |             |  o: Corner node
!          .......o             |
!            Wall |             |  This node is a corner.
!                 |             |
!                 ---------------
!
!       This is to simulate the actual domain shown below:
!      
!         -----------------------
! Inflow  |                     |
!         |                     |  o: Corner node
!         --------o             |
!            Wall |             |
!                 |             |
!                 ---------------
!      In effect, we're simulating this flow by a simplified
!      rectangular domain (easier to generate the grid).
!      So, an appropriate slip BC at the corner node needs to be applied,
!      which is "zero y-momentum", and that's all.
!
     if (i==2 .and. j==1) then
                 inode = bound(i)%bnode(j)
      node(inode)%u(3) = zero               ! Make sure zero y-momentum.
      node(inode)%w    = u2w(node(inode)%u) !Update primitive variables
      cycle bnodes_slip_wall ! That's all we neeed. Go to the next.
     endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      inode = bound(i)%bnode(j)
     n12(1) = bound(i)%bnx(j)
     n12(2) = bound(i)%bny(j)

     normal_mass_flux = node(inode)%u(2)*n12(1) + node(inode)%u(3)*n12(2)

     node(inode)%u(2) = node(inode)%u(2) - normal_mass_flux*n12(1)
     node(inode)%u(3) = node(inode)%u(3) - normal_mass_flux*n12(2)

     node(inode)%w    = u2w(node(inode)%u)

    end do bnodes_slip_wall

    write(*,*) " Finished eliminating the normal momentum on slip wall boundary ", i
    write(*,*)

   endif only_slip_wall

  end do bc_loop

 end subroutine eliminate_normal_mass_flux
!--------------------------------------------------------------------------------


!********************************************************************************
!* This subroutine computes the residual for a node-centered finite-volume method
!*
!* ------------------------------------------------------------------------------
!*  Input: the current solution
!*
!* Output: node(:)%res = the residual computed by the current solution.
!* ------------------------------------------------------------------------------
!*
!* Note: dU/dt + dF/dx + dG/dy = 0. Residuals are first computed as
!*       the integral of (dF/dx + dG/dy), and at the end negative sign is added
!*       so that we have dU/dt = Res at every node.
!* 
!********************************************************************************
 subroutine compute_residual_ncfv

 use constants   , only : p2, zero, one, two, four, half, third
 use my_main_data, only : nnodes, node, nedges, edge, nbound, bound, elm, &
                          rho_inf, u_inf, v_inf, p_inf, &
                          inviscid_flux, limiter_type
 implicit none

!Local variables
 real(p2), dimension(4) :: num_flux       !Numerical flux
 real(p2), dimension(4) :: wL, wR         !Left and right face values
 real(p2), dimension(4) :: dwL, dwR       !Slope at left and right nodes
 real(p2), dimension(4) :: dwm, dwp, dwij !Re-defined slopes to be limited
 real(p2), dimension(2) :: e12            !Unit edge vector
 real(p2), dimension(2) :: n12            !Unit face normal vector
 real(p2)               :: mag_e12        !Magnitude of the edge vector
 real(p2)               :: mag_n12        !Magnitude of the face-normal vector
 real(p2)               :: wsn            !Scaled maximum wave speed
 real(p2)               :: norm_momentum  !Normal component of the momentum
 real(p2), dimension(4) :: bfluxL, bfluxR !Boundary flux at left/right nodes
 integer                :: node1, node2   !Left and right nodes of each edge
 integer                :: boundary_elm   !Element adjacent to boundary face
 integer                :: n1, n2         !Left and right nodes of boundary face
 integer                :: i, j
 integer                :: ix=1, iy=2

!-------------------------------------------------------------------------
! Gradient Reconstruction for second-order accuracy

!  Initialization

    nodes0 : do i = 1, nnodes
     node(i)%gradw = zero
    end do nodes0

!  Perform gradient reconstruction

   nodes1 : do i = 1, nnodes
     call lsq01_2x2_gradients_nc(i)
    end do nodes1

!-------------------------------------------------------------------------
! Residual computation: interior fluxes

! Initialization
  nodes : do i = 1, nnodes
   node(i)%res = zero
   node(i)%wsn = zero
  end do nodes

! Flux computation across internal edges (to be accumulated in res(:))
!
!   node2              1. Extrapolate the solutions to the edge-midpoint
!       o                 from the nodes, n1 and n2.
!        \   face      2. Compute the numerical flux
!         \ -------c2  3. Add it to the residual for n1, and subtract it from
!        / \              the residual for n2.
!   face/   \ edge
!      /     o         Directed area is the sum of the left and the right faces.
!    c1    node1       Left/right face is defined by the edge-midpoint and
!                      the centroid of the left/right element.
!                      Directed area is positive in n1 -> n2
!
! (c1, c2: element centroids)
!
!--------------------------------------------------------------------------------
  edges : do i = 1, nedges
!--------------------------------------------------------------------------------

! Left and right nodes of the i-th edge

    node1 = edge(i)%n1  ! Left node of the edge
    node2 = edge(i)%n2  ! Right node of the edge
      n12 = edge(i)%dav ! This is the directed area vector (unit vector)
  mag_n12 = edge(i)%da  ! Magnitude of the directed area vector
      e12 = edge(i)%ev  ! This is the vector along the edge (uniti vector)
  mag_e12 = edge(i)%e   ! Magnitude of the edge vector (Length of the edge)

! Solution gradient projected along the edge
!
!  NOTE: The gradient is multiplied by the distance.
!        So, it is equivalent to the solution difference.

  dwL = (node(node1)%gradw(:,ix)*e12(ix) + node(node1)%gradw(:,iy)*e12(iy) )*half*mag_e12
  dwR = (node(node2)%gradw(:,ix)*e12(ix) + node(node2)%gradw(:,iy)*e12(iy) )*half*mag_e12

!  It is now limiter time!

!  (1) No limiter (good for smooth solutions)

   limiter : if (trim(limiter_type) == "none") then

!      Simple linear extrapolation
       wL = node(node1)%w + dwL
       wR = node(node2)%w - dwR

!  (2) UMUSCL-type limiters: simple 1D limiting.

   elseif (trim(limiter_type) == "vanalbada") then

!       In 1D: dwp = w_{j+1}-w_j, dwm = w_j-w_{j-1} => limited_slope = limiter(dwm,dwp)
!
!       We can do the same in 2D as follows.
!
!       In 2D:    dwp = w_{neighbor}-w_j, dwm = 2*(grad_w_j*edge)-dwp
!              => limited_slope = limiter(dwm,dwp)
!
!      NOTE: On a regular grid, grad_w_j*edge will be the central-difference,
!            so that the average (dwm+dwp)/2 will be the central-difference just like in 1D.

!     Edge derivative
      dwij = half*(node(node2)%w - node(node1)%w)

!     Left face value (wL) with the Van Albada limiter
      dwm  = two*dwL-dwij
      dwp  = dwij
       wL  = node(node1)%w + va_slope_limiter(dwm,dwp,mag_e12)

!     Right face value (wR) with the Van Albada limiter
      dwm  = -(two*dwR-dwij)
      dwp  = -dwij
       wR  = node(node2)%w + va_slope_limiter(dwm,dwp,mag_e12)

!  (3) No other limiters available.

   else

    write(*,*) " Invalid input for limiter_type = ", trim(limiter_type)
    write(*,*) " Choose none or vanalbada, and try again."
    write(*,*) " ... Stop."
    stop

   endif limiter

!  Compute the numerical flux for given wL and wR.

!  (1) Roe flux (carbuncle is expected for strong shocks)
   if     (trim(inviscid_flux)=="roe") then

     call roe(wL,wR,n12, num_flux,wsn)

!  (2) Rotated-RHLL flux (no carbuncle is expected for strong shocks)
   elseif (trim(inviscid_flux)=="rhll") then

     call rotated_rhll(wL,wR,n12, num_flux,wsn)

   else

    write(*,*) " Invalid input for inviscid_flux = ", trim(inviscid_flux)
    write(*,*) " Choose roe or rhll, and try again."
    write(*,*) " ... Stop."
    stop

   endif

!  Add the flux multiplied by the magnitude of the directed area vector to node1,
!  and accumulate the max wave speed quantity for use in the time step calculation.

     node(node1)%res = node(node1)%res  +  num_flux * mag_n12
     node(node1)%wsn = node(node1)%wsn  +       wsn * mag_n12

!  Subtract the flux multiplied by the magnitude of the directed area vector from node2,
!  and accumulate the max wave speed quantity for use in the time step calculation.
!
!  NOTE: Subtract because the outward face normal is -n12 for the node2.

     node(node2)%res = node(node2)%res  -  num_flux * mag_n12
     node(node2)%wsn = node(node2)%wsn  +       wsn * mag_n12

!--------------------------------------------------------------------------------
  end do edges
!--------------------------------------------------------------------------------

!-------------------------------------------------------------------------
! Close with the boundary flux using the element-based formula that is
! exact for linear fluxes (See Nishikawa AIAA2010-5093 for boundary weights
! that ensure the linear exactness for 2D/3D elements).
!
!      |  Interior Domain          |
!      |        .........          |
!      |        .       .          |
!      |        .       .          |
!      o--o--o-----o---------o--o--o  <- Boundary segment
!                  n1   |   n2
!                       v
!                     n12 (unit face normal vector)
!
! NOTE: We visit each boundary face, defined by the nodes n1 and n2,
!       and compute the flux across the boundary face: left half for node1,
!       and the right half for node2. In the above figure, the dots indicate
!       the control volume around the node n1. Clearly, the flux across the
!       left half of the face contributes to the node n1. Similarly for n2.
!
!
!--------------------------------------------------------------------------------
  bc_loop : do i = 1, nbound
!--------------------------------------------------------------------------------

!------------------------------------------------
!  BC: Upwind flux via freestream values
!
!      NOTE: If the final solution at the boundary node is far from
!            the freestream values, then the domain is probably is not large enough.

   bc : if (trim(bound(i)%bc_type) == "freestream") then

    bnodes_numerical_flux_via_freestream : do j = 1, bound(i)%nbfaces

         n1 = bound(i)%bnode(j  )  !Left node
         n2 = bound(i)%bnode(j+1)  !Right node
     n12(1) = bound(i)%bfnx(j)     !x-component of the unit face normal vector
     n12(2) = bound(i)%bfny(j)     !y-component of the unit face normal vector
    mag_e12 = bound(i)%bfn(j)*half !Half length of the boundary face, j.

!   1. Left node
      wL = node(n1)%w
      wR = (/ rho_inf, u_inf, v_inf, p_inf /)
      call roe(wL,wR,n12, num_flux,wsn)
      bfluxL = num_flux
      node(n1)%wsn = node(n1)%wsn + wsn*mag_e12

!   2. Right node
      wL = node(n2)%w
      wR = (/ rho_inf, u_inf, v_inf, p_inf /)
      call roe(wL,wR,n12, num_flux,wsn)
      bfluxR = num_flux
      node(n2)%wsn = node(n2)%wsn + wsn*mag_e12

!   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
      boundary_elm = bound(i)%belm(j)

     if     (elm(boundary_elm)%nvtx == 3) then !Triangle

      node(n1)%res = node(n1)%res + (5.0_p2*bfluxL + bfluxR)/6.0_p2*mag_e12
      node(n2)%res = node(n2)%res + (5.0_p2*bfluxR + bfluxL)/6.0_p2*mag_e12

     elseif (elm(boundary_elm)%nvtx == 4) then !Quad

      node(n1)%res = node(n1)%res + bfluxL*mag_e12
      node(n2)%res = node(n2)%res + bfluxR*mag_e12

     else

      write(*,*) " Error: Element is neither tria nor quad. Stop. "; stop

     endif

    end do bnodes_numerical_flux_via_freestream

!------------------------------------------------
!  BC: Solid body and Supersonic outflow
!
!      NOTE: Basically, simply compute the physical flux, which
!            can be done by calling the Roe flux with wR = wL.
!            It is equivalent to the interior-extrapolation condition.
!
!      NOTE: Tangency condition for solid body will be applied later.

   elseif (trim(bound(i)%bc_type) == "slip_wall" .or.    &
           trim(bound(i)%bc_type) == "outflow_supersonic") then

    bnodes_slip_wall : do j = 1, bound(i)%nbfaces

         n1 = bound(i)%bnode(j  ) ! Left node
         n2 = bound(i)%bnode(j+1) ! Right node
     n12(1) = bound(i)%bfnx(j)
     n12(2) = bound(i)%bfny(j)
    mag_e12 = bound(i)%bfn(j)*half

!   1. Left node
      wL = node(n1)%w
      wR = wL
      call roe(wL,wR,n12, num_flux,wsn)
      bfluxL = num_flux
      node(n1)%wsn = node(n1)%wsn + wsn*mag_e12

!   2. Right node
      wL = node(n2)%w
      wR = wL
      call roe(wL,wR,n12, num_flux,wsn)
      bfluxR = num_flux
      node(n2)%wsn = node(n2)%wsn + wsn*mag_e12

!   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
      boundary_elm = bound(i)%belm(j)

     if     (elm(boundary_elm)%nvtx == 3) then !Triangle

      node(n1)%res = node(n1)%res + (5.0_p2*bfluxL + bfluxR)/6.0_p2*mag_e12
      node(n2)%res = node(n2)%res + (5.0_p2*bfluxR + bfluxL)/6.0_p2*mag_e12

     elseif (elm(boundary_elm)%nvtx == 4) then !Quad

      node(n1)%res = node(n1)%res + bfluxL*mag_e12
      node(n2)%res = node(n2)%res + bfluxR*mag_e12

     else

      write(*,*) " Error: Element is neither tria nor quad. Stop. "; stop

     endif

    end do bnodes_slip_wall


!------------------------------------------------
!  BC: Subsonic Outflow - Fixed Back Pressure
!
!      NOTE: Fix the pressure as freestream pressure
!            on the right side of the face (outside the domain).
!            Assumption is that the outflow boundary is far from the body.

   elseif (trim(bound(i)%bc_type) == "outflow_back_pressure") then

    bnodes_outflow : do j = 1, bound(i)%nbfaces

         n1 = bound(i)%bnode(j  ) ! Left node
         n2 = bound(i)%bnode(j+1) ! Right node
     n12(1) = bound(i)%bfnx(j)
     n12(2) = bound(i)%bfny(j)
    mag_e12 = bound(i)%bfn(j)*half

!   1. Left node
      wL    = node(n1)%w
      wR    = wL
      wR(4) = p_inf !Fix the pressure
      call roe(wL,wR,n12, num_flux,wsn)
      bfluxL = num_flux
      node(n1)%wsn = node(n1)%wsn + wsn*mag_e12

!   2. Right node
      wL = node(n2)%w
      wR = wL
      wR(4) = p_inf !Fix the pressure
      call roe(wL,wR,n12, num_flux,wsn)
      bfluxR = num_flux
      node(n2)%wsn = node(n2)%wsn + wsn*mag_e12

!   3. Add contributions to the two nodes (See Nishikawa AIAA2010-5093)
      boundary_elm = bound(i)%belm(j)

     if     (elm(boundary_elm)%nvtx == 3) then !Triangle

      node(n1)%res = node(n1)%res + (5.0_p2*bfluxL + bfluxR)/6.0_p2*mag_e12
      node(n2)%res = node(n2)%res + (5.0_p2*bfluxR + bfluxL)/6.0_p2*mag_e12

     elseif (elm(boundary_elm)%nvtx == 4) then !Quad

      node(n1)%res = node(n1)%res + bfluxL*mag_e12
      node(n2)%res = node(n2)%res + bfluxR*mag_e12

     else

      write(*,*) " Error: Element is neither tria nor quad. Stop. "; stop

     endif

    end do bnodes_outflow

!------------------------------------------------

   endif bc

!--------------------------------------------------------------------------------
  end do bc_loop
!--------------------------------------------------------------------------------


!------------------------------------------------
!  BC: Solid body - Slip condition (Tangency condition).
!
!  NOTE: It is good to enforce this condition after everything else has been done.

  bc_loop2 : do i = 1, nbound

   only_slip_wall : if (trim(bound(i)%bc_type) == "slip_wall") then

    bnodes_slip_wall2 : do j = 1, bound(i)%nbnodes

         n1 = bound(i)%bnode(j)
     n12(1) = bound(i)%bnx(j)
     n12(2) = bound(i)%bny(j)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! THIS IS A SPECIAL TREATMENT FOR SHOCK DIFFRACTION PROBLEM.
! Same as in the subroutine "eliminate_normal_mass_flux" above.

  if (i==2 .and. j==1) then
   node(n1)%res(3) = zero  ! Make sure no updates to y-momentum.
   cycle bnodes_slip_wall2 ! That's all we neeed. Go to the next.
  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    Subtract the normal component of the mass flow for tangency, so that
!    normal mass flow through the boundary will not be created at nodes
!    after the solution update.

       norm_momentum = node(n1)%res(2)*n12(1) + node(n1)%res(3)*n12(2)
     node(n1)%res(2) = node(n1)%res(2) - norm_momentum*n12(1)
     node(n1)%res(3) = node(n1)%res(3) - norm_momentum*n12(2)

    end do bnodes_slip_wall2

   endif only_slip_wall

  end do bc_loop2

! Switch the residual sign.

  nodes3 : do i = 1, nnodes

   node(i)%res = - node(i)%res

  end do nodes3

 end subroutine compute_residual_ncfv
!--------------------------------------------------------------------------------


!********************************************************************************
!* This subroutine updates the solution.
!*
!* ------------------------------------------------------------------------------
!*  Input:       coeff = coefficient for RK time-stepping
!*                  dt = global time step (not used if local time stepping)
!*                 CFL = CFL number
!*          node(:)res = the residual
!*
!* Output:   node(:)u  = updated conservative variables 
!*           node(:)w  = updated primitive variables 
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine update_solution(coeff,dt, CFL)

 use constants   , only : p2, zero
 use my_main_data, only : nnodes, node, M_inf

 implicit none

 real(kind=p2), intent(in) :: coeff, dt, CFL
!Local variables
 integer :: i

!--------------------------------------------------------------------------------
  nodes : do i = 1, nnodes
!--------------------------------------------------------------------------------

!   Solution change based on the global time stepping 

     node(i)%du = (CFL*coeff*dt/node(i)%vol) * node(i)%res

!   Solution update

     node(i)%u  = node(i)%u + node(i)%du !Make changes
     node(i)%w  = u2w(node(i)%u)         !Update primitive variables

!--------------------------------------------------------------------------------
  end do nodes
!--------------------------------------------------------------------------------

! Check the density/pressure positivity.
!
!  NOTE: Well, let's just print the warning message.
!        We may stop immediately, though.

  nodes2 : do i = 1, nnodes

!  Check density
   if (node(i)%w(1) <= zero) then
    write(*,*) " Negative density detected in update(): rho = ", &
               node(i)%w(1), i, nnodes
    node(i)%w(1) = min(0.1_p2*M_inf**2, 1.0e-04_p2)
      node(i)%u  = w2u(node(i)%w)
   endif

!  Check pressure
   if (node(i)%w(4) <= zero) then
    write(*,*) " Negative pressure detected in update(): p = ", &
               node(i)%w(4), i, nnodes
    node(i)%w(4) = min(0.1_p2*M_inf**2, 1.0e-04_p2)
      node(i)%u  = w2u(node(i)%w)
   endif

  end do nodes2


 end subroutine update_solution
!--------------------------------------------------------------------------------


!********************************************************************************
!* -- vanAlbada Slope Limiter Function--
!*
!* 'A comparative study of computational methods in cosmic gas dynamics', 
!* Van Albada, G D, B. Van Leer and W. W. Roberts, Astronomy and Astrophysics,
!* 108, p76, 1982
!*
!* ------------------------------------------------------------------------------
!*  Input:   da, db     : two differences
!*
!* Output:   va_limiter : limited difference
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
  pure function va_slope_limiter(da, db, h)

  use constants   , only : p2, two, half, one
  real(p2), dimension(4), intent(in) :: da, db
  real(p2),               intent(in) :: h
  real(p2), dimension(4)             :: va_slope_limiter
  real(p2) :: eps2

  continue

         eps2 = (0.3_p2*h)**3

    va_slope_limiter = half*( sign(one,da*db) + one ) * &
           ( (db**2 + eps2)*da + (da**2 + eps2)*db )/(da**2 + db**2 + two*eps2)


  end function va_slope_limiter


!********************************************************************************
!* -- Roe's Flux Function with entropy fix---
!*
!* P. L. Roe, Approximate Riemann Solvers, Parameter Vectors and Difference
!* Schemes, Journal of Computational Physics, 43, pp. 357-372.
!*
!* NOTE: 3D version of this subroutine is available for download at
!*       http://cfdbooks.com/cfdcodes.html
!*
!* ------------------------------------------------------------------------------
!*  Input:   primL(1:5) =  left state (rhoL, uL, vL, pL)
!*           primR(1:5) = right state (rhoR, uR, vR, pR)
!*               njk(2) = Face normal (L -> R). Must be a unit vector.
!*
!* Output:    flux(1:5) = numerical flux
!*                  wsn = half the max wave speed
!*                        (to be used for time step calculations)
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine roe(primL, primR, njk,  flux, wsn)

 use constants   , only : p2, zero, one, two, half, fifth
 use my_main_data, only : gamma

 implicit none

!Input:
 real(p2), intent( in) :: primL(4), primR(4) ! [rho, u, v, p]_{L,R}
 real(p2), intent( in) :: njk(2)             ! Face normal, njk=[nx, ny]

!Output
 real(p2), intent(out) :: flux(4)
 real(p2), intent(out) :: wsn

!Local variables
 real(p2) :: nx, ny                  ! Normal vector
 real(p2) :: mx, my                  ! Tangent vector: mx*nx+my*ny = 0
 real(p2) :: uL, uR, vL, vR          ! Velocity components.
 real(p2) :: rhoL, rhoR, pL, pR      ! Primitive variables.
 real(p2) :: unL, unR, umL, umR      ! Normal and tangent velocities
 real(p2) :: aL, aR, HL, HR          ! Speeds of sound.
 real(p2) :: RT,rho,u,v,H,a,un, um   ! Roe-averages
 real(p2) :: drho,dun,dum,dp,LdU(4)  ! Wave strenghs
 real(p2) :: ws(4), Rv(4,4)          ! Wave speeds and right-eigevectors
 real(p2) :: fL(4), fR(4), diss(4)   ! Fluxes ad dissipation term
 real(p2) :: dws(4)                  ! User-specified width for entropy fix
 integer :: i, j

  nx = njk(1)
  ny = njk(2)

!Tangent vector (Do you like it? Actually, Roe flux can be implemented 
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  mx = -ny
  my =  nx

!Primitive and other variables.
!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
     unL = uL*nx+vL*ny
     umL = uL*mx+vL*my
      pL = primL(4)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
     unR = uR*nx+vR*ny
     umR = uR*mx+vR*my
      pR = primR(4)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR)

!First compute the Roe Averages
    RT = sqrt(rhoR/rhoL)
   rho = RT*rhoL
     u = (uL+RT*uR)/(one+RT)
     v = (vL+RT*vR)/(one+RT)
     H = (HL+RT* HR)/(one+RT)
     a = sqrt( (gamma-one)*(H-half*(u*u+v*v)) )
    un = u*nx+v*ny
    um = u*mx+v*my

!Wave Strengths
   drho = rhoR - rhoL 
     dp =   pR - pL
    dun =  unR - unL
    dum =  umR - umL

  LdU(1) = (dp - rho*a*dun )/(two*a*a)
  LdU(2) = rho*dum
  LdU(3) =  drho - dp/(a*a)
  LdU(4) = (dp + rho*a*dun )/(two*a*a)

!Wave Speed
  ws(1) = abs(un-a)
  ws(2) = abs(un)
  ws(3) = abs(un)
  ws(4) = abs(un+a)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
! only for the nonlinear fields.
  dws(1) = fifth
   if ( ws(1) < dws(1) ) ws(1) = half * ( ws(1)*ws(1)/dws(1)+dws(1) )
  dws(4) = fifth
   if ( ws(4) < dws(4) ) ws(4) = half * ( ws(4)*ws(4)/dws(4)+dws(4) )

!Right Eigenvectors
  Rv(1,1) = one    
  Rv(2,1) = u - a*nx
  Rv(3,1) = v - a*ny
  Rv(4,1) = H - un*a

  Rv(1,2) = zero
  Rv(2,2) = mx
  Rv(3,2) = my
  Rv(4,2) = um

  Rv(1,3) = one
  Rv(2,3) = u
  Rv(3,3) = v 
  Rv(4,3) = half*(u*u+v*v)

  Rv(1,4) = one
  Rv(2,4) = u + a*nx
  Rv(3,4) = v + a*ny
  Rv(4,4) = H + un*a

!Dissipation Term
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*LdU(j)*Rv(i,j)
   end do
  end do

!Compute the flux.
  fL(1) = rhoL*unL
  fL(2) = rhoL*unL * uL + pL*nx
  fL(3) = rhoL*unL * vL + pL*ny
  fL(4) = rhoL*unL * HL

  fR(1) = rhoR*unR
  fR(2) = rhoR*unR * uR + pR*nx
  fR(3) = rhoR*unR * vR + pR*ny
  fR(4) = rhoR*unR * HR

  flux = half * (fL + fR - diss)
  wsn = half*(abs(un) + a)  !Normal max wave speed times half

 end subroutine roe
!--------------------------------------------------------------------------------

!*****************************************************************************
!* -- Rotated-RHLL Flux Function ---
!*
!* H. Nishikawa and K. Kitamura, Very Simple, Carbuncle-Free, Boundary-Layer
!* Resolving, Rotated-Hybrid Riemann Solvers,
!* Journal of Computational Physics, 227, pp. 2560-2581, 2008.
!*
!* Robust Riemann solver for nonlinear instability (carbuncle).
!*
!* NOTE: 3D version of this subroutine is available for download at
!*       http://cfdbooks.com/cfdcodes.html
!*
!* ------------------------------------------------------------------------------
!*  Input:   primL(1:5) =  left state (rhoL, uL, vL, pL)
!*           primR(1:5) = right state (rhoR, uR, vR, pR)
!*               njk(2) = Face normal (L -> R). Must be a unit vector.
!*
!* Output:    flux(1:5) = numerical flux
!*                  wsn = half the max wave speed
!*                        (to be used for time step calculations)
!* ------------------------------------------------------------------------------
!*
!*****************************************************************************
 subroutine rotated_rhll(primL, primR, njk,  flux, wsn)

 use constants   , only : p2, zero, one, two, half, fifth
 use my_main_data, only : gamma, M_inf

 implicit none

!Input:
 real(p2), intent( in) :: primL(4), primR(4) ! [rho, u, v, p]_{L,R}
 real(p2), intent( in) :: njk(2)             ! Face normal, njk=[nx, ny]

!Output
 real(p2), intent(out) :: flux(4)
 real(p2), intent(out) :: wsn

!Local variables
 real(p2) :: nx, ny                  ! Normal vector
 real(p2) :: mx, my                  ! Tangent vector: mx*nx+my*ny = 0
 real(p2) :: uL, uR, vL, vR          ! Velocity components.
 real(p2) :: rhoL, rhoR, pL, pR      ! Primitive variables.
 real(p2) :: unL, unR, umL, umR      ! Normal and tangent velocities
 real(p2) :: aL, aR, HL, HR          ! Speeds of sound.
 real(p2) :: RT,rho,u,v,H,a,un, um   ! Roe-averages
 real(p2) :: drho,dun,dum,dp,LdU(4)  ! Wave strenghs
 real(p2) :: ws(4), Rv(9,9)          ! Wave speeds and right-eigevectors
 real(p2) :: abs_ws(4)
 real(p2) :: fL(4), fR(4), diss(4)   ! Fluxes ad dissipation term
 real(p2) :: dws(4)                  ! User-specified width for entropy fix
 integer  :: i, j
 real(p2) :: eps

 real(p2) :: nx1, ny1, nx2, ny2             ! Rotated normals, n1 and n2
 real(p2) :: alpha1, alpha2                 ! Projections of the new normals
 real(p2) :: abs_dq                         ! Magnitude of the velocity difference
 real(p2) :: SRp,SLm                        ! Wave speeds for the HLL part
 real(p2) :: temp

  nx = njk(1)
  ny = njk(2)

!Tangent vector (Do you like it? Actually, Roe flux can be implemented 
! without any tangent vector. See "I do like CFD, VOL.1" for details.)
  mx = -ny
  my =  nx

!Primitive and other variables.
!  Left state
    rhoL = primL(1)
      uL = primL(2)
      vL = primL(3)
     unL = uL*nx+vL*ny
     umL = uL*mx+vL*my
      pL = primL(4)
      aL = sqrt(gamma*pL/rhoL)
      HL = aL*aL/(gamma-one) + half*(uL*uL+vL*vL)
!  Right state
    rhoR = primR(1)
      uR = primR(2)
      vR = primR(3)
     unR = uR*nx+vR*ny
     umR = uR*mx+vR*my
      pR = primR(4)
      aR = sqrt(gamma*pR/rhoR)
      HR = aR*aR/(gamma-one) + half*(uR*uR+vR*vR)

!Primitive and oth
!Compute the flux.
   fL(1) = rhoL*unL
   fL(2) = rhoL*unL * uL + pL*nx
   fL(3) = rhoL*unL * vL + pL*ny
   fL(4) = rhoL*unL * HL

   fR(1) = rhoR*unR
   fR(2) = rhoR*unR * uR + pR*nx
   fR(3) = rhoR*unR * vR + pR*ny
   fR(4) = rhoR*unR * HR

!Define n1 and n2, and compute alpha1 and alpha2: (4.2) in the original paper.
!(NB: n1 and n2 may need to be frozen at some point during 
!     a steady calculation to fully make it converge. For time-accurate 
!     calculation, this is fine.)
! NB: For a boundary face, set (nx2,ny2)=(nx,ny), (nx1,ny1)=(-ny,nx).

    eps = 1.0e-12_p2 * M_inf
    abs_dq = sqrt( (uR-uL)**2 + (vR-vL)**2 )

  if ( abs_dq > eps) then
    nx1 = (uR-uL)/abs_dq
    ny1 = (vR-vL)/abs_dq
  else
    nx1 = -ny 
    ny1 =  nx
  endif

!  Rey = 1000.0_p2
!  temp = ( tanh(Rey*(abs_dq-eps)) - tanh(-Rey*eps) ) &
!        /( tanh(Rey*(   one-eps)) - tanh(-Rey*eps) )
!  nx1 = temp*(uR-uL)/(abs_dq + eps) + (one-temp)*(-ny)
!  ny1 = temp*(vR-vL)/(abs_dq + eps) + (one-temp)*( nx)

    alpha1 = nx * nx1 + ny * ny1 
!   To make alpha1 always positive.
      temp = sign(one,alpha1)
       nx1 = temp * nx1
       ny1 = temp * ny1
    alpha1 = temp * alpha1

! Take n2 as perpendicular to n1.
       nx2 = -ny1
       ny2 =  nx1
    alpha2 = nx * nx2 + ny * ny2
!   To make alpha2 always positive.
      temp = sign(one,alpha2)
       nx2 = temp * nx2
       ny2 = temp * ny2
    alpha2 = temp * alpha2

!Now we are going to compute the Roe flux with n2 as the normal
!and n1 as the tagent vector, with modified wave speeds (5.12)

     RT = sqrt(rhoR/rhoL)
    rho = RT*rhoL
      u = (uL + RT*uR)/(one + RT)
      v = (vL + RT*vR)/(one + RT)
      H = (HL + RT*HR)/(one + RT)
      a = sqrt( (gamma-one)*(H-half*(u*u+v*v)) )
     un = u*nx2+v*ny2
     um = u*nx1+v*ny1

!Wave Strengths (remember that n2 is the normal and n1 is the tangent.)
    unL = uL*nx2 + vL*ny2
    unR = uR*nx2 + vR*ny2
    umL = uL*nx1 + vL*ny1
    umR = uR*nx1 + vR*ny1

   drho = rhoR - rhoL 
     dp =   pR - pL
    dun =  unR - unL
    dum =  umR - umL

  LdU(1) = (dp - rho*a*dun )/(two*a*a)
  LdU(2) =  rho*dum
  LdU(3) =  drho - dp/(a*a)
  LdU(4) = (dp + rho*a*dun )/(two*a*a)

!Wave Speeds for Roe flux part.
    ws(1) = un-a
    ws(2) = un
    ws(3) = un
    ws(4) = un+a
  abs_ws  = abs(ws)

!Harten's Entropy Fix JCP(1983), 49, pp357-393:
!only for the nonlinear fields.
  dws(1) = fifth
   if (abs_ws(1)<dws(1)) abs_ws(1) = half*(abs_ws(1)*abs_ws(1)/dws(1)+dws(1))
  dws(4) = fifth
   if (abs_ws(4)<dws(4)) abs_ws(4) = half*(abs_ws(4)*abs_ws(4)/dws(4)+dws(4))

!HLL wave speeds, evaluated with [nx1,ny1] (=tangent wrt n2).
   SRp = max( zero, umR + aR, um + a)
   SLm = min( zero, umL - aL, um - a)

!Modified wave speeds for the Rotated-RHLL flux: (5.12) in the original paper.
   ws = alpha2*abs_ws - ( alpha2*(SRp+SLm)*ws + two*alpha1*SRp*SLm )/ (SRp-SLm)

!Right Eigenvectors: with n2 as normal and n1 as tangent.
  mx = nx1
  my = ny1

  Rv(1,1) = one    
  Rv(2,1) = u - a*nx2
  Rv(3,1) = v - a*ny2
  Rv(4,1) = H - a*un

  Rv(1,2) = zero
  Rv(2,2) = mx
  Rv(3,2) = my
  Rv(4,2) = um

  Rv(1,3) = one
  Rv(2,3) = u
  Rv(3,3) = v 
  Rv(4,3) = half*(u*u+v*v)

  Rv(1,4) = one
  Rv(2,4) = u + a*nx2
  Rv(3,4) = v + a*ny2
  Rv(4,4) = H + a*un

!Dissipation Term: Roe dissipation with the modified wave speeds.
  diss = zero
  do i=1,4
   do j=1,4
    diss(i) = diss(i) + ws(j)*LdU(j)*Rv(i,j)
   end do
  end do

!Compute the Rotated-RHLL flux.
  flux = (SRp*fL - SLm*fR)/(SRp-SLm) - half*diss

  wsn = half*( abs(un) + abs(um) + a)  !Normal max wave speed times half

 end subroutine rotated_rhll

!********************************************************************************
!* This subroutine computes the explicit time-step: the minimum dt over nodes.
!*
!* ------------------------------------------------------------------------------
!*  Input: node(i)%vol = Dual volume
!*         node(i)%wsn = Sum of the max wave speed multiplied by the face length
!*
!* Output:         dt  = global time step
!*         node(:)%dt  =  local time step
!* ------------------------------------------------------------------------------
!*
!* NOTE: Local time step is computed and stored at every node, but not used.
!*       For steady problems, it can be used to accelerate the convergence.
!*
!********************************************************************************
 subroutine compute_time_step(dt)

 use constants   , only : p2, half, one, two
 use my_main_data, only : nnodes, node

 implicit none

 real(p2), intent(out) ::  dt
!Local variables
 integer  :: i
 real(p2) :: dt_min

 dt_min = 1.0e+05_p2

!--------------------------------------------------------------------------------
  nodes : do i = 1, nnodes
!--------------------------------------------------------------------------------

! Local time step: dt = volume/sum(0.5*max_wave_speed*face_area).

    node(i)%dt = node(i)%vol / node(i)%wsn

! Keep the minimum dt

   if (i==1) dt_min = node(i)%dt
   dt_min = min( dt_min, node(i)%dt )

!--------------------------------------------------------------------------------
  end do nodes
!--------------------------------------------------------------------------------

! Global time-step

   dt = dt_min

 end subroutine compute_time_step
!--------------------------------------------------------------------------------

!********************************************************************************
!* This subroutine computes the residual norms: L1, L2, L_infty
!*
!* ------------------------------------------------------------------------------
!*  Input:  node(:)res = the residuals
!*
!* Output:  res_norm   = residual norms (L1, L2, Linf)
!* ------------------------------------------------------------------------------
!*
!* NOTE: It is not done here, but I advise you to keep the location of the
!*       maximum residual (L_inf).
!*
!********************************************************************************
 subroutine residual_norm(res_norm)

 use constants   , only : p2, zero, one
 use my_main_data, only : node, nnodes

 implicit none

 real(kind=p2), dimension(4,3), intent(out)   :: res_norm

!Local variables
 real(kind=p2), dimension(4) :: residual
 integer :: i

   res_norm(:,1) =  zero
   res_norm(:,2) =  zero
   res_norm(:,3) = - one

!--------------------------------------------------------------------------------
  nodes : do i = 1, nnodes
!--------------------------------------------------------------------------------
   residual = abs( node(i)%res/node(i)%vol )      !Divided residual
   res_norm(:,1) = res_norm(:,1)    + residual    !L1   norm
   res_norm(:,2) = res_norm(:,2)    + residual**2 !L2   norm
   res_norm(:,3) = max(res_norm(:,3), residual)   !Linf norm
!--------------------------------------------------------------------------------
  end do nodes
!--------------------------------------------------------------------------------

   res_norm(:,1) =      res_norm(:,1)/real(nnodes,p2)
   res_norm(:,2) = sqrt(res_norm(:,2)/real(nnodes,p2))

 end subroutine residual_norm
!--------------------------------------------------------------------------------


!********************************************************************************
!* Given the primitive variables at nodes {k} around node i, this computes the 
!* gradients (wx,wy) at node j by the unweighted least-squares method.
!*
!* ------------------------------------------------------------------------------
!*  Input:  node(inode)%w     = current primitive variables at node 'inode'
!*
!* Output:  node(inode)%gradw = gradients of the primitive variables
!* ------------------------------------------------------------------------------
!*  
!*  1. At node i, compute a vector b = \sum_k [ (xk-xi)*(uk-ui), (yk-yi)*(uk-ui) ]
!*  2. Compute the gradient by multiplying b by the inverse LSQ matrix that has
!*     been pre-computed by the subroutine lsq01_matrix() in the main:
!*             ux = inverse(1,1)*b(1) + inverse(1,2)*b(2)
!*             uy = inverse(2,1)*b(1) + inverse(2,2)*b(2)
!*
!********************************************************************************
 subroutine lsq01_2x2_gradients_nc(inode)

 use my_main_data, only : node
 use constants   , only : p2, zero

 implicit none

 integer, intent(in) :: inode

!Local variables
 integer  :: k, in, inghbr
 integer  :: ix, iy
 real(p2) :: dvar, dx, dy, b(2)

  ix=1 
  iy=2

! Loop over variables

    do k = 1, 4

       b = zero

!   Loop over neighbors

     do in = 1, node(inode)%nnghbrs
       inghbr = node(inode)%nghbr(in)

      dvar = node(inghbr)%w(k) - node(inode)%w(k) !primitive variables
        dx = node(inghbr)%x    - node(inode)%x
        dy = node(inghbr)%y    - node(inode)%y

        b(1) = b(1) + dvar*dx
        b(2) = b(2) + dvar*dy

     end do

 !  Multiply the inverse LSQ matrix to get the gradients

      node(inode)%gradw(k,ix) = node(inode)%lsq01_2x2_minv(ix,1)*b(1) &
                              + node(inode)%lsq01_2x2_minv(ix,2)*b(2)

      node(inode)%gradw(k,iy) = node(inode)%lsq01_2x2_minv(iy,1)*b(1) &
                              + node(inode)%lsq01_2x2_minv(iy,2)*b(2)

    end do

 end subroutine lsq01_2x2_gradients_nc
!--------------------------------------------------------------------------------

!********************************************************************************
!* --- Inverse Matrix for 2x2 Least-Squares Gradient Reconstruction ---
!*
!* Construct a matrix for the linear least-squares(LSQ) gradient reconstruction.
!* (unweighted LSQ; more accurate than weighted ones to my knowledge.)
!*
!* Note: it requires at least 2 non-colinear neighbors.
!*
!* Example: Consider constructing (ux,uy) at i with the following stencil.
!*
!*      3 o     o 2
!*         \   / 
!*          \ /
!*         i *-----o 1
!*          /|
!*         / |
!*        /  o 5      *: node in interest (i)
!*       o 4          o: neighbors (k = 1,2,3,4,5)
!*
!*  5 equations:
!*    (x1-xi)*ux + (y1-yi)*uy = (u1-ui)
!*    (x2-xi)*ux + (y2-yi)*uy = (u2-ui)
!*    (x3-xi)*ux + (y3-yi)*uy = (u3-ui)
!*    (x4-xi)*ux + (y4-yi)*uy = (u4-ui)
!*    (x5-xi)*ux + (y5-yi)*uy = (u5-ui)
!*
!*  This system is written in the matrix form:
!*
!*        A*x = b,  x=(ux,uy), A=5x2 matrix, b=5x1 matrix
!*
!*  The least-squares problem is
!*
!*      A^T*A*x = A^T*b, (T denotes the transpose: A^T=2x5 matrix)
!*  
!*  which is
!*
!*  [sum_k (xk-xi)^2]*ux       + [sum_k (xk-xi)*(yk-yi)]*uy = [sum_k (uk-ui)*(xk-xi)]
!*  [sum_k (xk-xi)*(yk-yi)]*ux + [sum_k (yk-yi)]*uy         = [sum_k (uk-ui)*(yk-yi)]
!*
!* This subroutine computes the inverse of (A^T*A) at every node (which depends
!* only on the grid geometry), so that the gradient at a node can be computed
!* by a matrix-vector multiplication, i.e., (A^T*A)^{-1}*(A^T*b), 
!* (only A^T*b needs to be re-computed).
!*
!* ------------------------------------------------------------------------------
!*  Input:  inode = node number
!*
!* Output:  node(inode)%lsq01_2x2_minv = inverse matrix for LSQ reconstruction
!* ------------------------------------------------------------------------------
!*
!********************************************************************************
 subroutine lsq01_2x2_matrix_nc(inode)

 use my_main_data, only : node
 use constants   , only : p2, zero, half, third, fourth

 implicit none

 integer, intent(in) :: inode
!Local variables
 real(p2) :: a(2,2), dx, dy, det
 integer :: k, inghbr

   a = zero

!  Loop over the neighbor nodes.
   do k = 1, node(inode)%nnghbrs
    inghbr = node(inode)%nghbr(k)

      dx = node(inghbr)%x - node(inode)%x
      dy = node(inghbr)%y - node(inode)%y

      a(1,1) = a(1,1) + dx*dx
      a(1,2) = a(1,2) + dx*dy

      a(2,1) = a(2,1) + dx*dy
      a(2,2) = a(2,2) + dy*dy

   end do

    det = a(1,1)*a(2,2) - a(1,2)*a(2,1)
    if (abs(det) < 1.0e-14) write(*,*) " Singular: LSQ det = ", det, " i=",inode

! OK, invert and store the inverse matrix:
     node(inode)%lsq01_2x2_minv(1,1) =  a(2,2)/det
     node(inode)%lsq01_2x2_minv(1,2) = -a(2,1)/det
     node(inode)%lsq01_2x2_minv(2,1) = -a(1,2)/det
     node(inode)%lsq01_2x2_minv(2,2) =  a(1,1)/det

 end subroutine lsq01_2x2_matrix_nc
!********************************************************************************
!*
!********************************************************************************

!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  w =    primitive variables (rho,     u,     v,     p)
!* Output:  u = conservative variables (rho, rho*u, rho*v, rho*E)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
 function w2u(w) result(u)

 use constants   , only : p2, one, half
 use my_main_data, only : gamma

 implicit none

 real(p2), dimension(4), intent(in) :: w

!Local variables
 real(p2), dimension(4)             :: u !output

  u(1) = w(1)
  u(2) = w(1)*w(2)
  u(3) = w(1)*w(3)
  u(4) = w(4)/(gamma-one)+half*w(1)*(w(2)*w(2)+w(3)*w(3))

 end function w2u
!--------------------------------------------------------------------------------

!********************************************************************************
!* Compute U from W
!*
!* ------------------------------------------------------------------------------
!*  Input:  u = conservative variables (rho, rho*u, rho*v, rho*E)
!* Output:  w =    primitive variables (rho,     u,     v,     p)
!* ------------------------------------------------------------------------------
!* 
!********************************************************************************
 function u2w(u) result(w)

 use constants   , only : p2, one, half
 use my_main_data, only : gamma

 implicit none

 real(p2), dimension(4), intent(in) :: u

!Local variables
 real(p2), dimension(4)             :: w !output

  w(1) = u(1)
  w(2) = u(2)/u(1)
  w(3) = u(3)/u(1)
  w(4) = (gamma-one)*( u(4) - half*w(1)*(w(2)*w(2)+w(3)*w(3)) )

 end function u2w
!--------------------------------------------------------------------------------


 end module euler_solver
