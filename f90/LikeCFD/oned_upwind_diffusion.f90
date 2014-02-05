!*******************************************************************************
!* One-dimensional O(h)-Time-Step Upwind Diffusion Scheme.
!*
!*  Problem: u_t = nu*u_{xx} + nu*pi^2*sin(pi*x), u(0)=u(1)=0,
!*
!*           which is solved in the first-order system form:
!*
!*           u_t = nu*p_x + nu*pi^2*sin(pi*x)
!*           p_t = (u_x - p)/Tr
!*
!*     This scheme is O(1/h) times faster than conventional diffusion schemes.
!*
!* Note: The first-order system is hyperbolic. The eigenvalues (wave speeds) are
!*       + nu/Lr, and - nu/Lr.
!*
!* Note: The first-order system is equivalent to the original diffusion equation
!*       in the steady state. The idea is to integrate the first-order system in
!*       time, instead of the original diffusion equation. Thereby, the time step
!*       is O(h), not O(h^2), and accurate solution gradient can be computed
!*       simultaneously. See Nishikawa-JCP2007 for details.
!*
!* Note: The upwind scheme employed here is a residual-distribution scheme.
!*       Other methods such as the finite-volume scheme can be employed instead.
!*       Just solve the first-order system by any method, you'll have O(h) time
!*       step and equally accurate solution gradient.
!*
!* Note: The program now has a common finite-difference scheme.
!*       Compare the performance of the upwind scheme and the common scheme for
!*       a large number of nodes (say, 512). The upwind scheme gets faster and
!*       faster as the grid gets finer. Remember the speed-up factor increases
!*       as O(1/h).
!*
!* H. Nishikawa, A First-Order System Approach for Diffusion Equation. I: 
!* Second-Order Residual Distribution Schemes, Journal of Computational Physics, 
!* 227, pp. 315-352, 2007
!* (All equation numbers refer to the above paper.)
!*
!* This program is written by Hiroaki Nishikawa, 10-31-2009
!* Revised on April 16, 2011
!* Revised on April 27, 2011 (Common scheme added for comparison)
!*******************************************************************************
 program oned_upwind_diffusion_scheme

 implicit none

  integer, parameter :: dp = selected_real_kind(15)
 real(dp), parameter :: zero=0.0_dp, one=1.0_dp, two=2.0_dp, half=0.5_dp, four=4.0_dp
 real(dp), parameter :: pi=3.141592653589793238_dp

    real(dp), dimension(:), allocatable ::      u,      p  ! Solution & Gradient
    real(dp), dimension(:), allocatable :: uexact, pexact  ! Exact solutions
    real(dp), dimension(:), allocatable ::  res_u,  res_p  ! Nodal Residuals
    real(dp) :: nu          ! Diffusion coefficient (constant)
    real(dp) :: h, xj, xjp1 ! Mesh spacing, nodal coordinates of node j and j+1
    real(dp) :: Lr, Tr      ! Length scale and relaxation time
    real(dp) :: dt          ! Time step which will be O(h), not O(h^2).
    real(dp) :: res_max     ! Maximum residual to check the convergence.
    real(dp) :: phiC(2)     ! Cell residual (fluctuation)
    real(dp) :: Bj(2,2)     ! Distribution matrix to distribute phiC to the left
    real(dp) :: Bjp1(2,2)   ! Distribution matrix to distribute phiC to the right
    real(dp) :: source_term ! Source term = nu*pi^2*sin(pi*x)
 integer :: nnodes          ! Total number of nodes
 integer :: k, j
 integer :: scheme_type     ! = 1 : Upwind scheme solving the hyperbolic system
                            ! = 2 : Common finite-difference scheme

! Initialization (values have no meaning; they will be overwritten later.)
 Lr = one
 Tr = one
 dt = one

!--------------------------------------------------------------------------------
! Diffusion coefficient

  nu = one

!--------------------------------------------------------------------------------
! Input the type of scheme

  do j = 1, 100

  write(*,*) " Type of scheme = ?"
  write(*,*) 
  write(*,*) "   1 -> Upwind scheme solving the hyperbolic system"
  write(*,*) "   2 -> Common finite-difference (Galerkin) scheme"
  read(*,*) scheme_type

  if (scheme_type /= 1 .and. scheme_type /= 2) then
   write(*,*) " Wrong number. Enter 1 or 2."
  else
   exit
  endif

  end do

!--------------------------------------------------------------------------------
! Input the number of nodes.

  write(*,*) " The number of nodes = ?"
  read(*,*) nnodes

!--------------------------------------------------------------------------------
! Allocate the arrays

   allocate(     u(nnodes),      p(nnodes)) ! Solution & Gradient
   allocate(uexact(nnodes), pexact(nnodes)) ! Exact solutions
   allocate( res_u(nnodes),  res_p(nnodes)) ! Nodal Residuals

!--------------------------------------------------------------------------------
! Set up                           h
!          Grid ->  o-----o-----o-----o-----o-----o-----o
!                  j=1          j     j+1                 j=nnodes

     h = one / real(nnodes-1,dp) ! Mesh spacing of uniform grid

   do j = 1, nnodes
           xj = real(j-1)*h      ! xj = x-coordinate of j-th node
    uexact(j) = sin(pi*xj)       ! Exact solution at xj
    pexact(j) = pi*cos(pi*xj)    ! Exact gradient at xj
         u(j) = zero             ! Initial solution
         p(j) = zero             ! Initial gradient
   end do

!--------------------------------------------------------------------------------
! Parameters

! 1. Upwind scheme:
  if (scheme_type == 1) then

    Lr = h*(one + one/sin(half*pi*h))/four ! Optimal formula for Lr: Eq.(3.57)
    Tr = Lr*Lr / nu                        ! Relaxation time: Eq.(2.28)
    dt = 0.99_dp * h/(nu/Lr)               ! Time step (CFL condition): Eq.(3.38)

!   Note: The time step is O(h), not O(h^2). The number of iterations to reach
!         the steady state will therefore be proportional to 1/h or equivalently
!         to nnodes. This is orders of magnitude faster than almost all conventional
!         diffusion schemes for which the number of iterations increases quadratically.

! 2. Common scheme
  elseif (scheme_type == 2) then

   dt = 0.99_dp * h*h/(two*nu)             ! Time step, typical O(h^2)

  endif

!--------------------------------------------------------------------------------
! Advance in time to reach the steady state

  time_loop : do k = 1, 10000000

!--- Step 1. Residual Computation (compute res_u and/or res_p):

!************************************************************************
! Option 1: Upwind scheme solving the hyperbolic system
!           Looks complicated?
!           Very expensive compared with the common scheme?
!           But it converges much much faster than the common scheme!
!************************************************************************
   scheme_choice : if (scheme_type == 1) then

!     Upwind distribution matrices, Eqs.(3.20) and (3.21):
        Bj(1,1) = half   ;    Bj(1,2) = half*Lr
        Bj(2,1) = half/Lr;    Bj(2,2) = half

      Bjp1(1,1) = half   ;  Bjp1(1,2) =-half*Lr
      Bjp1(2,1) =-half/Lr;  Bjp1(2,2) = half

!     Initialize the nodal residual arrays
      res_u = zero
      res_p = zero

!     Loop over cells: evaluate the cell-residual, phiC, and distribute it
!     to the left and right nodes in each cell.

      cell_loop : do j = 1, nnodes-1 ! j-th cell = [x(j), x(j+1)]

!      Compute the source term, nu*pi*pi*sin(pi*x), by the trapezoidal rule.
              xj   = real(j-1)*h
              xjp1 = real(j  )*h
       source_term = nu*pi*pi*half*( sin(pi*xj) + sin(pi*xjp1) )*h ! nu*pi^2*sin(pi*x)

!      Compute the cell-residual, phiC, as in Eq.(3.8):
           phiC(1) =  nu * ( p(j+1) - p(j) ) + source_term     ! nu*px + source
           phiC(2) = ( u(j+1)-u(j) - half*(p(j)+p(j+1))*h )/Tr ! (ux - p)/Tr

!      Distribute the cell-residual to the left(j) and the right(j+1) nodes.

!      (1)To the left node: j <-- Bj*phiC
        res_u(j  ) = res_u(j  ) +    Bj(1,1)*phiC(1)+  Bj(1,2)*phiC(2)
        res_p(j  ) = res_p(j  ) +    Bj(2,1)*phiC(1)+  Bj(2,2)*phiC(2)

!      (2)To the right node: Bjp1*phiC --> j+1 
        res_u(j+1) = res_u(j+1) +  Bjp1(1,1)*phiC(1)+Bjp1(1,2)*phiC(2)
        res_p(j+1) = res_p(j+1) +  Bjp1(2,1)*phiC(1)+Bjp1(2,2)*phiC(2)

      end do cell_loop

!************************************************************************
! Option 2: Common finite-difference scheme (equivalent to the Galerkin)
!           Yes, very very simple and cheap scheme, but SLOW to converge.
!************************************************************************
   else scheme_choice

    res_p = zero ! Gradient, p, is not computed.

    node_loop : do j = 2, nnodes-1 ! j-th interior node

           xj = real(j-1)*h
     res_u(j) = nu*(u(j+1)-two*u(j)+u(j-1))/h + nu*pi*pi*sin(pi*xj)*h

    end do node_loop


   end if scheme_choice

!-----------------------
!  Check convergence
!-----------------------

!     Compute the maximum nodal residual (to check the convergence)
      res_max = max( maxval(abs(res_u)), maxval(abs(res_p)) )

!     Stop if the tolerance (say, 1.0e-08) is reached.
      if ( res_max < 1.0e-08_dp ) exit time_loop

!     Display the max nodal residual at every 100 iterations
      if (mod(k,100) == 0) then
        write(*,'(a4,i10,a20,es12.5)') "Itr=", k, "   Max(nodal res) = ", res_max
      endif

!--- Step 2. Solution Updates, Eq.(3.11):
            u = u + (dt/h)*res_u
            p = p + (dt/h)*res_p

         u(1) = zero ! BC at x=0
    u(nnodes) = zero ! BC at x=1

  end do time_loop

!--------------------------------------------------------------------------------
!  Steady state is reached.

   write(*,'(a4,i10,a20,es12.5)') "Itr=", k-1, "   Max(nodal res) = ", res_max


!--------------------------------------------------------------------------------
! Display the results

! 1. Upwind scheme

  if (scheme_type == 1) then

!  L_infinity Errors:
   write(*,*)
   write(*,'(a25,es10.3)') " L_infinity(u-uexact) = ", maxval(abs(u-uexact))
   write(*,'(a25,es10.3)') " L_infinity(p-pexact) = ", maxval(abs(p-pexact))
   write(*,*)

   write(*,*) "-- Run again with nnodes = ", 2*(nnodes-1) + 1
   write(*,*) "   The errors will be about 1/4 of the above (i.e.,2nd-order accurate)."
   write(*,*) "   The number of iterations will be around 2 times the above = ",2*(k-1),")"
   write(*,*)

! Common scheme

  elseif (scheme_type == 2) then

!  L_infinity Errors:
   write(*,*)
   write(*,'(a25,es10.3)') " L_infinity(u-uexact) = ", maxval(abs(u-uexact))
   write(*,*)
   write(*,*) "-- Run again with nnodes = ", 2*(nnodes-1) + 1
   write(*,*) "   The errors will be about 1/4 of the above (i.e.,2nd-order accurate)."
   write(*,*) "   The number of iterations will be around 4 times the above = ",4*(k-1),")"
   write(*,*) "   Or try the upwind scheme for the same number of nodes (much faster!)."
   write(*,*)

  endif

 stop

 end program oned_upwind_diffusion_scheme
