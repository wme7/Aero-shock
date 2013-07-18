subroutine Compute_time_step
  use MD2D_Grid
  use RK_Var
  use Material_Var
  use Metric_Var
  implicit none
  !
  !---- Purpose: Calculate time step for the linear Maxwell system
  !	       normalized to unit speed of light.
  !
  !     Input  :	-
  !
  !     Output : dt     New time step - dt
  !
  !
  !
  !.... Declare local parameters
  integer:: i, j, k
  real(kind=8) :: dt_max
  real(kind=8) :: local_speed

  !
  !.... Loop over all domains
  dt_max=0.d0

  do DDK=1,TotNum_DM

     local_speed=dsqrt( (Lame_lambda(DDK)+2.d0*Lame_mu(DDK))/rho(DDK) )
     !write(*,*) local_speed
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)


     !... Calculate local maximum inverse timestep
     
     do j=0,ND2
        do i=0,ND1
           dt_max=max(dt_max,local_speed*dtrans(i,j,DDK))
        end do
     end do

  enddo
!
!.... Calculate maximum global inverse time step
  dt=CFL/dt_max
!      dt=CFL*dt_max
  return
end subroutine Compute_time_step
!
!======================================================================
