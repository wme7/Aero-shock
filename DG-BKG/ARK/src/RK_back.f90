module RK_Var
  implicit none
  !parameter for RK carpenter version
  real(kind=8), save             :: CFL
  real(kind=8), save             :: Time_Final
  real(kind=8), save, allocatable:: RK_para_a(:),RK_para_b(:),RK_para_c(:)
  real(kind=8), allocatable:: const_a_E(:,:),const_a_I(:,:)
  real(kind=8), allocatable:: const_b(:),const_b_hat(:),const_c(:)
  real(kind=8), save, allocatable:: alpha(:)
  real(kind=8), save             :: dt
  real(kind=8), save             :: time, time_loc 
  integer, save                  :: ks
  real(kind=8)                   :: dt_min,dt_local
  integer, save                  :: First_Comp, Last_Comp
  integer:: RK_Stage
  INTEGER, DIMENSION (1:2) :: ORDER2 = (/ 2, 1 /)
!  real(kind=8), allocatable:: A(:,:)
!  real(kind=8), allocatable:: c(:) ,b(:)

contains 

  !------------------------------------------------------------------------------

  subroutine Input_RK_parameter
    implicit none
    
    integer:: lid
    
    lid=20
    open(lid,file='RK_Para.in',form='formatted', status='old')
    read(lid,*) !==============================================!
    read(lid,*) CFL
    read(lid,*) Time_Final
    read(lid,*) !==============================================!
    
    close(lid) 
    
    return 
    
  end subroutine Input_RK_parameter
  
  !------------------------------------------------------------------------------
  !Calculation of the coefficients alpha for SSPRK
  subroutine Init_SSPRK(rk)
    implicit none
    integer:: rk,i,j,k,ierr
    real(kind=8)::temp

    allocate(alpha(1:rk), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate rk_alpha'
       stop
    endif

    alpha(1) = 1d0
    do i=1,rk
       do j=(i-1),1,-1
          temp = dble(j)
          temp = dble(1d0/temp)
          alpha(j+1) = temp*alpha(j)
       enddo
       temp = 1.d0
       do k=1,i
          temp = temp*dble(k)
       enddo
       alpha(i) = 1d0/dble(temp)
       alpha(1) = 1d0-sum(alpha(2:i))
    enddo
  end subroutine Init_SSPRK
    
  !----------------------------------------------------------------------------
  subroutine SSPRK_marching_2D(q,rk,dtt)
  use MD2D_Grid
  use Legendre
  implicit none
    integer:: rk,q
    real(kind=8):: dtt

    if (q .lt. rk) then  
       F_new = F_new + alpha(q+1)*(F_alt+dtt*F_t)
       F_alt = F_alt + dtt*F_t
    else
       F_new = F_new + alpha(rk)*dtt*F_t
    endif
!    write(6,*) "F_t",F_t(1:3,1,0:2,0:2,3)
  end subroutine SSPRK_marching_2D
  !------------------------------------------------------------------------------
  subroutine RK4s5_marching_2D_kinetic(q,rk,dtt)
    !--Declare subroutine arguments
  use MD2D_Grid
  use Legendre

    implicit none
    integer:: rk,q
    real(kind=8):: dtt

    
     F_new = RK_para_a(q)*F_new+ dtt*F_t
     F_alt = F_alt + RK_para_b(q)*F_new

!     q_tmp(0:N1) = RK_para_a(ks)*q_tmp(0:N1) &
!                          + dt*dqdt(0:N1,0:N2)    
!     q(0:N1,0:N2) = q(0:N1,0:N2) &
!                          + RK_para_b(ks)*q_tmp(0:N1,0:N2) 
    
    
    return
    !
  end subroutine RK4s5_marching_2D_kinetic

!=====================================================================

  subroutine Init_RK4S5()
    implicit none
    
    integer:: RK_Stage
    integer i,ierr
    
    ! This is RK4S5 of 2N storage version
    RK_Stage=5
    allocate(Rk_para_a(0:RK_Stage), &
             Rk_para_b(0:RK_Stage), &
             Rk_para_c(0:RK_Stage), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Can not allocate RK_para'
       write(*,*)'Abort!'
       stop
    endif
    !
    RK_para_a(1) =              0.0d0
    RK_para_a(2) =  -567301805773.0d0/ 1357537059087.0d0
    RK_para_a(3) = -2404267990393.0d0/ 2016746695238.0d0
    RK_para_a(4) = -3550918686646.0d0/ 2091501179385.0d0
    RK_para_a(5) = -1275806237668.0d0/  842570457699.0d0
    !
    RK_para_b(1) =  1432997174477.0d0/ 9575080441755.0d0
    RK_para_b(2) =  5161836677717.0d0/13612068292357.0d0
    RK_para_b(3) =  1720146321549.0d0/ 2090206949498.0d0
    RK_para_b(4) =  3134564353537.0d0/ 4481467310338.0d0
    RK_para_b(5) =  2277821191437.0d0/14882151754819.0d0
    !
    RK_para_c(1) =              0.0d0
    RK_para_c(2) =  1432997174477.0d0/ 9575080441755.0d0
    RK_para_c(3) =  2526269341429.0d0/ 6820363962896.0d0
    RK_para_c(4) =  2006345519317.0d0/ 3224310063776.0d0
    RK_para_c(5) =  2802321613138.0d0/ 2924317926251.0d0

  end subroutine Init_RK4S5

!==========================================================================

  subroutine RK4s5_marching_1D(Lead_Q,Dim_Q,q,q_tmp,dqdt)
    !--Declare subroutine arguments
    implicit none

    integer:: Lead_Q
    integer:: Dim_Q
    real(kind=8)::     q(0:Lead_Q)
    real(kind=8):: q_tmp(0:Lead_Q)
    real(kind=8)::  dqdt(0:Lead_Q)
    
    !--Declare local arguments
    integer:: N1,N2
    
    ! subroutine begin
    N1=Dim_Q;
    
     q_tmp(0:N1) = RK_para_a(ks)*q_tmp(0:N1) &
                          + dt*dqdt(0:N1)
    
         q(0:N1) = q(0:N1) &
                          + RK_para_b(ks)*q_tmp(0:N1) 
    
    
    return
    !
  end subroutine RK4s5_marching_1D


!==========================================================================

  subroutine RK4s5_marching_2D(Lead_Q,Dim_Q,q,q_tmp,dqdt)
    !--Declare subroutine arguments
    implicit none

    integer:: Lead_Q(1:2)
    integer:: Dim_Q(1:2)
    real(kind=8)::     q(0:Lead_Q(1),0:Lead_Q(2))
    real(kind=8):: q_tmp(0:Lead_Q(1),0:Lead_Q(2))
    real(kind=8)::  dqdt(0:Lead_Q(1),0:Lead_Q(2))
    
    !--Declare local arguments
    integer:: N1,N2
    
    ! subroutine begin
    N1=Dim_Q(1); N2=Dim_Q(2);
    
     q_tmp(0:N1,0:N2) = RK_para_a(ks)*q_tmp(0:N1,0:N2) &
                          + dt*dqdt(0:N1,0:N2)
    
         q(0:N1,0:N2) = q(0:N1,0:N2) &
                          + RK_para_b(ks)*q_tmp(0:N1,0:N2) 
    
    
    return
    !
  end subroutine RK4s5_marching_2D

!=====================================================================

  subroutine RK4s5_marching(Lead_Q,Dim_Q,q,q_tmp,dqdt)
    !--Declare subroutine arguments
    implicit none

    integer:: Lead_Q(1:3)
    integer:: Dim_Q(1:3)
    real(kind=8)::     q(0:Lead_Q(1),0:Lead_Q(2),0:Lead_Q(3))
    real(kind=8):: q_tmp(0:Lead_Q(1),0:Lead_Q(2),0:Lead_Q(3))
    real(kind=8)::  dqdt(0:Lead_Q(1),0:Lead_Q(2),0:Lead_Q(3))
    
    !--Declare local arguments
    integer:: N1,N2,N3
    
    ! subroutine begin
    N1=Dim_Q(1); N2=Dim_Q(2); N3=Dim_Q(3)
    
    q_tmp(0:N1,0:N2,0:N3) = RK_para_a(ks)*q_tmp(0:N1,0:N2,0:N3) &
                          + dt*dqdt(0:N1,0:N2,0:N3)
    
        q(0:N1,0:N2,0:N3) = q(0:N1,0:N2,0:N3) &
                          + RK_para_b(ks)*q_tmp(0:N1,0:N2,0:N3) 
    
    
    return
    !
  end subroutine RK4s5_marching
  
  
end module RK_VAR
