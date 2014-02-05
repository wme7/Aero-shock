subroutine Test001
  use universal_const    ! universal.f90
  use Legendre           ! Legendre.f90
  use MD2D_Grid          ! MD2D_Grid.f90
  use State_Var          ! State_Var.f90  
  implicit none
  integer:: i, j, k 
  integer:: ND1, ND2, ND3
  
     
  do DDK=1,TotNum_DM 
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)
!     write(*,*)ND1,ND2,ND3,pi
     do j=0,ND2
        do i=0,ND1
           
           v1(i,j,DDK) = dcos(pi*x1(i,j,DDK)) &
                       * dcos(pi*x2(i,j,DDK)) 
           v1(i,j,DDK) = dsin(pi*x1(i,j,DDK)) &
                       * dsin(pi*x2(i,j,DDK)) 
           
        enddo
     enddo
  enddo
  

  open(78,file='./Test/test001.dat')
  write(78,*) 'VARIABLES = "X", "Y", "V1", "V2"'
   
  do DDK=1,TotNum_DM
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)
     
     write(78,*) 'ZONE I=', ND1+1, ', J=', ND2+1,  'F=POINT'
     
     do j=0,ND2
        do i=0,ND1
           write(78, 10001) x1(i,j,DDK), x2(i,j,DDK), & 
                            v1(i,j,DDK), v2(i,j,DDK)
        enddo
     enddo
  enddo ! DDK

  close(78) 
10001 format(4e23.15)
  
end subroutine Test001




!!$subroutine Test100
!!$  use universal_const    ! universal.f90
!!$  use Legendre     ! Legendre.f90
!!$  use MD3D_Grid    ! MD3D_Grid.f90
!!$  use State_Var    ! State_Var.f90  
!!$  implicit none
!!$
!!$  
!!$  integer:: i,j,k
!!$  integer:: ND1, ND2, ND3
!!$  real(kind=8):: err_local, err_max
!!$  real(kind=8):: x,y,z
!!$     
!!$   
!!$  do DDK=1,1 
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$     !write(*,*)ND1,ND2,ND3,pi
!!$     do k=0,ND3
!!$        z=LGLCoord(k,ND3)
!!$        do j=0,ND2
!!$           y=LGLCoord(j,ND2)
!!$           do i=0,ND1
!!$              x=LGLCoord(i,ND1)
!!$
!!$              E1(i,j,k,DDK) = dcos(pi*x1(i,j,k,DDK)) &
!!$                            * dcos(pi*x2(i,j,k,DDK)) &
!!$                            * dcos(pi*x3(i,j,k,DDK)) 
!!$
!!$
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$     do k=0,ND3   ! Numerical derivatives
!!$        dE1_dt(0:ND1,0:ND2,k,DDK) = &
!!$           Matmul(Diff_xi1(0:ND1,0:ND1,ND1),E1(0:ND1,0:ND2,k,DDK))
!!$     enddo
!!$
!!$     do k=0,ND3   ! Exact derivatives
!!$        z=LGLCoord(k,ND3)
!!$        do j=0,ND2
!!$           y=LGLCoord(j,ND2)
!!$           do i=0,ND1
!!$              x=LGLCoord(i,ND1)
!!$
!!$              E1_tmp(i,j,k,DDK) = -pi*dsin(pi*x1(i,j,k,DDK)) &
!!$                                *     dcos(pi*x2(i,j,k,DDK)) &
!!$                                *     dcos(pi*x3(i,j,k,DDK)) 
!!$
!!$           enddo
!!$        enddo 
!!$     enddo
!!$
!!$     err_max=0.d0
!!$     err_max=maxval(abs(E1_tmp(0:ND1,0:ND2,0:ND3,DDK)-dE1_dt(0:ND1,0:ND2,0:ND3,DDK)))
!!$     write(*,*)'Error_inf:',err_max  
!!$
!!$  enddo
!!$
!!$end subroutine Test100
!!$
!!$
subroutine Test110
  use Legendre     ! Legendre.f90
  use MD2D_Grid    ! MD3D_Grid.f90
  use Metric_Var   ! Metric_Var.f90
  implicit none

  integer:: i,j,k

  do DDK=1,TotNum_DM
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)
     
     do j=0,ND2
        do i=0,ND1
           write(*,1000)i,j,DDK,Jacobin(i,j,DDK)
        enddo
     enddo
     
  enddo
1000  format(3i4,1f15.7) 
    
end subroutine Test110
!!$
!!$subroutine Test111
!!$  use universal_const ! universal.f90
!!$  use Legendre        ! Legendre.f90
!!$  use MD3D_Grid       ! MD3D_Grid.f90
!!$  use Metric_Var      ! Metric_Var.f90
!!$  use State_Var       ! State_Var.f90
!!$  implicit none
!!$
!!$  integer:: i,j,k
!!$  real(kind=8):: err_max
!!$  real(kind=8):: tmp(0:30,0:30,0:30)
!!$
!!$
!!$  do DDK=1,TotNum_DM
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$     
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              E1(i,j,k,DDK) = dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                            * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                            * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$              dE1_dt(i,j,k,DDK) = -0.5d0*pi*dsin(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                                *           dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                                *           dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$     ! differentiate wrt x1
!!$     tmp(0:ND1,0:ND2,0:ND3)=0.d0
!!$     do k=0,ND3
!!$        E1_tmp(0:ND1,0:ND2,k,DDK) = dxi1_dx1(0:ND1,0:ND2,k,DDK) * &
!!$          matmul(Diff_xi1(0:ND1,0:ND1,ND1), E1(0:ND1,0:ND2,k,DDK)) &
!!$                                + dxi2_dx1(0:ND1,0:ND2,k,DDK) * &
!!$          matmul(E1(0:ND1,0:ND2,k,DDK), Diff_xi2(0:ND2,0:ND2,ND2)) 
!!$
!!$        do i=0,ND3
!!$           tmp(0:ND1,0:ND2,i) = tmp(0:ND1,0:ND2,i) &
!!$                + Diff_xi2(k,i,ND3) * E1(0:ND1,0:ND2,k,DDK) * &
!!$                  dxi3_dx1(0:ND1,0:ND2,i,DDK)
!!$        enddo
!!$
!!$     enddo 
!!$     E1_tmp(0:ND1,0:ND2,0:ND3,DDK) = E1_tmp(0:ND1,0:ND2,0:ND3,DDK) + &
!!$           tmp(0:ND1,0:ND2,0:ND3)
!!$
!!$          
!!$
!!$     err_max = maxval( abs( E1_tmp(0:ND1,0:ND2,0:ND3,DDK) &
!!$                          - dE1_dt(0:ND1,0:ND2,0:ND3,DDK) ) )
!!$     write(*,*)'L_inf for x1 differentiation:',err_max
!!$
!!$  enddo
!!$
!!$
!!$  do DDK=1,TotNum_DM
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$     
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              E1(i,j,k,DDK) = dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$              dE1_dt(i,j,k,DDK) = -0.5d0*pi*dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dsin(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$     ! differentiate wrt x2
!!$     do k=0,ND3
!!$        E1_tmp(0:ND1,0:ND2,k,DDK) = dxi1_dx2(0:ND1,0:ND2,k,DDK) * &
!!$          matmul(Diff_xi1(0:ND1,0:ND1,ND1), E1(0:ND1,0:ND2,k,DDK)) &
!!$                                + dxi2_dx2(0:ND1,0:ND2,k,DDK) * &
!!$          matmul(E1(0:ND1,0:ND2,k,DDK), Diff_xi2(0:ND2,0:ND2,ND2)) 
!!$     enddo 
!!$
!!$     do j=0,ND2
!!$        E1_tmp(0:ND1,j,0:ND3,DDK) = E1_tmp(0:ND1,j,0:ND3,DDK) + &
!!$          dxi3_dx2(0:ND1,j,0:ND3,DDK) * &
!!$          matmul(E1(0:ND1,j,0:ND3,DDK), Diff_xi2(0:ND3,0:ND3,ND3))
!!$     enddo 
!!$          
!!$
!!$     err_max = maxval( abs( E1_tmp(0:ND1,0:ND2,0:ND3,DDK) &
!!$                          - dE1_dt(0:ND1,0:ND2,0:ND3,DDK) ) )
!!$     write(*,*)'L_inf for x2 differentiation:',err_max
!!$
!!$  enddo
!!$
!!$  do DDK=1,TotNum_DM
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$     
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              E1(i,j,k,DDK) = dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$              dE1_dt(i,j,k,DDK) = -0.5d0*pi*dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dsin(0.5d0*pi*x3(i,j,k,DDK)) 
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$     ! differentiate wrt x3
!!$     do k=0,ND3
!!$        E1_tmp(0:ND1,0:ND2,k,DDK) = dxi1_dx3(0:ND1,0:ND2,k,DDK) * &
!!$          matmul(Diff_xi1(0:ND1,0:ND1,ND1), E1(0:ND1,0:ND2,k,DDK)) &
!!$                                + dxi2_dx3(0:ND1,0:ND2,k,DDK) * &
!!$          matmul(E1(0:ND1,0:ND2,k,DDK), Diff_xi2(0:ND2,0:ND2,ND2)) 
!!$     enddo 
!!$
!!$!     E1(0:ND1,0:ND2,0:ND3,DDK,1)=0.d0
!!$     do j=0,ND2
!!$        E1_tmp(0:ND1,j,0:ND3,DDK) = E1_tmp(0:ND1,j,0:ND3,DDK) + &
!!$          dxi3_dx3(0:ND1,j,0:ND3,DDK) * &
!!$          matmul(E1(0:ND1,j,0:ND3,DDK), Diff_xi2(0:ND3,0:ND3,ND3))
!!$     enddo 
!!$          
!!$
!!$     err_max = maxval( abs( E1_tmp(0:ND1,0:ND2,0:ND3,DDK) &
!!$                          - dE1_dt(0:ND1,0:ND2,0:ND3,DDK) ) )
!!$     write(*,*)'L_inf for x3 differentiation:',err_max
!!$
!!$  enddo
!!$
!!$
!!$!1000  format(4i4,1f15.7) 
!!$
!!$end subroutine Test111
!!$
!!$
!!$subroutine Test112
!!$  use universal_const ! universal.f90
!!$  use Legendre        ! Legendre.f90
!!$  use MD3D_Grid       ! MD3D_Grid.f90
!!$!  use Metric_Var      ! Metric_Var.f90
!!$  use State_Var       ! State_Var.f90
!!$  implicit none
!!$
!!$  integer:: i,j,k
!!$  integer::ND1_max, ND2_max, ND3_max
!!$  real(kind=8):: err_max
!!$
!!$
!!$
!!$  ND1_max=PolyDegN_Max(1)
!!$  ND2_max=PolyDegN_Max(2)
!!$  ND3_max=PolyDegN_max(3)
!!$
!!$  do DDK=1,TotNum_Dm
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              E1(i,j,k,DDK) = dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$              dE1_dt(i,j,k,DDK) = -0.5d0*pi*dsin(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$           enddo
!!$        enddo
!!$     enddo
!!$     
!!$     call differentiate_xi1(PolyDegN_max(1:3),PolyDegN_DM(1:3,DDK),&
!!$                           E1(0,0,0,DDK),E1_tmp(0,0,0,DDK))
!!$
!!$     err_max=maxval(abs(E1_tmp(0:ND1,0:ND2,0:ND3,DDK)- &
!!$                        dE1_dt(0:ND1,0:ND2,0:ND3,DDK)))
!!$     write(*,*)'L_inf for xi1 differentiation:',err_max
!!$
!!$
!!$  enddo
!!$
!!$
!!$  do DDK=1,TotNum_Dm
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              E1(i,j,k,DDK) = dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$              dE1_dt(i,j,k,DDK) = -0.5d0*pi*dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dsin(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$           enddo
!!$        enddo
!!$     enddo
!!$     
!!$     call differentiate_xi2(PolyDegN_max(1:3),PolyDegN_DM(1:3,DDK),&
!!$                           E1(0,0,0,DDK),E1_tmp(0,0,0,DDK))
!!$
!!$     err_max=maxval(abs(E1_tmp(0:ND1,0:ND2,0:ND3,DDK)-&
!!$                        dE1_dt(0:ND1,0:ND2,0:ND3,DDK)))
!!$     write(*,*)'L_inf for xi2 differentiation:',err_max
!!$
!!$
!!$  enddo
!!$
!!$  do DDK=1,TotNum_Dm
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              E1(i,j,k,DDK) = dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$              dE1_dt(i,j,k,DDK) = -0.5d0*pi*dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dsin(0.5d0*pi*x3(i,j,k,DDK)) 
!!$           enddo
!!$        enddo
!!$     enddo
!!$     
!!$     call differentiate_xi3(PolyDegN_max(1:3),PolyDegN_DM(1:3,DDK),&
!!$                           E1(0,0,0,DDK),E1_tmp(0,0,0,DDK))
!!$
!!$     err_max=maxval(abs(E1_tmp(0:ND1,0:ND2,0:ND3,DDK) - &
!!$                        dE1_dt(0:ND1,0:ND2,0:ND3,DDK)))
!!$     write(*,*)'L_inf for xi3 differentiation:',err_max
!!$
!!$
!!$  enddo
!!$
!!$
!!$end subroutine Test112
!!$
!!$subroutine Test113
!!$  use universal_const ! universal.f90
!!$  use Legendre        ! Legendre.f90
!!$  use MD3D_Grid       ! MD3D_Grid.f90
!!$  use Metric_Var      ! Metric_Var.f90
!!$  use State_Var       ! State_Var.f90
!!$  implicit none
!!$
!!$  integer:: i,j,k
!!$  real(kind=8):: err_max
!!$  real(kind=8):: tmp(0:30,0:30,0:30)
!!$
!!$
!!$  do DDK=1,TotNum_DM
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$     
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              E1(i,j,k,DDK) = dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$              dE1_dt(i,j,k,DDK) = -0.5d0*pi*dsin(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$     ! differentiate wrt x1
!!$     call diff_x1(PolyDegN_Max(1:3),PolyDegN_DM(1:3,DDK),&
!!$                  E1(0,0,0,DDK), E1_tmp(0,0,0,DDK), DDK) 
!!$
!!$     err_max = maxval( abs( E1_tmp(0:ND1,0:ND2,0:ND3,DDK) &
!!$                          - dE1_dt(0:ND1,0:ND2,0:ND3,DDK) ) )
!!$     write(*,*)'L_inf for x1 differentiation:',err_max
!!$  enddo
!!$
!!$
!!$  do DDK=1,TotNum_DM
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$     
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              E1(i,j,k,DDK) = dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$              dE1_dt(i,j,k,DDK) = -0.5d0*pi*dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dsin(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$     
!!$     ! differentiate wrt x1
!!$     call diff_x2(PolyDegN_Max(1:3),PolyDegN_DM(1:3,DDK),&
!!$                  E1(0,0,0,DDK), E1_tmp(0,0,0,DDK), DDK) 
!!$     err_max = maxval( abs( E1_tmp(0:ND1,0:ND2,0:ND3,DDK) &
!!$                          - dE1_dt(0:ND1,0:ND2,0:ND3,DDK) ) )
!!$     write(*,*)'L_inf for x2 differentiation:',err_max
!!$     
!!$  enddo
!!$
!!$
!!$  do DDK=1,TotNum_DM
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$     
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              E1(i,j,k,DDK) = dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x3(i,j,k,DDK)) 
!!$              dE1_dt(i,j,k,DDK) = -0.5d0*pi*dcos(0.5d0*pi*x1(i,j,k,DDK)) &
!!$                              * dcos(0.5d0*pi*x2(i,j,k,DDK)) &
!!$                              * dsin(0.5d0*pi*x3(i,j,k,DDK)) 
!!$           enddo
!!$        enddo
!!$     enddo
!!$
!!$     ! differentiate wrt x3
!!$     call diff_x3(PolyDegN_Max(1:3),PolyDegN_DM(1:3,DDK),&
!!$                  E1(0,0,0,DDK), E1_tmp(0,0,0,DDK), DDK) 
!!$     err_max = maxval( abs( E1_tmp(0:ND1,0:ND2,0:ND3,DDK) &
!!$                          - dE1_dt(0:ND1,0:ND2,0:ND3,DDK) ) )
!!$     write(*,*)'L_inf for x3 differentiation:',err_max
!!$          
!!$
!!$     err_max = maxval( abs( E1(0:ND1,0:ND2,0:ND3,DDK,1) &
!!$                          - E1(0:ND1,0:ND2,0:ND3,DDK,2) ) )
!!$     write(*,*)'L_inf for x3 differentiation:',err_max
!!$
!!$  enddo
!!$
!!$
!!$
!!$!1000  format(4i4,1f15.7) 
!!$
!!$end subroutine Test113
!!$
!!$
!!$subroutine Test200 
!!$  use MD3D_Grid       ! MD3D_Grid.f90
!!$  use Metric_Var      ! Metric_Var.f90
!!$  use NorVec_Var      ! NorVec_Var.f90
!!$  implicit none
!!$
!!$  integer:: i,j,k
!!$  integer:: lid  
!!$
!!$  lid=81
!!$  open(lid,file='NorVec.dat')
!!$  write(lid,*)'VARIABLES = "X", "Y", "Z", "N1", "N2", "N3", "ABS N"'
!!$
!!$  do DDK=1,TotNum_DM
!!$
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$
!!$     do surf_num=1,6
!!$
!!$        select case(surf_num)
!!$        case(1) 
!!$           write(lid,*)'ZONE J=', ND2+1, ', K=', ND3+1,  'F=POINT'
!!$           do k=0,ND3
!!$              do j=0,ND2
!!$                 write(lid,1001)x1(ND1,j,k,DDK), x2(ND1,j,k,DDK), x3(ND1,j,k,DDK), &
!!$                      NorVec_x1(j,k,surf_num,DDK),NorVec_x2(j,k,surf_num,DDK),&
!!$                      NorVec_x3(j,k,surf_num,DDK),NorVec_mg(j,k,surf_num,DDK)
!!$                 
!!$              enddo
!!$           enddo
!!$           !write(lid,*)
!!$
!!$        case(2) 
!!$           write(lid,*)'ZONE J=', ND2+1, ', K=', ND3+1,  'F=POINT'
!!$           do k=0,ND3
!!$              do j=0,ND2
!!$                 write(lid,1001)x1(0,j,k,DDK),x2(0,j,k,DDK),x3(0,j,k,DDK), &
!!$                      NorVec_x1(j,k,surf_num,DDK),&
!!$                      NorVec_x2(j,k,surf_num,DDK),NorVec_x3(j,k,surf_num,DDK),&
!!$                      NorVec_mg(j,k,surf_num,DDK)
!!$                 
!!$              enddo
!!$           enddo
!!$!           pause
!!$
!!$        case(3)
!!$           write(lid,*)'ZONE I=', ND1+1, ', K=', ND3+1,  'F=POINT'
!!$           do k=0,ND3
!!$              do i=0,ND1
!!$                 write(lid,1001)x1(i,ND2,k,DDK),x2(i,ND2,k,DDK),x3(i,ND2,k,DDK), &
!!$                      NorVec_x1(i,k,surf_num,DDK),&
!!$                      NorVec_x2(i,k,surf_num,DDK),NorVec_x3(i,k,surf_num,DDK),&
!!$                      NorVec_mg(i,k,surf_num,DDK)
!!$                 
!!$              enddo
!!$           enddo
!!$
!!$        case(4)
!!$           write(lid,*)'ZONE I=', ND1+1, ', K=', ND3+1,  'F=POINT'
!!$           do k=0,ND3
!!$              do i=0,ND1
!!$                 write(lid,1001)x1(i,0,k,DDK),x2(i,0,k,DDK),x3(i,0,k,DDK), &
!!$                      NorVec_x1(i,k,surf_num,DDK),&
!!$                      NorVec_x2(i,k,surf_num,DDK),NorVec_x3(i,k,surf_num,DDK),&
!!$                      NorVec_mg(i,k,surf_num,DDK)
!!$                 
!!$              enddo
!!$           enddo
!!$!           pause
!!$
!!$        case(5)
!!$           write(lid,*)'ZONE I=', ND1+1, ', J=', ND2+1,  'F=POINT'
!!$           do j=0,ND2
!!$              do i=0,ND1
!!$                 write(lid,1001)x1(i,j,ND3,DDK),x2(i,j,ND3,DDK),x3(i,j,ND3,DDK), &
!!$                      NorVec_x1(i,j,surf_num,DDK),&
!!$                      NorVec_x2(i,j,surf_num,DDK),NorVec_x3(i,j,surf_num,DDK),&
!!$                      NorVec_mg(i,j,surf_num,DDK)
!!$                 
!!$              enddo
!!$           enddo
!!$!           pause
!!$
!!$        case(6)
!!$           write(lid,*)'ZONE I=', ND1+1, ', J=', ND2+1,  'F=POINT'
!!$           do j=0,ND2
!!$              do i=0,ND1
!!$                 write(lid,1001)x1(i,j,0,DDK),x2(i,j,0,DDK),x3(i,j,0,DDK), &
!!$                      NorVec_x1(i,j,surf_num,DDK),&
!!$                      NorVec_x2(i,j,surf_num,DDK),NorVec_x3(i,j,surf_num,DDK),&
!!$                      NorVec_mg(i,j,surf_num,DDK)
!!$                 
!!$              enddo
!!$           enddo
!!$!           pause
!!$        
!!$        end select
!!$        
!!$     enddo ! surf_num
!!$
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           do i=0,ND1
!!$              if(Jacobin(i,j,k,DDK) .lt. 0.d0) then 
!!$                 write(*,*)'jacobin < 0',i,j,k, Jacobin(i,j,k,DDK)
!!$              endif
!!$           enddo
!!$        enddo
!!$     enddo
!!$  
!!$  enddo !DDK
!!$
!!$  close(lid)
!!$
!!$
!!$1000 format(3i4,4f15.7)
!!$1001 format(7f15.7)
!!$end subroutine Test200
!!$
!!$
subroutine Test220
  use universal_const    ! universal.f90
  use Legendre     ! Legendre.f90
  use MD2D_Grid    ! MD3D_Grid.f90
  use State_Var    ! State_Var.f90  
  implicit none
  integer:: i,j,k
  integer:: ND1, ND2, ND3
  
     
  open(78,file='./Test/Initial_Field.dat')
  write(78,*) 'VARIABLES = "X", "Y", "V1", "V2", "T11", "T12", "T22"'
  
  do DDK=1,TotNum_DM
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)
     
     write(78,*) 'ZONE I=', ND1+1, ', J=', ND2+1,  'F=POINT'
     
     do j=0,ND2
        do i=0,ND1
           write(78, 10001) x1(i,j,DDK), x2(i,j,DDK), &
                            v1(i,j,DDK), v2(i,j,DDK), &
                  T11(i,j,DDK), T12(i,j,DDK), T22(i,j,DDK) 
        enddo
     enddo
  enddo ! DDK
  close(78) 
10001 format(9e23.15)

!!$  do DDK=1,TotNum_DM
!!$     ND1=PolyDegN_DM(1,DDK)
!!$     ND2=PolyDegN_DM(2,DDK)
!!$     ND3=PolyDegN_DM(3,DDK)
!!$     do k=0,ND3
!!$        do j=0,ND2
!!$           i=ND1
!!$              write(*,1001) x1(i,j,k,DDK), x2(i,j,k,DDK), x3(i,j,k,DDK), &
!!$                   E1(i,j,k,DDK), E2(i,j,k,DDK), E3(i,j,k,DDK), &
!!$                   H1(i,j,k,DDK), H2(i,j,k,DDK), H3(i,j,k,DDK)
!!$              1001 format(9f8.4) 
!!$          
!!$        enddo
!!$     enddo
!!$  enddo ! DDK
!!$    
     
end subroutine Test220
