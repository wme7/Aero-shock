subroutine Init_Normal_Vector
  use Metric_Var
  use NorVec_Var
  use MD2D_Grid
  use Legendre
  implicit none

  integer:: i,j

  ! subroutine begin
  ! allocate memory for Normal Vector Veriables

  call alloc_mem_norvec

  ! Compute Normal Vectors of all the surface of a cube

  do DDK = 1,TotNum_DM

     ND1=PND1
     ND2=PND2
     
     do Edge_Num=1,4
        
        select case(Edge_Num)
        case(1) ! n = - grad xi_2
              
           NorVec_mg(0:ND1,Edge_Num,DDK) = dsqrt( &
                dxi2_dx1(0:ND1,0,DDK)**2 + dxi2_dx2(0:ND1,0,DDK)**2 )
           
           NorVec_x1(0:ND1,Edge_Num,DDK) =  - &
                ( dxi2_dx1(0:ND1,0,DDK) / NorVec_mg(0:ND1,Edge_Num,DDK) )
           
           NorVec_x2(0:ND1,Edge_Num,DDK) =  - &
                ( dxi2_dx2(0:ND1,0,DDK) / NorVec_mg(0:ND1,Edge_Num,DDK) )
                   
        case(2) ! n = + grad xi_1

           NorVec_mg(0:ND2,Edge_Num,DDK) = dsqrt( &
                dxi1_dx1(ND1,0:ND2,DDK)**2 + dxi1_dx2(ND1,0:ND2,DDK)**2 )
           
           NorVec_x1(0:ND2,Edge_Num,DDK) =   &
                ( dxi1_dx1(ND1,0:ND2,DDK) / NorVec_mg(0:ND2,Edge_Num,DDK) )
           
           NorVec_x2(0:ND2,Edge_Num,DDK) =   &
                ( dxi1_dx2(ND1,0:ND2,DDK) / NorVec_mg(0:ND2,Edge_Num,DDK) )
           
        case(3) ! n = + grad xi_2

           NorVec_mg(0:ND1,Edge_Num,DDK) = dsqrt( &
                dxi2_dx1(0:ND1,ND2,DDK)**2 + dxi2_dx2(0:ND1,ND2,DDK)**2 )
           
           NorVec_x1(0:ND1,Edge_Num,DDK) =    &
                ( dxi2_dx1(0:ND1,ND2,DDK) / NorVec_mg(0:ND1,Edge_Num,DDK) )
           
           NorVec_x2(0:ND1,Edge_Num,DDK) =    &
                ( dxi2_dx2(0:ND1,ND2,DDK) / NorVec_mg(0:ND1,Edge_Num,DDK) )

        case(4) ! n = -grad xi_1

           NorVec_mg(0:ND2,Edge_Num,DDK) = dsqrt( &
                dxi1_dx1(0,0:ND2,DDK)**2 + dxi1_dx2(0,0:ND2,DDK)**2 )
           
           NorVec_x1(0:ND2,Edge_Num,DDK) = - &
                ( dxi1_dx1(0,0:ND2,DDK) / NorVec_mg(0:ND2,Edge_Num,DDK) )
           
           NorVec_x2(0:ND2,Edge_Num,DDK) = - &
                ( dxi1_dx2(0,0:ND2,DDK) / NorVec_mg(0:ND2,Edge_Num,DDK) )

        end select

!!$        select case(Edge_Num)
!!$        case(1,3)
!!$
!!$           do i=0,ND1
!!$              write(*,10000)i,DDK,NorVec_x1(i,Edge_Num,DDK),NorVec_x2(i,Edge_Num,DDK),&
!!$                        NorVec_mg(i,Edge_Num,DDK)
!!$           enddo
!!$
!!$        case(2,4)
!!$
!!$           do i=0,ND2
!!$              write(*,10000)i,DDK,NorVec_x1(i,Edge_Num,DDK),NorVec_x2(i,Edge_Num,DDK),&
!!$                        NorVec_mg(i,Edge_Num,DDK)
!!$           enddo
!!$
!!$        end select
        
     enddo

  enddo

10000 format(2i4,3f15.7)
             
end subroutine Init_Normal_Vector

!============================================================================
subroutine DG2D_Initial(Time_final)
  use Metric_Var
  use NorVec_Var
  use MD2D_Grid
  use Legendre
  use universal_const
  use State_Var

  implicit none

  integer:: i,j,k,ierr,ii,jj
!  real(kind=8) :: a_vec(1:2),
  real(kind=8) ::Time_final
  real(kind=8),allocatable :: b_ini(:,:),Pi_Pm(:,:),Pi_Pj(:,:)
  real(kind=8),allocatable :: temp(:,:),tempsum(:),func(:,:,:)
  real(kind=8),allocatable :: LGLWeights_2D(:,:)
  real(kind=8),allocatable :: ValOfPolyNatGrids(:) ! ,cs_1d(:),w_1d(:)

  allocate( b_ini(0:PDeg1,0:PDeg2),Pi_Pm(0:PDeg1,0:PDeg2), &
            Pi_Pj(0:PND1,0:PND2),func(0:PND1,0:PND2,1:TotNum_DM), &
            temp(0:PND1,0:PND2),tempsum(0:PND1), &
            LGLWeights_2D(0:PND1,0:PND1) ,stat=ierr )
  
  LGLWeights_2D = LGLWeights_Grid_xi1*LGLWeights_Grid_xi2

  if (ierr .ne. 0) then
     write(*,*)"Can't allocate u"
     stop
  endif

if (I_prob .ge. 1) then
allocate( ValOfPolyNatGrids(1:IGH))
write(6,*) IGH

GH=0d0

call ZEHEGA(IGH,GH,ValOfPolyNatGrids)
call WEHEGA(IGH,GH,GHW)

write(6,*) GH
write(6,*) GHW

Vx=GH
Vy=GH

!GH(1)=a_vec(1)
!GH(2)=a_vec(2)
!GH(2:IGH)=a_vec(2)

write(6,*) GH
deallocate( ValOfPolyNatGrids)

endif
  
  !setup initial condition
if (I_prob .ge. 2) then
  do DDK = 1,TotNum_DM
     do i = 0,PND1
        do j = 0,PND2
           func(i,j,DDK) = dcos(2*pi*(x1(i,j,DDK)+x2(i,j,DDK)))
        enddo
     enddo
  enddo

else
  do DDK = 1,TotNum_DM
     do i = 0,PND1
        do j = 0,PND2
           func(i,j,DDK) = dsin(2*pi*(x1(i,j,DDK)+x2(i,j,DDK)))
        enddo
     enddo
  enddo
endif



  do i = 0,PDeg1
     do j = 0,PDeg2
        Pi_Pm(i,j) = sum(Leg_Grid_xi1(:,0,i)*Leg_Grid_xi1(:,0,j)*LGLWeights_Grid_xi1(:,0))
     enddo
  enddo

  do DDK = 1,TotNum_DM
!     B(1,DDK) = a_vec(1)*dx2_dxi2(0,0,DDK)-a_vec(2)*dx1_dxi2(0,0,DDK)
!     B(2,DDK) = -a_vec(1)*dx2_dxi1(0,0,DDK)+a_vec(2)*dx1_dxi1(0,0,DDK) 
     
     Bq(2,2,DDK) = dx2_dxi2(0,0,DDK)
     Bq(2,1,DDK) = dx2_dxi1(0,0,DDK)
     Bq(1,2,DDK) = dx1_dxi2(0,0,DDK)
     Bq(1,1,DDK) = dx1_dxi1(0,0,DDK)
!     B(1,DDK) = a_vec(1)*Bq(2,2,DDK)-a_vec(2)*Bq(1,2,DDK)
!     B(2,DDK) = -a_vec(1)*Bq(2,1,DDK)+a_vec(2)*Bq(1,1,DDK)

  do i = 1,IGH
     do j = 1,IGH
     B_gen(1,i,j) = GH(i)*Bq(2,2,DDK)-GH(j)*Bq(1,2,DDK)
     B_gen(2,i,j) = -GH(i)*Bq(2,1,DDK)+GH(j)*Bq(1,1,DDK)
     enddo
  enddo   


do ii=1,IGH
do jj=1,IGH
     b_ini = 0.0d0
     do i = 0,PDeg1
        do j = 0,PDeg2
           Pi_Pj = Leg_Grid_xi1(:,:,i)*Leg_Grid_xi2(:,:,j)
           temp(:,:) = func(:,:,DDK)*Pi_Pj*Jacobian(:,:,DDK)*LGLWeights_2D
           do k = 0,PND1
              tempsum(k) = sum(temp(k,:))
           enddo
           b_ini(i,j) = sum(tempsum)
           A1(i,j,DDK) = Pi_Pm(i,i)*Pi_Pm(j,j)*Jacobian(0,0,DDK)
           A1(i,j,DDK) = dble(1/A1(i,j,DDK))
        enddo   
     enddo
     F_new(ii,jj,:,:,DDK) = b_ini*A1(:,:,DDK)
  enddo
enddo
enddo

  do i = 0,PDeg1
     do j = 0,PDeg2
        tempsum = Leg_Grid_xi1(:,0,i)*DLeg_Grid_xi1(:,0,j)*LGLWeights_Grid_xi1(:,0)
        B_tal1x(i,j) = sum(tempsum)

        tempsum = Leg_Grid_xi2(0,:,j)*Leg_Grid_xi2(0,:,j)*LGLWeights_Grid_xi2(0,:)
        B_tal1y(i,j) = sum(tempsum)

        tempsum = Leg_Grid_xi1(:,0,i)*Leg_Grid_xi1(:,0,i)*LGLWeights_Grid_xi1(:,0)
        B_tal2x(i,j) = sum(tempsum)

     enddo
  enddo
  write(6,*) F_new(1,1,0:PDeg1,0:PDeg2,1)
  write(6,*) F_new(1,2,0:PDeg1,0:PDeg2,1)

  deallocate( b_ini,Pi_Pj,Pi_Pm,temp,tempsum,func,LGLWeights_2D )

write(6,*) "End of DG2D Initial"
end subroutine DG2D_Initial
!============================================================================
subroutine DG2D_Edge()!(a_vec)
  use Metric_Var
  use NorVec_Var
  use MD2D_Grid
  use Legendre
use State_Var

  use universal_const
  implicit none

  integer:: i,j,ierr,edg_cho
!  real(kind=8) :: a_vec(1:2),
  real(kind=8) :: Nor1_mns,Nor2_mns,Normg,a_n
  real(kind=8),allocatable :: LGLWeights_xi1(:),Pj_t(:,:)

  call alloc_mem_DG2D_Edge

  allocate( LGLWeights_xi1(0:PND1),Pj_t(0:PND1,0:PDeg1), stat=ierr )
  
  if (ierr .ne. 0) then
     write(*,*)"Can't allocate DG_edge"
     stop
  endif

  LGLWeights_xi1(:) = LGLWeights_Grid_xi1(0:PND1,0)
  do i=0,PDeg1
     Pj_t(0:PND1,i) = Leg_Grid_xi1(0:PND1,0,i)
  enddo

!  do DDK = 1,TotNum_DM
DDK=1
do i=1,IGH
do j=1,IGH
    !------------ edge 1 ----------------

    edg_cho = 1
    Nor1_mns = NorVec_x1(0,edg_cho,DDK)
    Nor2_mns = NorVec_x2(0,edg_cho,DDK)
    Normg = NorVec_mg(0,edg_cho,DDK)
!    a_n = a_vec(1)*Nor1_mns + a_vec(2)*Nor2_mns
    a_n = GH(i)*Nor1_mns + GH(j)*Nor2_mns  

!    E(1,edg_cho,DDK) 
E_gen(1,edg_cho,i,j)= 0.5*( a_n - abs(a_n) )*Jacobian(0,0,DDK)*Normg
!    E(2,edg_cho,DDK) 
E_gen(2,edg_cho,i,j)= 0.5*( a_n + abs(a_n) )*Jacobian(0,0,DDK)*Normg

    !------------ edge 2 ----------------

    edg_cho = 2
    Nor1_mns = NorVec_x1(0,edg_cho,DDK)
    Nor2_mns = NorVec_x2(0,edg_cho,DDK)
    Normg = NorVec_mg(0,edg_cho,DDK)
!    a_n = a_vec(1)*Nor1_mns + a_vec(2)*Nor2_mns
    a_n = GH(i)*Nor1_mns + GH(j)*Nor2_mns

!    E(1,edg_cho,DDK) 
E_gen(1,edg_cho,i,j) = 0.5*( a_n - abs( a_n ) )*Jacobian(0,0,DDK)*Normg
!    E(2,edg_cho,DDK)
E_gen(2,edg_cho,i,j) = 0.5*( a_n + abs( a_n ) )*Jacobian(0,0,DDK)*Normg

    !------------ edge 3 ----------------

    edg_cho = 3
    Nor1_mns = NorVec_x1(0,edg_cho,DDK)
    Nor2_mns = NorVec_x2(0,edg_cho,DDK)
    Normg = NorVec_mg(0,edg_cho,DDK)
!    a_n = a_vec(1)*Nor1_mns + a_vec(2)*Nor2_mns
    a_n = GH(i)*Nor1_mns + GH(j)*Nor2_mns

!    E(1,edg_cho,DDK)
E_gen(1,edg_cho,i,j) = 0.5*( a_n - abs( a_n ) )*Jacobian(0,0,DDK)*Normg
!    E(2,edg_cho,DDK)
E_gen(2,edg_cho,i,j) = 0.5*( a_n + abs( a_n ) )*Jacobian(0,0,DDK)*Normg

    !------------ edge 4 ----------------

    edg_cho = 4
    Nor1_mns = NorVec_x1(0,edg_cho,DDK)
    Nor2_mns = NorVec_x2(0,edg_cho,DDK)
    Normg = NorVec_mg(0,edg_cho,DDK)
!    a_n = a_vec(1)*Nor1_mns + a_vec(2)*Nor2_mns
    a_n = GH(i)*Nor1_mns + GH(j)*Nor2_mns

!    E(1,edg_cho,DDK)
E_gen(1,edg_cho,i,j) = 0.5*( a_n - abs( a_n ) )*Jacobian(0,0,DDK)*Normg
!    E(2,edg_cho,DDK)
E_gen(2,edg_cho,i,j) = 0.5*( a_n + abs( a_n ) )*Jacobian(0,0,DDK)*Normg
  
enddo
enddo
!enddo
    
do i = 0,PDeg1
   do j = 0,PDeg2
      
      C1x_Edge1(i,j) = sum( Pj_t(:,j)*Pj_t(:,j)*LGLWeights_xi1 )

      C2x_Edge1(i,j) = sum( Pj_t(:,j)*Pj_t(:,j)*LGLWeights_xi1 )

      C2x_Edge1(i,j) = ((-1)**(i))*C2x_Edge1(i,j)

      C1x_Edge3(i,j) = sum( Pj_t(:,i)*Pj_t(:,i)*LGLWeights_xi1 )

      C1x_Edge3(i,j) = ((-1)**(j))*C1x_Edge3(i,j)

      C2x_Edge3(i,j) = sum( Pj_t(:,i)*Pj_t(:,i)*LGLWeights_xi1 )

   enddo

   Cy_Edge1(1,i) = (-1)**(i)
  
   Cy_Edge2(i,1) = 1

   Cy_Edge3(1,i) = 1

   Cy_Edge4(i,1) = (-1)**(i)        

enddo

  deallocate( LGLWeights_xi1,Pj_t)
  
end subroutine DG2D_Edge

!============================================================================
