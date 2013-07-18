subroutine Init_Normal_Vector
  use Metric_Var
  use NorVec_Var
  use MD2D_Grid
  implicit none

  integer i,j
  ! subroutine begin
  ! allocate memory for Normal Vector Veriables

  call alloc_mem_norvec(maxval(PolyDegN_Max(1:2)),TotNum_DM)


  ! Compute Normal Vectors of all the surface of a cube

  do DDK=1, TotNum_DM

     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)
     
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
