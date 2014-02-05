
subroutine Init_Penalty_Parameters!(ttau_scale)
  use MD2D_Grid
  use Legendre 
  use Char_Var
  implicit none
!  real(kind=8) :: ttau_scale

  do DDK=1,TotNum_DM 
     ND1=PND1
     ND2=PND2     

     ! Edge 1    
     tau_Edge(0:ND1,1,DDK) = ttau_scale/2.d0/LGLWeights( 0 ,ND2)

     ! Edge 2
     tau_Edge(0:ND2,2,DDK) = ttau_scale/2.d0/LGLWeights(ND1,ND1)

     ! Edge 3    
     tau_Edge(0:ND1,3,DDK) = ttau_scale/2.d0/LGLWeights(ND2,ND2)

     ! Edge 4
     tau_Edge(0:ND2,4,DDK) = ttau_scale/2.d0/LGLWeights( 0 ,ND1)     
          
  enddo
  

  return

end subroutine Init_Penalty_Parameters

!====================================================================!

!=======================================================================!

subroutine Compute_Characteristic
  use MD2D_Grid
  use Char_Var
  use State_Var
  implicit none
  integer:: DIF
  integer:: i


  do DDK=1,TotNum_DM 
     ND1=PND1 
     ND2=PND2 

     do Edge_Num=1,4
        ND=ND1; if(mod(Edge_Num,2) .eq. 0) ND=ND2
        
        select case(Edge_Num)
        case(1,3)
           
           if (Edge_Num .eq. 1) DIF=0
           if (Edge_Num .eq. 3) DIF=ND2
           
           ! colllect field variables
           do i=0,ND
              ! collect field variables
              R(1,i,Edge_Num,DDK) =  u(i,DIF,DDK)

           enddo
           
        case(2,4)

           if (Edge_Num .eq. 2) DIF=ND1
           if (Edge_Num .eq. 4) DIF=0
            
           do i=0,ND
              ! collect field variables
              R(1,i,Edge_Num,DDK) =  u(DIF,i,DDK)
              
           enddo
        
        end select
                
     enddo

  enddo

end subroutine Compute_Characteristic

!=======================================================================

subroutine Impose_Penalty_BC
  use MD2D_Grid
  use State_Var
  use Char_Var
  use NorVec_Var
  implicit none
  ! Declare local variables
  integer:: DI1, DI2, DIF
  integer:: i,j

  real(kind=8):: PBC        !   SI B (R-RBC)
!  ! subroutine begins
!  
  do DDK=1,TotNum_DM
     ND1=PND1 
     ND2=PND2

     do Edge_Num=1,4,3

        DDK_Connect = DM_Connect(1, Edge_Num, DDK)
       Edge_Connect = DM_Connect(2, Edge_Num, DDK)
         Patch_Type = DM_Connect(3, Edge_Num, DDK)
         

        ND=ND1; if( mod(Edge_Num,2) .eq. 0) ND=ND2

        select case(Edge_Num)
        case(1,3)

           if (Edge_Num .eq. 1) DIF=0
           if (Edge_Num .eq. 3) DIF=ND2

        case(2,4)

           if (Edge_Num .eq. 2) DIF=ND1
           if (Edge_Num .eq. 4) DIF=0
           
        end select

           select case(Patch_Type)
           case(1) ! same direction
              
              do i=0,ND
              
              PBC = - (R(1,i,Edge_Num,DDK) - R(1,i,Edge_Connect,DDK_Connect)) *&
                      tau_edge(i,Edge_Num,DDK)*NorVec_mg(i,Edge_Num,DDK)
                      
!              write(*,*) DDK,Edge_Num,"PBC=", PBC

              select case(Edge_Num)
              case(1,3)

                  du_dt(i,DIF,DDK) =  du_dt(i,DIF,DDK)+PBC

              case(2,4)

                  du_dt(DIF,i,DDK) =  du_dt(DIF,i,DDK)+PBC
              
              end select ! Edge_Num

              enddo ! i

           case(-1) ! opposite direction
              
              do i=0,ND
              
              PBC = - (R(1,i,Edge_Num,DDK) - R(1,ND-i,Edge_Connect,DDK_Connect)) &
                           * tau_edge(i,Edge_Num,DDK)

!              write(*,*) DDK,Edge_Num,"PBC=", PBC
              select case(Edge_Num)
              case(1,3)

                  du_dt(i,DIF,DDK) =  du_dt(i,DIF,DDK)+PBC

              case(2,4)

                  du_dt(DIF,i,DDK) =  du_dt(DIF,i,DDK)+PBC
              
              end select ! Edge_Num

              enddo ! i

           end select !Patch_type BC case4

     enddo ! Edge_Num

  enddo ! DDK

end subroutine Impose_Penalty_BC

