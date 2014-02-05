subroutine Init_BC_Variables
  use MD2D_Grid
  use Char_Var
  implicit none
  integer:: N_max
  
  N_max=maxval(PolyDegN_Max(1:2));
  call alloc_mem_Char_Var(N_max,TotNum_DM)
  
  call Init_BC_Type

  call Init_Penalty_Parameters

  call Init_S_Operator  ! check this one

end subroutine Init_BC_Variables


Subroutine Init_BC_Type
  use MD2D_Grid
  use Char_Var
  implicit none
  integer:: lid, Num_DM_BC, ierr
  

  lid=81
  open(lid,file='BC.in',form='formatted', status='old')
  read(lid,*) !'==============================================='
  read(lid,*) Num_DM_BC

  if (Num_DM_BC .ne. TotNum_DM) then 
     write(*,*)'Message from BC_Pack.f90'
     write(*,*)'Inconsistent specificition of Total Number of BC'
     write(*,*)'Stop!'
     close(lid)
     stop
  endif

  read(lid,*) !'===============================================' 
  do DDK =1,TotNum_DM 
     read(lid,*) ! Domain Num
     read(lid,*) !'---------------------------------------------'
     read(lid,1000)BC_type(1,DDK),BC_type(2,DDK),BC_type(3,DDK),BC_type(4,DDK)
     read(lid,*) !'============================================='
  enddo 

  close(lid)

  return 
  
1000 format(4i4)

end subroutine Init_BC_Type

!====================================================================!

subroutine Init_Penalty_Parameters
  use MD2D_Grid
  use Legendre 
  use Char_Var
  implicit none

  do DDK=1,TotNum_DM 
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)     

     ! Edge 1    
     tau_Edge(0:ND1,1,DDK) = tau_scale/2.d0/LGLWeights( 0 ,ND2)

     ! Edge 2
     tau_Edge(0:ND2,2,DDK) = tau_scale/2.d0/LGLWeights(ND1,ND1)

     ! Edge 3    
     tau_Edge(0:ND1,3,DDK) = tau_scale/2.d0/LGLWeights(ND2,ND2)

     ! Edge 4
     tau_Edge(0:ND2,4,DDK) = tau_scale/2.d0/LGLWeights( 0 ,ND1)     
          
  enddo

  return

end subroutine Init_Penalty_Parameters

!====================================================================!

subroutine Init_S_Operator
  use MD2D_Grid
  use NorVec_Var
  use Char_Var
  implicit none

  integer:: i

  real(kind=8):: n1, n2 
  real(kind=8):: n1sq, n1n2, n2sq
  real(kind=8):: a_p, a_m, a_pm  

  do DDK=1,TotNum_DM 
     
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)     
     
     do Edge_Num=1,4 
        
        ND=ND1; if (mod(Edge_Num,2) .eq. 0) ND=ND2
        
        do i=0,ND

           n1 = NorVec_x1(i,Edge_Num,DDK)
           n2 = NorVec_x2(i,Edge_Num,DDK)
           n1sq = n1*n1
           n1n2 = n1*n2
           n2sq = n2*n2

           a_p = dsqrt(1.d0+n1*n2); a_m = dsqrt(1.d0-n1*n2)
           a_pm= a_p * a_m 
            
           S(1:5,1:5,i,Edge_Num,DDK) = 0.5d0 * reshape( & 
            (/ 0.d0 ,  0.d0, 2.d0*n2sq/a_pm, -2.d0*n1n2/a_pm, 2.d0*n1sq/a_pm, &
               1.d0 , -1.d0,       -n1/a_m ,    (n1-n2)/a_m ,        n2/a_m , &
              -1.d0 ,  1.d0,       -n1/a_m ,    (n1-n2)/a_m ,        n2/a_m , &
              -1.d0 , -1.d0,        n1/a_p ,    (n1+n2)/a_p ,        n2/a_p , &
               1.d0 ,  1.d0,        n1/a_p ,    (n1+n2)/a_p ,        n2/a_p   &
            /) , (/5,5/) )

           SI(1:5,1:5,i,Edge_Num,DDK) = transpose(S(1:5,1:5,i,Edge_Num,DDK))

           BDO(1:5,i,Edge_Num,DDK) = &
                   (/           0.d0              ,&
                                0.d0              ,&
                      NorVec_mg(i,Edge_num,DDK)   ,&
                                0.d0              ,&
                      NorVec_mg(i,Edge_num,DDK) /)  

        enddo ! i

     enddo ! Edge_Num

  enddo! DDK
           
end subroutine Init_S_Operator

!=======================================================================!

subroutine Compute_Characteristic
  use MD2D_Grid
  use Char_Var
  use State_Var
  implicit none
  integer:: DIF
  integer:: i


  do DDK=1,TotNum_DM 
     ND1=PolyDegN_DM(1,DDK) 
     ND2=PolyDegN_DM(2,DDK) 

     do Edge_Num=1,4
        ND=ND1; if(mod(Edge_Num,2) .eq. 0) ND=ND2
        
        select case(Edge_Num)
        case(1,3)
           
           if (Edge_Num .eq. 1) DIF=0
           if (Edge_Num .eq. 3) DIF=ND2
           
           ! colllect field variables
           do i=0,ND

              ! collect field variables
              RFV(1,i,Edge_Num,DDK) =  v1(i,DIF,DDK)
              RFV(2,i,Edge_Num,DDK) =  v2(i,DIF,DDK)
              RFV(3,i,Edge_Num,DDK) = T11(i,DIF,DDK)
              RFV(4,i,Edge_Num,DDK) = T12(i,DIF,DDK)
              RFV(5,i,Edge_Num,DDK) = T22(i,DIF,DDK)

              ! compute characteristic
              R(1:5,i,Edge_Num,DDK)=&
                matmul( SI(1:5,1:5,i,Edge_Num,DDK),&
                       RFV(1:5,i,Edge_Num,DDK))

           enddo
           
        case(2,4)

           if (Edge_Num .eq. 2) DIF=ND1
           if (Edge_Num .eq. 4) DIF=0
            
           do i=0,ND

              ! collect field variables
              RFV(1,i,Edge_Num,DDK) =  v1(DIF,i,DDK)
              RFV(2,i,Edge_Num,DDK) =  v2(DIF,i,DDK)
              RFV(3,i,Edge_Num,DDK) = T11(DIF,i,DDK)
              RFV(4,i,Edge_Num,DDK) = T12(DIF,i,DDK)
              RFV(5,i,Edge_Num,DDK) = T22(DIF,i,DDK)

              ! compute characteristic
              R(1:5,i,Edge_Num,DDK)=&
                matmul( SI(1:5,1:5,i,Edge_Num,DDK),&
                       RFV(1:5,i,Edge_Num,DDK))
              
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
  use Material_Var
  implicit none
  ! Declare local variables
  integer:: DI1, DI2, DIF
  integer:: i,j

  real(kind=8):: R_Diff(1:5)     !        (R-RBC)
  real(kind=8):: B_R_Diff(1:5)   !      B (R-RBC)
  real(kind=8):: PBC(1:5)        !   SI B (R-RBC)
  

  ! subroutine begins
  
  do DDK=1,TotNum_DM
     ND1=PolyDegN_DM(1,DDK) 
     ND2=PolyDegN_DM(2,DDK)

     do Edge_Num=1,4
        
        ND=ND1; if( mod(Edge_Num,2) .eq. 0) ND=ND2

        select case(Edge_Num)
        case(1,3)

           if (Edge_Num .eq. 1) DIF=0
           if (Edge_Num .eq. 3) DIF=ND2

        case(2,4)

           if (Edge_Num .eq. 2) DIF=ND1
           if (Edge_Num .eq. 4) DIF=0
           
        end select

        select case(BC_Type(Edge_Num,DDK))
        case(2) ! Non-slip rigid solid boundary condition: V1=0, V2=0
                ! R3 = R2, R5 = R4
        
           do i=0,ND
              
              RBC(1,i,Edge_Num,DDK) = R(1,i,Edge_Num,DDK)
              RBC(2,i,Edge_Num,DDK) = R(2,i,Edge_Num,DDK)
              RBC(3,i,Edge_Num,DDK) = R(2,i,Edge_Num,DDK)
              RBC(4,i,Edge_Num,DDK) = R(4,i,Edge_Num,DDK)
              RBC(5,i,Edge_Num,DDK) = R(4,i,Edge_Num,DDK)
 

              R_Diff(1:5) = R(1:5,i,Edge_Num,DDK) - RBC(1:5,i,Edge_Num,DDK)

              B_R_Diff(1:5)=matmul(S(1:5,1:5,i,Edge_Num,DDK),&
                                 BDO(1:5,i,Edge_Num,DDK)*R_Diff(1:5))

              PBC(1:5) = - (matmul(H_inv(1:5,1:5,i,Edge_Num,DDK),&
                                  B_R_Diff(1:5))) &
                           * tau_edge(i,Edge_Num,DDK)

              select case(Edge_Num)
              case(1,3)

                  dv1_dt(i,DIF,DDK) =  dv1_dt(i,DIF,DDK)+PBC(1)
                  dv2_dt(i,DIF,DDK) =  dv2_dt(i,DIF,DDK)+PBC(2)
                 dT11_dt(i,DIF,DDK) = dT11_dt(i,DIF,DDK)+PBC(3)
                 dT12_dt(i,DIF,DDK) = dT12_dt(i,DIF,DDK)+PBC(4)
                 dT22_dt(i,DIF,DDK) = dT22_dt(i,DIF,DDK)+PBC(5)

              case(2,4)

                  dv1_dt(DIF,i,DDK) =  dv1_dt(DIF,i,DDK)+PBC(1)
                  dv2_dt(DIF,i,DDK) =  dv2_dt(DIF,i,DDK)+PBC(2)
                 dT11_dt(DIF,i,DDK) = dT11_dt(DIF,i,DDK)+PBC(3)
                 dT12_dt(DIF,i,DDK) = dT12_dt(DIF,i,DDK)+PBC(4)
                 dT22_dt(DIF,i,DDK) = dT22_dt(DIF,i,DDK)+PBC(5)
              
              end select ! Edge_Num

           enddo ! i

        case(1) ! stress free boundary condition
                ! R3 = - R2, R5 = -R4

           do i=0,ND
              
              RBC(1,i,Edge_Num,DDK) =  R(1,i,Edge_Num,DDK)
              RBC(2,i,Edge_Num,DDK) =  R(2,i,Edge_Num,DDK)
              RBC(3,i,Edge_Num,DDK) = -R(2,i,Edge_Num,DDK)
              RBC(4,i,Edge_Num,DDK) =  R(4,i,Edge_Num,DDK)
              RBC(5,i,Edge_Num,DDK) = -R(4,i,Edge_Num,DDK)
 

              R_Diff(1:5) = R(1:5,i,Edge_Num,DDK) - RBC(1:5,i,Edge_Num,DDK)

              B_R_Diff(1:5)=matmul(S(1:5,1:5,i,Edge_Num,DDK),&
                                 BDO(1:5,i,Edge_Num,DDK)*R_Diff(1:5))

              PBC(1:5) = - (matmul(H_inv(1:5,1:5,i,Edge_Num,DDK),&
                                  B_R_Diff(1:5))) &
                           * tau_edge(i,Edge_Num,DDK)

              select case(Edge_Num)
              case(1,3)

                  dv1_dt(i,DIF,DDK) =  dv1_dt(i,DIF,DDK)+PBC(1)
                  dv2_dt(i,DIF,DDK) =  dv2_dt(i,DIF,DDK)+PBC(2)
                 dT11_dt(i,DIF,DDK) = dT11_dt(i,DIF,DDK)+PBC(3)
                 dT12_dt(i,DIF,DDK) = dT12_dt(i,DIF,DDK)+PBC(4)
                 dT22_dt(i,DIF,DDK) = dT22_dt(i,DIF,DDK)+PBC(5)

              case(2,4)

                  dv1_dt(DIF,i,DDK) =  dv1_dt(DIF,i,DDK)+PBC(1)
                  dv2_dt(DIF,i,DDK) =  dv2_dt(DIF,i,DDK)+PBC(2)
                 dT11_dt(DIF,i,DDK) = dT11_dt(DIF,i,DDK)+PBC(3)
                 dT12_dt(DIF,i,DDK) = dT12_dt(DIF,i,DDK)+PBC(4)
                 dT22_dt(DIF,i,DDK) = dT22_dt(DIF,i,DDK)+PBC(5)
              
              end select ! Edge_Num

           enddo ! i
           

        case(4) ! Welded interface boundary conditions
                !  (1)      (2)    (1)      (2)  
                ! R_3  = - R_2,   R_5  = - R_4

            DDK_Connect = DM_Connect(1, Edge_Num, DDK)
           Edge_Connect = DM_Connect(2, Edge_Num, DDK)
             Patch_Type = DM_Connect(3, Edge_Num, DDK)

           select case(Patch_Type)
           case(1) ! same direction
              
              do i=0,ND
              
              RBC(1,i,Edge_Num,DDK) =  R(1,i,Edge_Connect,DDK_Connect)
              RBC(2,i,Edge_Num,DDK) =  R(2,i,Edge_Connect,DDK_Connect)
              RBC(3,i,Edge_Num,DDK) = -R(2,i,Edge_Connect,DDK_Connect)
              RBC(4,i,Edge_Num,DDK) =  R(4,i,Edge_Connect,DDK_Connect)
              RBC(5,i,Edge_Num,DDK) = -R(4,i,Edge_Connect,DDK_Connect)
 

              R_Diff(1:5) = R(1:5,i,Edge_Num,DDK) - RBC(1:5,i,Edge_Num,DDK)

              B_R_Diff(1:5)=matmul(S(1:5,1:5,i,Edge_Num,DDK),&
                                 BDO(1:5,i,Edge_Num,DDK)*R_Diff(1:5))

              PBC(1:5) = - (matmul(H_inv(1:5,1:5,i,Edge_Num,DDK),&
                                  B_R_Diff(1:5))) &
                           * tau_edge(i,Edge_Num,DDK)

              select case(Edge_Num)
              case(1,3)

                  dv1_dt(i,DIF,DDK) =  dv1_dt(i,DIF,DDK)+PBC(1)
                  dv2_dt(i,DIF,DDK) =  dv2_dt(i,DIF,DDK)+PBC(2)
                 dT11_dt(i,DIF,DDK) = dT11_dt(i,DIF,DDK)+PBC(3)
                 dT12_dt(i,DIF,DDK) = dT12_dt(i,DIF,DDK)+PBC(4)
                 dT22_dt(i,DIF,DDK) = dT22_dt(i,DIF,DDK)+PBC(5)

              case(2,4)

                  dv1_dt(DIF,i,DDK) =  dv1_dt(DIF,i,DDK)+PBC(1)
                  dv2_dt(DIF,i,DDK) =  dv2_dt(DIF,i,DDK)+PBC(2)
                 dT11_dt(DIF,i,DDK) = dT11_dt(DIF,i,DDK)+PBC(3)
                 dT12_dt(DIF,i,DDK) = dT12_dt(DIF,i,DDK)+PBC(4)
                 dT22_dt(DIF,i,DDK) = dT22_dt(DIF,i,DDK)+PBC(5)
              
              end select ! Edge_Num

              enddo ! i

           case(-1) ! opposite direction
              
              do i=0,ND
              
              RBC(1,i,Edge_Num,DDK) =  R(1,ND-i,Edge_Connect,DDK_Connect)
              RBC(2,i,Edge_Num,DDK) =  R(2,ND-i,Edge_Connect,DDK_Connect)
              RBC(3,i,Edge_Num,DDK) = -R(2,ND-i,Edge_Connect,DDK_Connect)
              RBC(4,i,Edge_Num,DDK) =  R(4,ND-i,Edge_Connect,DDK_Connect)
              RBC(5,i,Edge_Num,DDK) = -R(4,ND-i,Edge_Connect,DDK_Connect)

              R_Diff(1:5) = R(1:5,i,Edge_Num,DDK) - RBC(1:5,i,Edge_Num,DDK)

              B_R_Diff(1:5) = matmul(S(1:5,1:5,i,Edge_Num,DDK),&
                                   BDO(1:5,i,Edge_Num,DDK)*R_Diff(1:5))

              PBC(1:5) = - (matmul(H_inv(1:5,1:5,i,Edge_Num,DDK),&
                                    B_R_Diff(1:5))) &
                           * tau_edge(i,Edge_Num,DDK)

              select case(Edge_Num)
              case(1,3)

                  dv1_dt(i,DIF,DDK) =  dv1_dt(i,DIF,DDK)+PBC(1)
                  dv2_dt(i,DIF,DDK) =  dv2_dt(i,DIF,DDK)+PBC(2)
                 dT11_dt(i,DIF,DDK) = dT11_dt(i,DIF,DDK)+PBC(3)
                 dT12_dt(i,DIF,DDK) = dT12_dt(i,DIF,DDK)+PBC(4)
                 dT22_dt(i,DIF,DDK) = dT22_dt(i,DIF,DDK)+PBC(5)

              case(2,4)

                  dv1_dt(DIF,i,DDK) =  dv1_dt(DIF,i,DDK)+PBC(1)
                  dv2_dt(DIF,i,DDK) =  dv2_dt(DIF,i,DDK)+PBC(2)
                 dT11_dt(DIF,i,DDK) = dT11_dt(DIF,i,DDK)+PBC(3)
                 dT12_dt(DIF,i,DDK) = dT12_dt(DIF,i,DDK)+PBC(4)
                 dT22_dt(DIF,i,DDK) = dT22_dt(DIF,i,DDK)+PBC(5)
              
              end select ! Edge_Num

              enddo ! i

           end select

        case(3)




        case(5)






        end select ! BC_Type

     enddo ! Edge_Num

  enddo ! DDK

end subroutine Impose_Penalty_BC

