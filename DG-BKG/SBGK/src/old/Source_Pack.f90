subroutine Source_Field(demo_case,t_coor)
  use State_Var
  use MD2D_Grid
  implicit none
  integer::demo_case
  real(kind=8):: t_coor
  
  select case(demo_case)
  case(0) 
     ! do nothing
     
  case(1) 
     call src_001(t_coor)
     
  end select
  
  
  return 
end subroutine Source_Field

subroutine src_001(t_coor)
  use State_Var
  use MD2D_Grid
  use Material_Var
  implicit none
  real(kind=8):: t_coor
  integer::i,j
  real(kind=8)::x_coor, y_coor, z_coor
  real(kind=8)::EL_demo001
  
  
  do DDK=1,TotNum_DM 
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)
     
     
     do j=0,ND2
        do i=0,ND1 
           x_coor=x1(i,j,DDK)
           y_coor=x2(i,j,DDK)
           
           fs1(i,j,DDK)=EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK), &
                x_coor,y_coor,t_coor,6)

           fs2(i,j,DDK)=EL_demo001(rho(DDK),Lame_lambda(DDK),Lame_mu(DDK), &
                x_coor,y_coor,t_coor,7)
           
        enddo
     enddo
  enddo
  
  return

end subroutine src_001
