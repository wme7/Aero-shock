subroutine Flux_Cal
  use MD2D_Grid
  use State_Var
  use Material_Var
  implicit none

  
  ! compute flux
  do DDK=1, TotNum_DM 
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)
     
     !---------------------------------------------------------!
     ! dv1_dt = 1/rho * (dT11_dx1 + dT12_dx2)
     
     ! compute dT11_dx1
     call diff_x1(PolyDegN_Max,PolyDegN_DM,&
          T11(0,0,DDK),dqdx(0,0),DDK)
     
     dv1_dt(0:ND1,0:ND2,DDK) = dv1_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)/rho(DDK)
     
     ! compute dT12_dx2
     call diff_x2(PolyDegN_Max,PolyDegN_DM,&
          T12(0,0,DDK),dqdx(0,0),DDK)
     
     dv1_dt(0:ND1,0:ND2,DDK) = dv1_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)/rho(DDK)

     ! multiply 
           
     !---------------------------------------------------------!
     ! dv2_dt = 1/rho * ( dT12_dx1 + dT22_dx2)
     
     ! compute dH1_dx3
     call diff_x1(PolyDegN_Max,PolyDegN_DM,&
          T12(0,0,DDK),dqdx(0,0),DDK)
     
     dv2_dt(0:ND1,0:ND2,DDK) = dv2_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)/rho(DDK)
     
     ! compute dH3_dx1
     call diff_x2(PolyDegN_Max,PolyDegN_DM,&
          T22(0,0,DDK),dqdx(0,0),DDK)
           
     dv2_dt(0:ND1,0:ND2,DDK) = dv2_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)/rho(DDK)
           
     !---------------------------------------------------------!
     ! dT11_dt = (lambda+2*mu)*dv1_dx1 + lambda *dv2_dx2

     ! dT22_dt =  lambda*dv1_dx1 + (lambda+2*mu)*dv2_dx2
     
     ! compute dv1_dx1
     call diff_x1(PolyDegN_Max,PolyDegN_DM,&
          v1(0,0,DDK),dqdx(0,0),DDK)
     
     dT11_dt(0:ND1,0:ND2,DDK) = dT11_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2) * (Lame_lambda(DDK)+2.d0*Lame_mu(DDK))

     dT22_dt(0:ND1,0:ND2,DDK) = dT22_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2) * (Lame_lambda(DDK))
     
     ! compute dv2_dx2
     call diff_x2(PolyDegN_Max,PolyDegN_DM,&
          v2(0,0,DDK),dqdx(0,0),DDK)
     
     dT11_dt(0:ND1,0:ND2,DDK) = dT11_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2) * (Lame_lambda(DDK))

     dT22_dt(0:ND1,0:ND2,DDK) = dT22_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2) * (Lame_lambda(DDK)+2.d0*Lame_mu(DDK))
     
     
     !---------------------------------------------------------!
     ! dT12_dt = mu * (dv2_dx1+dv1_dx2) 
     
     ! compute dv2_dx1
     call diff_x1(PolyDegN_Max,PolyDegN_DM,&
          v2(0,0,DDK),dqdx(0,0),DDK)
     
     dT12_dt(0:ND1,0:ND2,DDK) = dT12_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)*Lame_mu(DDK)
     
     ! compute dv1_dx2
     call diff_x2(PolyDegN_Max,PolyDegN_DM,&
          v1(0,0,DDK),dqdx(0,0),DDK)
     
     dT12_dt(0:ND1,0:ND2,DDK) = dT12_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)*Lame_mu(DDK)
     
  enddo

     
end subroutine Flux_Cal

                                 
