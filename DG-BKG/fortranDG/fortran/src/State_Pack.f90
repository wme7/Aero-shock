subroutine Flux_Cal
  use MD2D_Grid
  use State_Var
  use Material_Var
  use RK_Var
  use HOFD
  implicit none

  ! compute flux
  !$OMP PARALLEL DO DEFAULT(none) &
  SHARED(v1, v2, T11,T12, T22, dv1_dt,dv2_dt,dT11_dt, dT12_dt, dT22_dt , du1_dt, du2_dt,&
  TotNum_DM,PolyDegN_DM,Nxi_max,rho,Lame_mu, Lame_lambda,fs1,fs2) &
  PRIVATE(ND1,ND2,dqdx)
  do DDK=1, TotNum_DM 
     ND1=PolyDegN_DM(1,DDK)
     ND2=PolyDegN_DM(2,DDK)
     
     !---------------------------------------------------------!
     ! dv1_dt = 1/rho * (dT11_dx1 + dT12_dx2)
     
     ! compute dT11_dx1
     call diff_x1(Nxi_max(1:2),PolyDegN_DM(1:2,DDK),&
          T11(0,0,DDK),dqdx(0,0),DDK)                ! Diff_Pack.f90
!!     write(*,*) T11(0:ND1,0:ND2,DDK)

     dv1_dt(0:ND1,0:ND2,DDK) = dv1_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)/rho(DDK)
     
     ! compute dT12_dx2
     call diff_x2(Nxi_max(1:2),PolyDegN_DM(1:2,DDK),&
          T12(0,0,DDK),dqdx(0,0),DDK)
!!     write(*,*) T12(0:ND1,0:ND2,DDK)     
     dv1_dt(0:ND1,0:ND2,DDK) = dv1_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)/rho(DDK)
          
     dv1_dt(0:ND1,0:ND2,DDK)=dv1_dt(0:ND1,0:ND2,DDK) + &
                fs1(0:ND1,0:ND2,DDK)/rho(DDK)


!!***********
!     dv1_dt(0:ND1,0:ND2,DDK) =  dv1_dt(0:ND1,0:ND2,DDK) &
!          + 0.d0/rho(DDK)

!**********************************************   

!     dv1_dt(0:ND1,0:ND2,DDK) = &
!                  (-2.d0*lame_mu(DDK)*pi*pi/omega*dcos(pi*x1(0:ND1,0:ND2,DDK)) &
!                   *dsin(pi*x2(0:ND1,0:ND2,DDK))*dsin(omega*time_loc))/rho(DDK)   
!!     write(*,*) 'x1x2', x1(0,0,DDK), x2(0,0,DDK), time_loc
!!     write(*,*) 'sin', dsin(pi*x2(0,0,DDK)), dsin(omega*time_loc)
!!     write(*,*) 'v1', dv1_dt(0:ND1,0:ND2,DDK)
!!     pause

     
     ! multiply 
           
     !---------------------------------------------------------!
     ! dv2_dt = 1/rho * ( dT12_dx1 + dT22_dx2)
     
     ! compute dH1_dx3
     call diff_x1(Nxi_max(1:2),PolyDegN_DM(1:2,DDK),&
          T12(0,0,DDK),dqdx(0,0),DDK)
     
     dv2_dt(0:ND1,0:ND2,DDK) = dv2_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)/rho(DDK)
     
     ! compute dH3_dx1
     call diff_x2(Nxi_max(1:2),PolyDegN_DM(1:2,DDK),&
          T22(0,0,DDK),dqdx(0,0),DDK)
           
     dv2_dt(0:ND1,0:ND2,DDK) = dv2_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)/rho(DDK)
          
     dv2_dt(0:ND1,0:ND2,DDK)=dv2_dt(0:ND1,0:ND2,DDK) + &
                fs2(0:ND1,0:ND2,DDK)/rho(DDK)
          
!!***********
!     dv2_dt(0:ND1,0:ND2,DDK) =  dv2_dt(0:ND1,0:ND2,DDK) &
!          + (2.d0*lame_mu(DDK)*pi*pi/omega*dsin(pi*x1(0:ND1,0:ND2,DDK)) &
!                   *dcos(pi*x2(0:ND1,0:ND2,DDK))*dsin(omega*time_loc))/rho(DDK)
!

!     dv2_dt(0:ND1,0:ND2,DDK) = &
!                  (2.d0*lame_mu(DDK)*pi*pi/omega*dsin(pi*x1(0:ND1,0:ND2,DDK)) &
!                   *dcos(pi*x2(0:ND1,0:ND2,DDK))*dsin(omega*time_loc))/rho(DDK)
!                    
!!     write(*,*) 'x1x2', x1(0,0,DDK), x2(0,0,DDK), time_loc
!!     write(*,*) 'sin', dsin(pi*x2(0,0,DDK)), dsin(omega*time_loc)
!!     write(*,*) 'v2', dv2_dt(0:ND1,0:ND2,DDK)
!!     pause
           
     !---------------------------------------------------------!
     ! dT11_dt = (lambda+2*mu)*dv1_dx1 + lambda *dv2_dx2

     ! dT22_dt =  lambda*dv1_dx1 + (lambda+2*mu)*dv2_dx2
     
     ! compute dv1_dx1
     call diff_x1(Nxi_max(1:2),PolyDegN_DM(1:2,DDK),&
          v1(0,0,DDK),dqdx(0,0),DDK)

!     write(*,*) 'dv1dx1_error', &
!                maxval(dqdx(0:ND1,0:ND2) &
!                  +pi*dsin(pi*x1(0:ND1,0:ND2,DDK))*dsin(pi*x2(0:ND1,0:ND2,DDK)))
     
     dT11_dt(0:ND1,0:ND2,DDK) = dT11_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2) * (Lame_lambda(DDK)+2.d0*Lame_mu(DDK))

     dT22_dt(0:ND1,0:ND2,DDK) = dT22_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2) * (Lame_lambda(DDK))
     
     ! compute dv2_dx2
     call diff_x2(Nxi_max(1:2),PolyDegN_DM(1:2,DDK),&
          v2(0,0,DDK),dqdx(0,0),DDK)

!     write(*,*) 'dv2dx2_error', &
!                maxval(dqdx(0:ND1,0:ND2) &
!                  -pi*dsin(pi*x1(0:ND1,0:ND2,DDK))*dsin(pi*x2(0:ND1,0:ND2,DDK)))
     
     dT11_dt(0:ND1,0:ND2,DDK) = dT11_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2) * (Lame_lambda(DDK))

     dT22_dt(0:ND1,0:ND2,DDK) = dT22_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2) * (Lame_lambda(DDK)+2.d0*Lame_mu(DDK))
!!***********
!     dT11_dt(0:ND1,0:ND2,DDK) = dT11_dt(0:ND1,0:ND2,DDK) &
!          +lame_lambda(DDK) &
!          *pi*dsin(pi*x1(0:ND1,0:ND2,DDK))*dsin(pi*x2(0:ND1,0:ND2,DDK)) &
!          *dcos(omega*time_loc)
!     
!     
!     dT22_dt(0:ND1,0:ND2,DDK) = dT22_dt(0:ND1,0:ND2,DDK) &
!          +(lame_lambda(DDK)+2.d0*lame_mu(DDK)) &
!          *pi*dsin(pi*x1(0:ND1,0:ND2,DDK))*dsin(pi*x2(0:ND1,0:ND2,DDK)) &
!          *dcos(omega*time_loc)

     

!     dT11_dt(0:ND1,0:ND2,DDK) = &
!                  (lame_lambda(DDK)+2.d0*lame_mu(DDK)) &
!                  *(-pi*dsin(pi*x1(0:ND1,0:ND2,DDK))*dsin(pi*x2(0:ND1,0:ND2,DDK)) &
!                  *dcos(omega*time_loc)) &
!                  +lame_lambda(DDK) &
!                    *pi*dsin(pi*x1(0:ND1,0:ND2,DDK))*dsin(pi*x2(0:ND1,0:ND2,DDK)) &
!                  *dcos(omega*time_loc)


!    dT22_dt(0:ND1,0:ND2,DDK) = &
!                  lame_lambda(DDK)*(-pi*dsin(pi*x1(0:ND1,0:ND2,DDK))*dsin(pi*x2(0:ND1,0:ND2,DDK)) &
!                  *dcos(omega*time_loc)) &
!                  +(lame_lambda(DDK)+2.d0*lame_mu(DDK)) &
!                  *pi*dsin(pi*x1(0:ND1,0:ND2,DDK))*dsin(pi*x2(0:ND1,0:ND2,DDK)) &
!                  *dcos(omega*time_loc)
!!     write(*,*) 'x1x2', x1(0,0,DDK), x2(0,0,DDK), time_loc
!!     write(*,*) 'sin', dsin(pi*x2(0,0,DDK)), dsin(omega*time_loc)
!!     write(*,*) 'T22', dT22_dt(0:ND1,0:ND2,DDK)
!!     pause

    
     !---------------------------------------------------------!
     ! dT12_dt = mu * (dv2_dx1+dv1_dx2) 
     
     ! compute dv2_dx1
     call diff_x1(Nxi_max(1:2),PolyDegN_DM(1:2,DDK),&
          v2(0,0,DDK),dqdx(0,0),DDK)
       
!     write(*,*) 'dv2dx1_error', &
!               maxval(dqdx(0:ND1,0:ND2) &
!               +pi*dcos(pi*x1(0:ND1,0:ND2,DDK))*dcos(pi*x2(0:ND1,0:ND2,DDK)))
     
     dT12_dt(0:ND1,0:ND2,DDK) = dT12_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)*Lame_mu(DDK)
     
     ! compute dv1_dx2
     call diff_x2(Nxi_max(1:2),PolyDegN_DM(1:2,DDK),&
          v1(0,0,DDK),dqdx(0,0),DDK)
          
!     write(*,*) 'dv1dx2_error', &
!                maxval(dqdx(0:ND1,0:ND2) &
!                  -pi*dcos(pi*x1(0:ND1,0:ND2,DDK))*dcos(pi*x2(0:ND1,0:ND2,DDK)))
                  
     dT12_dt(0:ND1,0:ND2,DDK) = dT12_dt(0:ND1,0:ND2,DDK) &
          + dqdx(0:ND1,0:ND2)*Lame_mu(DDK)
!!***********
!     dT12_dt(0:ND1,0:ND2,DDK) = dT12_dt(0:ND1,0:ND2,DDK) &
!          +lame_mu(DDK)*(pi*dcos(pi*x1(0:ND1,0:ND2,DDK))*dcos(pi*x2(0:ND1,0:ND2,DDK)) &
!          *dcos(omega*time_loc))
!

     
!    dT12_dt(0:ND1,0:ND2,DDK) = &
!                  lame_mu(DDK)*(pi*dcos(pi*x1(0:ND1,0:ND2,DDK))*dsin(pi*x2(0:ND1,0:ND2,DDK)) &
!                  *dcos(omega*time_loc) &
!                  -pi*dcos(pi*x1(0:ND1,0:ND2,DDK))*dcos(pi*x2(0:ND1,0:ND2,DDK)) &
!                  *dcos(omega*time_loc))
         
!!     write(*,*) 'x1x2', x1(0,0,DDK), x2(0,0,DDK), time_loc
!!     write(*,*) 'sin', dsin(pi*x2(0,0,DDK)), dsin(omega*time_loc)
!!     write(*,*) 'T22', dT22_dt(0:ND1,0:ND2,DDK)
!!     pause

        du1_dt = v1
        
        du2_dt = v2

     
  enddo
  !$OMP END PARALLEL DO  
end subroutine Flux_Cal

                                 
