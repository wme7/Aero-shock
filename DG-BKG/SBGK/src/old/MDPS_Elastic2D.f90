Program MDPS_Elastic2D
use universal_const    ! universal.f90
use MD2D_Grid    ! MD3D_Grid.f90
use Legendre     ! Legendre.f90
use State_Var    ! State_Var.f90
use Metric_Var   ! Metric_Var.f90
use RK_Var       ! RK_Var.f90
use Material_Var     
!!$use Filter_Var
!!$implicit none
!!$
!!$integer i,j,k 
!!$real(kind=8):: x_coor, y_coor, z_coor
!!$real(kind=8):: err_loc(1:6), err_max(1:6)
real(kind=8):: EL_demo001
!!$real(kind=8):: error
!!$integer:: err_count
!!$
!!$
!!$
! Program Begin 

! Initialize stage
! ---- Initialize universal constants
  call init_universal_const
! 
! ---- Initialize Domain Parameter 
  call Input_Domain_Grid_Parameters()        !-- MD3D_Grid.f90
  write(*,*)'Complete Initializing Grid Parameters'
  call Geometric_Parameters_On_Screen( )     !-- MD3D_Grid.f90

! ---- Initialize Legendre Pseudospectral Module
!      1. Gauss-Lobatto-Legendre Grid
!      2. Gauss-Lobatto-Legendre Quadurature Weights
!      3. Legenedre Differentiation Matrix
!
  call Init_LGL_Diff_Matrix_23D(PolyDegN_max(1), &
                                PolyDegN_max(2), &
                                0  )  !-- Legendre.f90
  write(*,*)'Complete Initializing Legendre Pseudospectral Operators'

! ---- Initialize computational nodes in physical space
!      by transfinite blending function

  call Init_Physical_Grid_2D !-- Grid2D_Pack.f90
!
  call Init_Patch_Direction !-- MD3D_Grid.f90
   write(*,*)'Complete Initializing Patching Direction'

! ---- Initialize EM wave field Variables

  call Init_State_Variables(PolyDegN_Max(1:2),TotNum_DM)  !-- State_Var.f90
  write(*,*)'Complete Initializing State Variables'

  call Test001 !-- Test_Pack.f90
   
!!$!  call Test100
!!$
! ---- Initialize Material Parameters
  call Init_Material_Parameters !-- Material_Pack.f90
  call Init_Hinv_Symmetrizer
!!$
!!$
!!$
! ---- Initialize Metric
  call Init_Metric !-- Metric_Pack.f90
  write(*,*)'Complete Initializing Metric Variables'

!  call Test110

!!$!  call Test111
!!$
!!$!  call Test112
!!$
!!$!  call Test113
!!$
! --- Initialize Normal Vector

  call Init_Normal_Vector !-- BC_Pack.f90
  write(*,*)'Complete Initializing Normal Vector'
!!$!  call Test200
!!$
! --- Initialize BC Parameters
  call Init_BC_Variables !-- BC_Pack.f90
  write(*,*)'Complete Initializing Variables for imposing Boundary Condition'
!!$
!!$
!!$! --- Initialize filter parameters
!!$!  call Init_Filter_Oper(maxval(PolyDegN_Max(1:3)))
!!$
!!$
! --- Initialize Time Integration
  call Input_RK_parameter ! RK_Var.f90
  call Init_RK4S5()       ! RK_Var.f90
  write(*,*)'Complete Initializing RK Parameters'
!!$  
!!$
  call Compute_time_step !-- Time_Pack.f90
  write(*,*)'Complete Initializing time step' 
   
!!$  dt= CFL/dsqrt(2.d0)/dble(PolyDegN_DM(1,1))**2
  write(*,*)dt,CFL/dsqrt(2.d0)/dble(PolyDegN_DM(1,1))**2 

  time = 0.d0
! Initialize Initialize field
  call Initial_Field(1)  !IC_Pack.f90   
  write(*,*)'Complete Initializing Initial Condition'
  call Test220
  write(*,*)'here 220'


! RK adavncing
  First_Comp=1
  Last_Comp=0
  err_count=0
  do while (time .lt. Time_final) 
!     write(*,*)'iter'
     if ( abs(time-0.d0) .ge. 1E-12) then 
        First_Comp=0
     endif

     if ( (time+dt) .ge. Time_final) then 
        dt = Time_final-time
        Last_Comp=1
     endif 

     ! loop over ks
     do ks=1,5
        
        time_loc=time+dt*RK_para_c(ks)
 !       write(*,*)'time_local',ks,time_local
        !Initialize flux
        !------------------------------------------------------
        do DDK=1, TotNum_DM 
           ND1=PolyDegN_DM(1,DDK)
           ND2=PolyDegN_DM(2,DDK)


            dv1_dt(0:ND1,0:ND2,DDK) = 0.d0 
            dv2_dt(0:ND1,0:ND2,DDK) = 0.d0
        
           dT11_dt(0:ND1,0:ND2,DDK) = 0.d0 
           dT12_dt(0:ND1,0:ND2,DDK) = 0.d0
           dT22_dt(0:ND1,0:ND2,DDK) = 0.d0

        enddo
        !------------------------------------------------------        

        ! compute penalty boundary condition
        call Compute_Characteristic  ! BC_Pack.f90
        call Impose_Penalty_BC       ! BC_Pack.f90
!        call Impose_Characteristic_BC


        ! Compute Flux 
        Call Flux_Cal ! State_Pack.f90
        

        ! compute Source
        call src_001(time_loc)
!        write(*,*)'src'
       
   
        ! March one intermediate Time Step
        do DDK=1, TotNum_DM 
!           write(*,*)'DDK',DDK
           ND1=PolyDegN_DM(1,DDK)
           ND2=PolyDegN_DM(2,DDK)

           ! UPGRADE IN TIME BY RK4S5
           dv1_dt(0:ND1,0:ND2,DDK)=dv1_dt(0:ND1,0:ND2,DDK) + &
                fs1(0:ND1,0:ND2,DDK)/rho(DDK)

           dv2_dt(0:ND1,0:ND2,DDK)=dv2_dt(0:ND1,0:ND2,DDK) + &
                fs2(0:ND1,0:ND2,DDK)/rho(DDK)
!           write(*,*)'Complete src'           
           
           call RK4s5_marching_2D(PolyDegN_Max(1),PolyDegN_DM(1,DDK),  &
                               v1(0,0,DDK), &
                           v1_tmp(0,0,DDK), &
                           dv1_dt(0,0,DDK))
!           write(*,*)'Complete v1'
           call RK4s5_marching_2D(PolyDegN_Max(1),PolyDegN_DM(1,DDK),  &
                               v2(0,0,DDK), &
                           v2_tmp(0,0,DDK), &
                           dv2_dt(0,0,DDK)) 
!           write(*,*)'Complete v2'
           call RK4s5_marching_2D(PolyDegN_Max(1),PolyDegN_DM(1,DDK),  &
                              T11(0,0,DDK), &
                          T11_tmp(0,0,DDK), &
                          dT11_dt(0,0,DDK))
!           write(*,*)'Complete T11'
           call RK4s5_marching_2D(PolyDegN_Max(1),PolyDegN_DM(1,DDK),  &
                              T12(0,0,DDK), &
                          T12_tmp(0,0,DDK), &
                          dT12_dt(0,0,DDK))
!           write(*,*)'Complete T12' 
           call RK4s5_marching_2D(PolyDegN_Max(1),PolyDegN_DM(1,DDK),  &
                              T22(0,0,DDK), &
                          T22_tmp(0,0,DDK), &
                          dT22_dt(0,0,DDK))
!           write(*,*)'Complete T22'
           ! assign correct value
           
        enddo !DDK
!        write(*,*)'DDK'

     enddo ! ks
     time = time + dt
!     write(*,*)time

     err_count=err_count+1
     if ( (err_count .eq. 10) .or. &
          (First_Comp .eq. 1) .or. &
          (Last_Comp .eq. 1) ) then 
        call Error_001
        err_count = 0
     endif
!!$     err_max(1:6)=0.d0
!!$     do DDK=1,TotNum_DM
!!$        ND1=PolyDegN_DM(1,DDK)
!!$        ND2=PolyDegN_DM(2,DDK)
!!$        ND3=PolyDegN_DM(3,DDK)
!!$
!!$        do k=0,ND3
!!$           do j=0,ND2 
!!$              do i=0,ND1
!!$
!!$                 x_coor=x1(i,j,k,DDK)
!!$                 y_coor=x2(i,j,k,DDK)
!!$                 z_coor=x3(i,j,k,DDK)
!!$         
!!$                 err_loc(1)= abs( E1(i,j,k,DDK) - &
!!$                               EM_demo001(x_coor,y_coor,z_coor,time,1))
!!$                 err_loc(2)= abs( E2(i,j,k,DDK) - &
!!$                               EM_demo001(x_coor,y_coor,z_coor,time,2))
!!$                 err_loc(3)= abs( E3(i,j,k,DDK) - &
!!$                               EM_demo001(x_coor,y_coor,z_coor,time,3))
!!$                 err_loc(4)= abs( H1(i,j,k,DDK) - &
!!$                               EM_demo001(x_coor,y_coor,z_coor,time,-1))
!!$                 err_loc(5)= abs( H2(i,j,k,DDK) - &
!!$                               EM_demo001(x_coor,y_coor,z_coor,time,-2))
!!$                 err_loc(6)= abs( H3(i,j,k,DDK) - &
!!$                               EM_demo001(x_coor,y_coor,z_coor,time,-3)) 
!!$                 
!!$                 err_max(1)=max(err_max(1),err_loc(1))
!!$                 err_max(2)=max(err_max(2),err_loc(2))
!!$                 err_max(3)=max(err_max(3),err_loc(3))
!!$                 err_max(4)=max(err_max(4),err_loc(4))
!!$                 err_max(5)=max(err_max(5),err_loc(5))
!!$                 err_max(6)=max(err_max(6),err_loc(6))
!!$
!!$              enddo
!!$           enddo
!!$        enddo
!!$
!!$     enddo
        
!!$     write(*,1000)time, err_max(1:6)



  enddo ! do while


!!$  call Field_Final
!!$
!!$
!!$1000 format(7e12.4)



end Program MDPS_Elastic2D
