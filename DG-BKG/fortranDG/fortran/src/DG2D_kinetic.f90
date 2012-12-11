Program DG2D_kinetic
use universal_const    ! universal.f90
use Legendre           ! Legendre.f90
use RK_Var
use MD2D_Grid
use State_Var
use Metric_Var
use NorVec_Var
use Kinetic_Var
use Char_Var
implicit none

!real(kind=8) :: time_start,time_end,sqTDM,ttau_scale,dx,CCFL
!real(kind=8) :: a_vec(1:2),
real(kind=8) :: error_L2,error_max,dtt
real(kind=8) :: a_max
integer :: ierr,q,lid,rk, sqrTDM, sqrTDM2,sqrTDM4, sqrTDM8

! --- Program Begin ---
! --- Initialize stage ---
! --- Initialize universal constants ---

  call init_universal_const				      !-- universal.f90

  call Input_Domain_Grid_Parameters			      !-- MD2D_Grid.f90



  lid=22
  open(22,file='./Parameter/in.in',form='formatted', status='old')
  read(lid,*) !==============================================!
  read(lid,*) Time_final
  read(lid,*) a_vec(1)
  read(lid,*) a_vec(2)
  read(lid,*) !==============================================!
 
! --- Initialize Domain Parameter ---
! --- Initialize Legendre Pseudospectral Module ---
! --- 1. Gauss-Lobatto-Legendre Grid ---
! --- 2. Gauss-Lobatto-Legendre Quadurature Weights ---
! --- 3. Legenedre Differentiation Matrix ---

  call alloc_mem_u(PDeg1,PDeg2,PND1,PND2,TotNum_DM)           !-- Legendre.f90

  call Init_LGL_Diff_Matrix_23D(PND1,PDeg1)      	      !-- Legendre.f90 
  
  call Init_Physical_Grid_2D 				      !-- Grid2D_Pack.f90

  call Init_intept_Grid_2D                                    !-- Grid2D_Pack.f90

  call Init_Metric					      !-- Metric_Pack.f90

  call Init_Normal_Vector  			   	      !-- NorVec_Pack.f90

  call alloc_mem_Char_Var               	              !-- Char_Var.f90
   
  call Init_Penalty_Parameters(ttau_scale)                    !-- BC_Pack.f90

  call Init_Kinetic_Variables(PolyNodes,TotNum_DM)              !-- Kinetic_Var

! --- set up initial condition in DG2D_Initial ---
! --- Initialize Time Integration ---

!  call DG2D_Initial(a_vec,Time_final)                         !-- NorVec_Pack.f90
!  call DG2D_Initial(Time_final)                         !-- NorVec_Pack.f90
  call DG2D_Kinetic_Initial(Time_final)                         !-- Kinetic_Pack.f90

  a_max=maxval(GH)
  write(6,*) "a_max",a_max
  sqTDM = sqrt(dble(TotNum_DM))

  sqrTDM=nint(sqTDM)
  sqrTDM2=sqrTDM/2
  sqrTDM4=sqrTDM/4
  sqrTDM8=sqrTDM/8

  ttau_scale = dble(2/sqTDM)
  dx = 2/sqTDM  
  CCFL = dble(2*(2*PDeg1+1))*a_max
  CCFL = dble(1/CCFL) !*0.9d0
  dtt = dx*CCFL
  rk = PDeg1 + 1


!  call DG2D_Edge(a_vec)                                       !-- NorVec_Pack.f90
  call DG2D_Edge()                                       !-- NorVec_Pack.f90

  call Init_SSPRK(rk)                                         !-- RK_Var.f90

  write(*,*)"dt = ",dtt,"Deg=",PDeg1

! --- Time  advancing ---
  
  time = 0.0
  call cpu_time(time_start)

  do while (time .lt. Time_final)
     if ( (time+dtt) .ge. Time_final) then 
        dtt = Time_final-time        
     endif 

     F_alt = F_new
     F_new = alpha(1)*F_new
       CALL CALDOM

     do q = 1,rk

        CALL EQUILIBRIUM

        call DG2D_flux                                        ! Diff.Pack.f90
 
        call SSPRK_marching_2D(q,rk,dtt)                      ! RK_Var.f90
!write(6,*) "Q=",q
        CALL EQUILIBRIUM

       CALL CALDOM 
     enddo


!        write(6,*) "Dm 28",x1(PND1,PND2,28),x2(PND1,PND2,28),"Density",R_loc(PND1,PND2,28)
!        write(6,*) "Dm 29",x1(0,PND2,29),x2(0,PND2,29),"Density",R_loc(0,PND2,29)
!        write(6,*) "Dm 36",x1(PND1,0,36),x2(PND1,0,36),"Density",R_loc(PND1,0,36)
!        write(6,*) "Dm 37",x1(0,0,37),x2(0,0,37),"Density",R_loc(0,0,37)

DDK=sqrTDM*(sqrTDM2-1)+sqrTDM2
DDK=sqrTDM*(3*sqrTDM8-1)+sqrTDM8*3
write(*,777) DDK,x1(PND1,PND2,DDK),x2(PND1,PND2,DDK),R_loc(PND1,PND2,DDK)
write(*,777) DDK+1,x1(0,PND2,DDK+1),x2(0,PND2,DDK+1),R_loc(0,PND2,DDK+1)
write(*,777) DDK+sqrTDM,x1(PND1,0,DDK+sqrTDM),x2(PND1,0,DDK+sqrTDM),R_loc(PND1,0,DDK+sqrTDM)
write(*,777) DDK+sqrTDM+1,x1(0,0,DDK+sqrTDM+1),x2(0,0,DDK+sqrTDM+1),R_loc(0,0,DDK+sqrTDM+1)
777 format (1X,'Domain #:',I7,2X,'at x=:',F7.4,2X,'at y=:',F7.4,2X,'Density=',F7.4)

     time = time + dtt
write(6,*) "time=",time

  enddo ! do while

  call cpu_time(time_end)

  call DG_compute_error(error_L2,error_max,Time_final) !,a_vec)  ! Diff.Pack.f90

  write(*,*)"error_L2=",error_L2,"error_max=",error_max
  write(*,"(a,1x,1f12.6,1x,a)") "Fortran, cost",time_end-time_start,"s"
  
end Program DG2D_kinetic