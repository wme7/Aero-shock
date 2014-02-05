subroutine DG_Init_Filter_Oper(Deg_Max)
  use Filter_Var
  use Legendre
  implicit none
  integer :: Deg_Max
  
  integer :: ierr, i, j ,l ,k, ND1
  
  !allocate variables for LGL_grid and initialize it
  !   ND1=PolyDegN_Max(1); ND2=PolyDegN_Max(2); ND_Max=max(ND1,ND2)
  ! Take the max of PolyDegreeN_max as leading dimention to set up
  ! the legendre gauss labatto grid and its differential matrix for
  ! each degree
  !   call Alloc_Mem_LG_grid(ND_Max)
  !   call Alloc_Mem_Leg_Diff_Mat(ND_Max)
  
  call Read_In_Filter_Profile
    allocate ( filter_sigma(0:Deg_Max))
!  call alloc_mem_filter_variables(Deg_Max)

!!$
!!$  ! set up LGL grid and differential matrix computational space
!!$
!!$  
!!$  do ND1=1,Deg_Max  ! compute filter for each degree
     
     call filter_profile(Deg_Max)
     write(6,*) "filter",filter_sigma
!!$     call init_LGL_grids(ND1)	! find the LGL grid for each degree
!!$     ! LG_weights(0:ND1) is also computed
!!$     
!!$     do l=0,ND1
!!$        
!!$        do i=0,ND1
!!$           gamma_N(i)=2.d0/(2.d0*dble(i)+1.d0)
!!$        enddo
!!$        gamma_N(ND1)=2.d0/dble(ND1)
!!$        
!!$        do k=0,ND1
!!$           call VALEPO(l,LG_Grids(k),PnAtgrid(k),dPnAtgrid(k),d2PnAtgrid(k))
!!$        enddo
!!$        
!!$        do i=0,ND1
!!$           ValPolyAtGrid(i,l)=PnAtgrid(i)
!!$!           	     write(*,*)'here',i,PnAtgrid(i) !,LG_weights(i)
!!$        enddo
!!$!        	  pause
!!$     enddo
     
!!$     Q(0:ND1,0:ND1)=ValPolyAtGrid(0:ND1,0:ND1)
!!$     Qt(0:ND1,0:ND1)=transpose(ValPolyAtGrid(0:ND1,0:ND1))
!!$     
!!$     do j=0,ND1
!!$        do i=0,ND1
!!$           
!!$           Q(i,j)  = filter_sigma(j)*Q(i,j)/gamma_N(j)
!!$           Qt(i,j) = LG_weights(j)*Qt(i,j)
!!$           !	     write(*,*)i,j,Q(i,j),Qt(i,j)
!!$           
!!$        enddo
!!$     enddo
!!$     
!!$     filter_oper_xi1(0:ND1,0:ND1,ND1)=matmul(Q(0:ND1,0:ND1),Qt(0:ND1,0:ND1))
!!$     filter_oper_xi2(0:ND1,0:ND1,ND1)= &
!!$          transpose(filter_oper_xi1(0:ND1,0:ND1,ND1))
     
!     do j=0,ND1
!        do i=0,ND1
!           write(*,*)i,j,filter_oper_xi1(i,j,ND1)
!        enddo
!     enddo
!     pause
!!$  enddo
!!$  
!!$  call Diff_x_Filter(Deg_Max)
!!$  !   deallocate()
  
  return
  
end subroutine DG_Init_Filter_Oper

subroutine Init_Filter_Oper(Deg_Max)
  use Filter_Var
  use Legendre
  implicit none
  integer :: Deg_Max
  
  integer :: ierr, i, j ,l ,k, ND1
  
  !allocate variables for LGL_grid and initialize it
  !   ND1=PolyDegN_Max(1); ND2=PolyDegN_Max(2); ND_Max=max(ND1,ND2)
  ! Take the max of PolyDegreeN_max as leading dimention to set up
  ! the legendre gauss labatto grid and its differential matrix for
  ! each degree
  !   call Alloc_Mem_LG_grid(ND_Max)
  !   call Alloc_Mem_Leg_Diff_Mat(ND_Max)
  
  call Read_In_Filter_Profile

  call alloc_mem_filter_variables(Deg_Max)


  ! set up LGL grid and differential matrix computational space

  
  do ND1=1,Deg_Max  ! compute filter for each degree
     
     call filter_profile(ND1)
     
     call init_LGL_grids(ND1)	! find the LGL grid for each degree
     ! LG_weights(0:ND1) is also computed
     
     do l=0,ND1
        
        do i=0,ND1
           gamma_N(i)=2.d0/(2.d0*dble(i)+1.d0)
        enddo
        gamma_N(ND1)=2.d0/dble(ND1)
        
        do k=0,ND1
           call VALEPO(l,LG_Grids(k),PnAtgrid(k),dPnAtgrid(k),d2PnAtgrid(k))
        enddo
        
        do i=0,ND1
           ValPolyAtGrid(i,l)=PnAtgrid(i)
!           	     write(*,*)'here',i,PnAtgrid(i) !,LG_weights(i)
        enddo
!        	  pause
     enddo
     
     Q(0:ND1,0:ND1)=ValPolyAtGrid(0:ND1,0:ND1)
     Qt(0:ND1,0:ND1)=transpose(ValPolyAtGrid(0:ND1,0:ND1))
     
     do j=0,ND1
        do i=0,ND1
           
           Q(i,j)  = filter_sigma(j)*Q(i,j)/gamma_N(j)
           Qt(i,j) = LG_weights(j)*Qt(i,j)
           !	     write(*,*)i,j,Q(i,j),Qt(i,j)
           
        enddo
     enddo
     
     filter_oper_xi1(0:ND1,0:ND1,ND1)=matmul(Q(0:ND1,0:ND1),Qt(0:ND1,0:ND1))
     filter_oper_xi2(0:ND1,0:ND1,ND1)= &
          transpose(filter_oper_xi1(0:ND1,0:ND1,ND1))
     
!     do j=0,ND1
!        do i=0,ND1
!           write(*,*)i,j,filter_oper_xi1(i,j,ND1)
!        enddo
!     enddo
!     pause
  enddo
  
  call Diff_x_Filter(Deg_Max)
  !   deallocate()
  
  return
  
end subroutine Init_Filter_Oper

subroutine Diff_x_Filter(Deg_Max)
  use Legendre
  use Filter_Var
  implicit none
  integer::Deg_Max
  integer::ND

  do ND=1,Deg_Max
     Q(0:ND,0:ND)=matmul(Diff_xi1(0:ND,0:ND,ND),&
                  filter_oper_xi1(0:ND,0:ND,ND))

     Diff_xi1(0:ND,0:ND,ND)=Q(0:ND,0:ND)
  enddo

  do ND=1,Deg_Max
     Q(0:ND,0:ND)=matmul(filter_oper_xi2(0:ND,0:ND,ND),&
                                Diff_xi2(0:ND,0:ND,ND))
       
     Diff_xi2(0:ND,0:ND,ND)=Q(0:ND,0:ND)
  enddo
  

end subroutine Diff_x_Filter

subroutine Smoothing_Field
use DG_Var
use Filter_Var
use State_Var
use MD3D_Grid
implicit none
integer::l1,l2,l3
real(kind=8):: o

  do DDK=1, TotNum_DM 

     ND1=DOF_DM(1,DDK)
     ND2=DOF_DM(2,DDK)
do l1=0,ND1
do l2=0,ND2
o=filter_sigma(l1)*filter_sigma(l2)*filter_sigma(l3)
E1(l1,l2,l3,DDK)=E1(l1,l2,l3,DDK)*o
E2(l1,l2,l3,DDK)=E2(l1,l2,l3,DDK)*o
E3(l1,l2,l3,DDK)=E3(l1,l2,l3,DDK)*o
enddo
enddo

enddo

end subroutine Smoothing_Field

!!$subroutine Smoothing_Field(PML_On_Off)
!!$  use Ctrl_Domain
!!$  use Domain_Variables
!!$  use PML_Variables
!!$  use filter
!!$  integer PML_On_Off
!!$  
!!$  integer side_num
!!$  
!!$  if (PML_On_Off .eq. 0 ) return
!!$  
!!$  ! filter solution if needed
!!$  do PML_DDK=1,TotNumPMLDM
!!$     
!!$     DDK=PMLDM_2_GlobalDM_Table(PML_DDK)
!!$     ND1=DM_PolyDegN(1,DDK); ND2=DM_PolyDegN(2,DDK)
!!$     ! filter the computed solution
!!$     do side_num=1,4
!!$        if (DMBnd_RBC_type(side_num,DDK) .eq. -1) then
!!$           
!!$           if ( mod(side_num,2) .eq. 0 )  then
!!$              ! filter in z
!!$              Q=0.d0
!!$              Q(0:ND1,0:ND2)=matmul(smoothing_oper_xi1(0:ND1,0:ND1,ND1),&
!!$                   DMGrd_Ey(0:ND1,0:ND2,DDK))
!!$              DMGrd_Ey(0:ND1,0:ND2,DDK)=Q(0:ND1,0:ND2)
!!$              
!!$              Q=0.d0
!!$              Q(0:ND1,0:ND2)=matmul(smoothing_oper_xi1(0:ND1,0:ND1,ND1),&
!!$                   DMGrd_Hx(0:ND1,0:ND2,DDK))
!!$              DMGrd_Hx(0:ND1,0:ND2,DDK)=Q(0:ND1,0:ND2)
!!$              
!!$              Q=0.d0
!!$              Q(0:ND1,0:ND2)=matmul(smoothing_oper_xi1(0:ND1,0:ND1,ND1),&
!!$                   DMGrd_Hz(0:ND1,0:ND2,DDK))
!!$              DMGrd_Hz(0:ND1,0:ND2,DDK)=Q(0:ND1,0:ND2)
!!$              
!!$           endif
!!$           
!!$           if ( mod(side_num,2) .eq. 1)  then
!!$              ! filter in x
!!$           endif
!!$           
!!$        endif
!!$     enddo
!!$     
!!$  enddo
!!$end subroutine Smoothing_Field
