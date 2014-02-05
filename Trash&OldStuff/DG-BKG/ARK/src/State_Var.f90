!===============================================================!
!                                                               !
! This is a module defining variables for 2D Elastic Waves      !
!                                                               !
!===============================================================!
module State_Var  !Elastic Waves
implicit none
!--------------<Cartesian Coordinate>---------------------
real(kind=8), save, allocatable::  uij(:,:,:,:,:) !(0:Deg_Max,1:TotNum_DM)
real(kind=8), save, allocatable::  uij_tmp(:,:,:,:,:) !(0:Deg_Max,1:TotNum_DM,0:level)
real(kind=8), save, allocatable::  duij_dt(:,:,:,:,:) !(0:Deg_Max,1:TotNum_DM,0:level)
real(kind=8), save, allocatable::  u(:,:,:) !(0:Deg_Max,1:TotNum_DM)
real(kind=8), save, allocatable::  u_tmp(:,:,:) !(0:Deg_Max,1:TotNum_DM,0:level)
real(kind=8), save, allocatable::  du_dt(:,:,:) !(0:Deg_Max,1:TotNum_DM,0:level)

real(kind=8), save, allocatable:: dqdx(:,:) !temporal storage 
!real(kind=8), save, allocatable:: GH(:,:),GHW(:,:)
real(kind=8), save::a_vec(2)
integer :: I_prob !1: Old; 2:Periodic; 3:Neumann BC
integer :: I_case ! Configuration
character(len=2) :: NumCase

!integer:: IGH=20

!---------------------------------------------------------------------------
! Declare state variable q
! Where q is the unknown of the following PDE
!
!                q_t = L q
!

!================================================================================
!
!    [Name]    :: q(:,:,:)
!     ^^^^
!    [Size]    :: (0:Degree1_Max,0:Degree2_Max,1:TotNum_DM)
!     ^^^^
!    [Purpose] :: To Store the value of the PDE State Variable q 
!     ^^^^^^^
!    [Detail]  ::
!     ^^^^^^
!
!================================================================================

contains

!------------------------------------------------------------------------------

  subroutine Init_State_Variables(Deg_Max,TotNumDomain)
    implicit none
    ! declare subroutine arguments
    integer:: Deg_Max(2)                ! Pass in the Nxi_max
    integer:: TotNumDomain
   

    integer:: ND1,ND2 
    integer:: ierr
   
    ! Allocate memory for the PDE state vector (unknown variable to be solve)
    ! by subroutine alloc_mem_state_variables(Degree_Max,TotNumDomain,Storage)
    ! provided by State_Var Module.
    ! 
    ! alloc_mem_state_variables(Degree_Max,TotNumDomain,Storage) !State_Var
    ! storage = 0 => one level 
    ! storage = 1 => two level
    ! storage = 2 => three level
    
    ND1=Deg_Max(1); ND2=Deg_Max(2);
    allocate(      u(0:ND1,0:ND2,1:TotNumDomain), &
               u_tmp(0:ND1,0:ND2,1:TotNumDomain), &
               du_dt(0:ND1,0:ND2,1:TotNumDomain), &
              stat=ierr)

    if ( ierr .ne. 0 ) then
       write(*,*)'Cannot allocate memory for state variable q'
       write(*,*)'Abort!'
       stop
    endif

     u=0.d0
     u_tmp=0.d0
     du_dt=0.d0

    ND1=maxval(Deg_Max(1:2))
    allocate(dqdx(0:ND1,0:ND1), stat=ierr)
    if (ierr .ne. 0) then 
       write(*,*)' Message from State_Var.f90'
       write(*,*)' Cannot allocate memory for variable dqdx'
       write(*,*)' Abort!'
       stop
    endif 
    dqdx=0.d0

    return 
    
  end subroutine Init_State_Variables

!-------------------------------------------------------------------------------  
  
end module State_Var
