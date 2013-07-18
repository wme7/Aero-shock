!===============================================================!
!                                                               !
! This is a module defining variables for 2D Elastic Waves      !
!                                                               !
!===============================================================!
module State_Var  !Elastic Waves
implicit none

real(kind=8), save, allocatable::  v1(:,:,:) !(0:Deg_Max,1:TotNum_DM)
real(kind=8), save, allocatable::  v2(:,:,:) 
real(kind=8), save, allocatable:: T11(:,:,:) !(0:Deg_Max,1:TotNum_DM)
real(kind=8), save, allocatable:: T12(:,:,:)
real(kind=8), save, allocatable:: T22(:,:,:) !(0:Deg_Max,1:TotNum_DM)

real(kind=8), save, allocatable::  v1_tmp(:,:,:) !(0:Deg_Max,1:TotNum_DM,0:level)
real(kind=8), save, allocatable::  v2_tmp(:,:,:)
real(kind=8), save, allocatable:: T11_tmp(:,:,:) !(0:Deg_Max,1:TotNum_DM)
real(kind=8), save, allocatable:: T12_tmp(:,:,:)
real(kind=8), save, allocatable:: T22_tmp(:,:,:) !(0:Deg_Max,1:TotNum_DM)

real(kind=8), save, allocatable::  dv1_dt(:,:,:) !(0:Deg_Max,1:TotNum_DM,0:level)
real(kind=8), save, allocatable::  dv2_dt(:,:,:)
real(kind=8), save, allocatable:: dT11_dt(:,:,:)
real(kind=8), save, allocatable:: dT12_dt(:,:,:) !(0:Deg_Max,1:TotNum_DM,0:level)
real(kind=8), save, allocatable:: dT22_dt(:,:,:)

real(kind=8), save, allocatable::  fs1(:,:,:)
real(kind=8), save, allocatable::  fs2(:,:,:)

real(kind=8), save, allocatable:: dqdx(:,:) !temporal storage 
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
    integer:: Deg_Max(2)                ! Pass in the PolyDegN_Max
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
    allocate(      v1(0:ND1,0:ND2,1:TotNumDomain), &
                   v2(0:ND1,0:ND2,1:TotNumDomain), &
                  T11(0:ND1,0:ND2,1:TotNumDomain), &
                  T12(0:ND1,0:ND2,1:TotNumDomain), &
                  T22(0:ND1,0:ND2,1:TotNumDomain), &
               v1_tmp(0:ND1,0:ND2,1:TotNumDomain), &
               v2_tmp(0:ND1,0:ND2,1:TotNumDomain), &
              T11_tmp(0:ND1,0:ND2,1:TotNumDomain), &
              T12_tmp(0:ND1,0:ND2,1:TotNumDomain), &
              T22_tmp(0:ND1,0:ND2,1:TotNumDomain), &
               dv1_dt(0:ND1,0:ND2,1:TotNumDomain), &
               dv2_dt(0:ND1,0:ND2,1:TotNumDomain), &
              dT11_dt(0:ND1,0:ND2,1:TotNumDomain), &
              dT12_dt(0:ND1,0:ND2,1:TotNumDomain), &
              dT22_dt(0:ND1,0:ND2,1:TotNumDomain), &
                  fs1(0:ND1,0:ND2,1:TotNumDomain), &
                  fs2(0:ND1,0:ND2,1:TotNumDomain), &    
              stat=ierr)

    if ( ierr .ne. 0 ) then
       write(*,*)'Cannot allocate memory for state variable q'
       write(*,*)'Abort!'
       stop
    endif

     v1=0.d0;  v2=0.d0; 
    T11=0.d0; T12=0.d0; T22=0.d0;

     v1_tmp=0.d0;  v2_tmp=0.d0;
    T11_tmp=0.d0; T12_tmp=0.d0; T22_tmp=0.d0;

     dv1_dt=0.d0;  dv2_dt=0.d0;
    dT11_dt=0.d0; dT12_dt=0.d0; dT22_dt=0.d0;

    fs1=0.d0; fs2=0.d0

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
