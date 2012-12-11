Module Metric_Var !Metric 2D
  implicit none

  real(kind=8), save, allocatable :: DM_Metric(:,:,:) 
  real(kind=8), save, allocatable :: Jacobian(:,:,:)   !((0:N_max)^3,TotNum_DM)
  real(kind=8), save, allocatable :: dxi1(:,:,:)      !((0:N_max)^3,TotNum_DM) 
  real(kind=8), save, allocatable :: dxi2(:,:,:)      

  ! Declare variable dxi/dx
  real(kind=8), save, allocatable :: dxi1_dx1(:,:,:)
  real(kind=8), save, allocatable :: dxi1_dx2(:,:,:)
  real(kind=8), save, allocatable :: dxi2_dx1(:,:,:)
  real(kind=8), save, allocatable :: dxi2_dx2(:,:,:)

  ! Declare variable dx/dxi
  real(kind=8), save, allocatable :: dx1_dxi1(:,:,:)
  real(kind=8), save, allocatable :: dx1_dxi2(:,:,:)
  real(kind=8), save, allocatable :: dx2_dxi1(:,:,:)
  real(kind=8), save, allocatable :: dx2_dxi2(:,:,:)
  
  ! Declare Grid distortion variables
  real(kind=8), save, allocatable :: dtrans(:,:,:)        !(0:N_max,TotNumDomain)

contains
  
  subroutine alloc_mem_metric_variables
    use MD2D_Grid
    implicit none
    integer:: Deg1_Max, Deg2_Max, TotNumDomain

    integer:: ierr
    Deg1_Max = PND1
    Deg2_Max = PND2
    TotNumDomain = TotNum_DM
    
    allocate( &
         DM_Metric(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
           Jacobian(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
              dxi1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
              dxi2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
          dxi1_dx1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
          dxi1_dx2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
          dxi2_dx1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
          dxi2_dx2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
          dx1_dxi1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
          dx1_dxi2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
          dx2_dxi1(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
          dx2_dxi2(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), &
            dtrans(0:Deg1_Max,0:Deg2_Max,1:TotNumDomain), stat=ierr)

    if (ierr .ne. 0) then 
       write(*,*)'Cannot allocate memory for Metric Variables'
       write(*,*)'Abort!'
       stop
    endif

    DM_Metric=0.d0
    Jacobian=0.d0

    dxi1=0.d0; dxi2=0.d0;

    dxi1_dx1=0.d0; dxi1_dx2=0.d0;
    dxi2_dx1=0.d0; dxi2_dx2=0.d0;

    dx1_dxi1=0.d0; dx1_dxi2=0.d0;
    dx2_dxi1=0.d0; dx2_dxi2=0.d0;

    dtrans=0.d0

    return 
  end subroutine alloc_mem_metric_variables


end module Metric_Var
