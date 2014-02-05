module NorVec_Var
  implicit none
  !Declare Module Variables

  real(kind=8),save,allocatable :: NorVec_x1(:,:,:),NorVec_x2(:,:,:),NorVec_mg(:,:,:)
  real(kind=8),save,allocatable :: C1x_Edge1(:,:),C2x_Edge1(:,:)
  real(kind=8),save,allocatable :: C1x_Edge3(:,:),C2x_Edge3(:,:)
  real(kind=8),save,allocatable :: Cy_Edge1(:,:),Cy_Edge2(:,:),Cy_Edge3(:,:),Cy_Edge4(:,:)
  real(kind=8),save,allocatable :: E_gen(:,:,:,:)
!  real(kind=8),save,allocatable :: E(:,:,:)
!real(kind=8), save, allocatable:: GH(:,:,:),GHW(:,:,:)
!integer:: IGH=20

contains
!=================================================================
  subroutine alloc_mem_norvec
    use MD2D_Grid
    implicit none
    integer:: N_max
    integer:: Num_Domain

    integer:: ierr
    N_max = PND1
    Num_Domain = TotNum_DM

    ! subroutine begins
    allocate (NorVec_x1(0:N_max,1:4,1:Num_Domain), &
              NorVec_x2(0:N_max,1:4,1:Num_Domain), &
              NorVec_mg(0:N_max,1:4,1:Num_Domain), stat=ierr )

    if (ierr .ne. 0) then
       write(*,*)'Error message from NorVec_Var.f90:'
       write(*,*)'Cannot allocate memory for Normal Vector Variables'
       write(*,*)'abort!'
    endif

    NorVec_x1=0.d0
    NorVec_x2=0.d0
    NorVec_mg=0.d0
   
    return
  end subroutine alloc_mem_norvec

!=================================================================

  subroutine alloc_mem_DG2D_Edge
    use Legendre
    use MD2D_Grid
    use State_Var
    implicit none
    integer:: ierr

    allocate( C1x_Edge1(0:PDeg1,0:PDeg2),C2x_Edge1(0:PDeg1,0:PDeg2), &
            C1x_Edge3(0:PDeg1,0:PDeg2),C2x_Edge3(0:PDeg1,0:PDeg2), &
            Cy_Edge1(1:1,0:PDeg1),Cy_Edge2(0:PDeg2,1:1), &
            Cy_Edge3(1:1,0:PDeg1),Cy_Edge4(0:PDeg2,1:1), &
            stat=ierr )
!            E(1:2,1:4,1:TotNum_DM), stat=ierr )

    if (ierr .ne. 0) then
       write(*,*)'Cannot allocate memory for DG_edge'
       stop
    endif

   allocate(E_gen(1:2,1:4,1:IGH,1:IGH), stat=ierr )
    if (ierr .ne. 0) then
       write(*,*)'Cannot allocate memory for DG_edge'
       stop
    endif

!   allocate(GH(1:IGH,1:IGH),GHW(1:IGH,1:IGH), stat=ierr)
!    if (ierr .ne. 0) then
!       write(*,*)' Message from State_Var.f90'
!       write(*,*)' Cannot allocate memory for variable GH'
!       write(*,*)' Abort!'
!       stop
!    endif


  end subroutine alloc_mem_DG2D_Edge
!=================================================================

end module NorVec_Var
