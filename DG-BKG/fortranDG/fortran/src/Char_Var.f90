module Char_Var
  implicit none
  !Declare BC type boundary conditions
  integer, save, allocatable:: BC_Type(:,:)
  integer, save, allocatable:: Connect_Table(:,:,:)

!  !Declare Module Variables
  real(kind=8), save, allocatable::     R(:,:,:,:)
   
  ! Penalty Variables
  real(kind=8), parameter :: tau_scale=2.d0
  real(kind=8), save, allocatable:: tau_edge(:,:,:)
  
 contains
!=================================================================
  subroutine alloc_mem_Char_Var
    use MD2D_Grid
    implicit none

    integer:: ierr,N_max,Num_Domain

    N_max = PND1
    Num_Domain = TotNum_DM

    ! subroutine begins
    allocate ( BC_Type(4,Num_Domain), &
               Connect_Table(2,4,Num_Domain), stat=ierr )
    if ( ierr .ne. 0 ) then 
       write(*,*)'Error message from Char_Var.f90:'
       write(*,*)'Cannot allocate memory for BC_Type Variables'
       write(*,*)'Abort!'
       stop
    endif
    BC_Type=0

    allocate (     R(1,0:N_max,4,Num_Domain), stat=ierr )

    if (ierr .ne. 0) then
       write(*,*)'Error message from Char_Var.f90:'
       write(*,*)'Cannot allocate memory for Characteristic Variables'
       write(*,*)'Abort!'
       stop
    endif

      R=0.d0

    allocate (tau_edge(0:N_max,4,Num_Domain), stat=ierr)
    if (ierr .ne. 0) then 
       write(*,*)'Error message from Char_Var.f90:'
       write(*,*)'Cannot allocate memory for Penalty Parameters'
       write(*,*)'Abort!'
       stop
    endif

    tau_edge=0.d0

    return

  end subroutine alloc_mem_Char_Var
!=================================================================
end module Char_Var
