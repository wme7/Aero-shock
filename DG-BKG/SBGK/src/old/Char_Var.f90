module Char_Var
  implicit none
  !Declare BC type boundary conditions
  integer, save, allocatable:: BC_Type(:,:)
  integer, save, allocatable:: Connect_Table(:,:,:)

  !Declare Module Variables
  real(kind=8), save, allocatable::   RFV(:,:,:,:)
  real(kind=8), save, allocatable::     R(:,:,:,:)
  real(kind=8), save, allocatable::   RBC(:,:,:,:)
   
  ! Penalty Variables
  real(kind=8), parameter :: tau_scale=4.d0
  real(kind=8), save, allocatable:: tau_edge(:,:,:)
  

  ! S, SI and BDO opeartor
  real(kind=8), save, allocatable::   S(:,:,:,:,:)
  real(kind=8), save, allocatable::  SI(:,:,:,:,:)
  real(kind=8), save, allocatable:: BDO(:,:,:,:)

contains
!=================================================================
  subroutine alloc_mem_Char_Var(N_max,Num_Domain)
    implicit none
    integer:: N_max
    integer:: Num_Domain

    integer:: ierr


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

    allocate (   RFV(5,0:N_max,4,Num_Domain), &
                   R(5,0:N_max,4,Num_Domain), &
                 RBC(5,0:N_max,4,Num_Domain), &
                 stat=ierr )

    if (ierr .ne. 0) then
       write(*,*)'Error message from Char_Var.f90:'
       write(*,*)'Cannot allocate memory for Characteristic Variables'
       write(*,*)'Abort!'
       stop
    endif

    RFV=0.d0
      R=0.d0
    RBC=0.d0

    allocate (tau_edge(0:N_max,4,Num_Domain), stat=ierr)
    if (ierr .ne. 0) then 
       write(*,*)'Error message from Char_Var.f90:'
       write(*,*)'Cannot allocate memory for Penalty Parameters'
       write(*,*)'Abort!'
       stop
    endif

    tau_edge=0.d0

    allocate (   S(5,5,0:N_max,4,Num_Domain), &
                SI(5,5,0:N_max,4,Num_Domain), &
                 BDO(5,0:N_max,4,Num_Domain), &
               stat = ierr ) 
    if (ierr .ne. 0) then 
       write(*,*)'Error message from Char_Var.f90:'
       write(*,*)'Cannot allocate memory for Penalty Parameters'
       write(*,*)'Abort!'
       stop
    endif

    S=0.d0; SI=0.d0; BDO=0.d0

    return

  end subroutine alloc_mem_Char_Var
!=================================================================
end module Char_Var
