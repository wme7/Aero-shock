module NorVec_Var
  implicit none
  !Declare Module Variables
  real(kind=8), save, allocatable:: NorVec_x1(:,:,:)
  real(kind=8), save, allocatable:: NorVec_x2(:,:,:)
  real(kind=8), save, allocatable:: NorVec_mg(:,:,:)

contains
!=================================================================
  subroutine alloc_mem_norvec(N_max,Num_Domain)
    implicit none
    integer:: N_max
    integer:: Num_Domain

    integer:: ierr


    ! subroutine begins
    allocate (NorVec_x1(0:N_max,4,Num_Domain), &
              NorVec_x2(0:N_max,4,Num_Domain), &
              NorVec_mg(0:N_max,4,Num_Domain), stat=ierr )

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
end module NorVec_Var
