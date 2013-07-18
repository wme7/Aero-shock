Module material_Var
  implicit none

  real(kind=8), save, allocatable:: rho(:)           ! density
  real(kind=8), save, allocatable:: Lame_mu(:)       ! Lame_eps(TotNum_DM)
  real(kind=8), save, allocatable:: Lame_lambda(:)   ! Lame_lambda(TotNum_DM)
  real(kind=8), save, allocatable:: H_inv(:,:,:,:,:) ! 
  
contains
  subroutine alloc_mem_material(Deg_Max, Num_Domain)
    implicit none 
    ! Declare subroutine arguments
    integer:: Deg_Max, Num_Domain
    
    ! Declare local arguments
    integer:: ierr
    
    ! subroutine begins
    allocate(        rho(1:Num_Domain), &
                 lame_mu(1:Num_Domain), &
             lame_lambda(1:Num_Domain), stat=ierr)

    if (ierr .ne. 0) then
       write(*,*)'Message from Material.f90'
       write(*,*)'Can not allocate memory for material variables'
       write(*,*)'Abort!'
       stop
    endif

    rho=0.d0
    lame_mu=0.d0
    lame_lambda=0.d0 

    allocate(H_inv(5,5,0:Deg_Max,4,Num_Domain), stat=ierr)
    if (ierr .ne. 0) then
       write(*,*)'Message from Material.f90'
       write(*,*)'Can not allocate memory for symmetrizer variables'
       write(*,*)'Abort!'
       stop
    endif

    H_inv=0.d0
 

  end subroutine alloc_mem_material
       
end module material_Var
