module State_Var
  implicit none
  !Variables for the equation
! For Implicit-Explicit RK
   real(kind=8), allocatable:: F_s(:,:,:,:) !(NV,nx,pp,stage)
   real(kind=8), allocatable:: F_ns(:,:,:,:)!(NV,nx,pp,stage)

   real(kind=8), allocatable:: F(:,:,:) !(NV,nx,pp)
   real(kind=8), allocatable:: FEQ(:,:,:) !(NV,nx,pp)
   real(kind=8), allocatable:: Fold(:,:,:) !(NV,nx,pp)
   real(kind=8), allocatable:: F_new(:,:,:) !(NV,nx,pp)
   real(kind=8), allocatable:: FS(:,:,:) !(NV,nx,pp)
! Physical Quantities
   real(kind=8), allocatable:: x(:,:) !(nx,pp)
   real(kind=8), allocatable:: Rs(:,:) !(nx,pp)
   real(kind=8), allocatable:: Ps(:,:) !(nx,pp)
   real(kind=8), allocatable:: Us(:,:) !(nx,pp)
   real(kind=8), allocatable:: Ts(:,:) !(nx,pp)
   real(kind=8), allocatable:: Zs(:,:) !(nx,pp)
   real(kind=8), allocatable:: ET(:,:) !(nx,pp)
   real(kind=8), allocatable:: AV(:,:) !(nx,pp)
   real(kind=8), allocatable:: VIS(:,:) !(nx,pp)

! Temperal Variables for DG scheme
!   real(kind=8), allocatable:: FR(:,:) !(pp,1)
!   real(kind=8), allocatable:: FU(:,:) !(pp,1)
!   real(kind=8), allocatable:: FC(:,:) !(pp,1)
!   real(kind=8), allocatable:: FN(:,:) !(pp,1)
! Discrete Ordinate Method
   real(kind=8), allocatable:: V(:) !(NV)
   real(kind=8), allocatable:: VW(:) !(NV)
   integer:: NV,IT,BC_type, nx_g,pp_g
   real(kind=8) :: dx_g, dt_g
contains 

  !------------------------------------------------------------------------------

  !------------------------------------------------------------------------------
  ! Initialize of the Variables
  subroutine Init_Var(NV,nnx,pp,rks)
    implicit none
    integer :: NV,nnx,pp,rks
    integer i,ierr
    
    allocate(F(1:NV,1:nnx,1:pp), &
FEQ(1:NV,1:nnx,1:pp), &
Fold(1:NV,1:nnx,1:pp), &
F_new(1:NV,1:nnx,1:pp), &
FS(1:NV,1:nnx,1:pp), &
Rs(1:nnx,1:pp), &
Ps(1:nnx,1:pp), &
Us(1:nnx,1:pp), &
Ts(1:nnx,1:pp), &
Zs(1:nnx,1:pp), &
ET(1:nnx,1:pp), &
AV(1:nnx,1:pp), &
VIS(1:nnx,1:pp), &
x(1:nnx,1:pp), &
F_s(1:NV,1:nnx,1:pp,rks), &
F_ns(1:NV,1:nnx,1:pp,rks), &
       stat=ierr)

    F_s=0d0
    F_ns=0d0


    if (ierr .ne. 0) then
       write(*,*)'Can not allocate State Variables'
       write(*,*)'Abort!'
       stop
    endif
allocate(V(NV),VW(NV),stat=ierr)
     if (ierr .ne. 0) then
       write(*,*)'Can not allocate Variables for Discrete Ordinate Methods'
       write(*,*)'Abort!'
       stop
    endif


write(6,*) "Init_Var",NV,nnx,pp

return
!============ coeffice of ARK4_6 ============
    
  end subroutine Init_Var

  !------------------------------------------------------------------------------
end module State_VAR
