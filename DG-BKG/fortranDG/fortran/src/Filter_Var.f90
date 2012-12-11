module Filter_Var
  implicit none
  real(kind=8), save :: cutoff
  integer, save :: filter_order,filter_switch
  real(kind=8), save :: filter_alpha
  real(kind=8), save, allocatable :: filter_sigma(:)
  real(kind=8), save, allocatable :: ValPolyAtGrid(:,:)
  real(kind=8), save, allocatable :: filter_oper_xi1(:,:,:)
  real(kind=8), save, allocatable :: filter_oper_xi2(:,:,:)
  real(kind=8), save, allocatable :: Q(:,:), Qt(:,:)
  real(kind=8), allocatable :: PnAtgrid(:),dPnAtgrid(:),d2PnAtgrid(:)
  real(kind=8), allocatable :: gamma_N(:)
  integer :: DegCutOff

contains

  !=================================================================!

  subroutine Read_In_Filter_Profile
    implicit none
    integer Filter_OnOff
    
    integer :: lid
    
    lid = 21
    open(lid,file='filter.in', form='formatted', status='old')
    read(lid,*) 
    read(lid,*) filter_switch !on=1, off=0
    read(lid,*)
    read(lid,*)filter_order
    read(lid,*)cutoff		! number between [0 1 ]
    read(lid,*)
    close(lid)
    
    return
    
  end subroutine Read_In_Filter_Profile

  !=================================================================!  

  subroutine filter_profile(PolyDeg)
    implicit none
    integer :: PolyDeg
    !
    !declare local arguments
    integer :: i,j
    real(kind=8) :: machine_zero, tmp
    
    machine_zero=1.0E-14
    filter_alpha=-log(machine_zero)
    DegCutOff=int(CutOff*dble(PolyDeg))
       write(*,*)filter_order,CutOff,PolyDeg!; pause
    do i=0,PolyDeg
       
       if ( i .gt. DegCutOff ) then
          tmp=(dble(i-DegCutOff)/dble(PolyDeg-DegCutOff))**filter_order
          filter_sigma(i)=dexp(-filter_alpha*tmp)
       else
          filter_sigma(i)=1.d0
       endif
       
       !	  filter_sigma(i)=1.d0
    enddo
    
    do i=0, PolyDeg
       if ( filter_sigma(i) .le. 1.0E-14 )  filter_sigma(i)=0.d0
    enddo
    
!!$       do i=0,PolyDeg
          write(*,*) "Filter",filter_alpha,filter_sigma
!!$       enddo
!!$    
!!$       pause
    return
    
  end subroutine filter_profile

  !=================================================================!
  
  subroutine alloc_mem_filter_variables(DegMax)
    implicit none
    integer :: DegMax
    !
    integer :: ierr
    !
    
    allocate ( filter_sigma(0:DegMax), &
         ValPolyAtGrid(0:DegMax,0:DegMax), &
         filter_oper_xi1(0:DegMax,0:DegMax,0:DegMax), &
         filter_oper_xi2(0:DegMax,0:DegMax,0:DegMax), stat=ierr )
    
    if ( ierr .ne. 0 ) then
       write(*,*)'Can not allocate memory for filter variables'
       write(*,*)'Abort!'
       stop
    endif
    
    filter_sigma=0.d0
    ValPolyAtGrid=0.d0
    filter_oper_xi1=0.d0
    filter_oper_xi2=0.d0
    
    allocate (  Q(0:DegMax,0:DegMax), &
         Qt(0:DegMax,0:DegMax), stat=ierr )
    if (ierr .ne. 0) then
       write(*,*) 'Can not allocate memory for Qt'
       write(*,*) 'abort!'
       stop
    endif
    Q=0.d0; Qt=0.d0

    allocate(  PnAtgrid(0:DegMax), &
              dPnAtgrid(0:DegMax), &
             d2PnAtgrid(0:DegMax), &
                gamma_N(0:DegMax), stat=ierr)
  
  if (ierr .ne. 0) then
     write(*,*) 'Can not allocate memory for ValPolyDegNAtgrid'
     write(*,*) 'abort!'
     stop
  endif
  
  PnAtgrid=0.d0; dPnAtgrid=0.d0; d2PnAtgrid=0.d0
  gamma_N=0.d0;
    
    
    return
    
  end subroutine alloc_mem_filter_variables
  
end module Filter_Var
