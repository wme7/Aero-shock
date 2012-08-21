!===============================================================!
!                                                               !
! This is a module defining variables for 2D Elastic Waves      !
!                                                               !
!===============================================================!
module Kinetic_Var  ! Kinetic
implicit none
!--------------<Cartesian Coordinate>---------------------
! PDeg : dof
! PND  : quad points
!
real(kind=8), save, allocatable::  SR(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)
real(kind=8), save, allocatable::  SUx(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)
real(kind=8), save, allocatable::  SUy(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)
real(kind=8), save, allocatable::  SE(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)

real(kind=8), save, allocatable::  R_loc(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)
real(kind=8), save, allocatable::  Ux(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)
real(kind=8), save, allocatable::  Uy(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)
real(kind=8), save, allocatable::  ET(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)

real(kind=8), save, allocatable::  P(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)
real(kind=8), save, allocatable::  T(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)
real(kind=8), save, allocatable::  Z(:,:,:) ! (0:PND1, 0:PND2,1:TotNum_DM)

!real(kind=8), save, allocatable::  AV(:,:,:) ! (0:Pdeg1, 0:PDeg2,1:TotNum_DM)

real(kind=8), save:: Tau = 0.01d0 !RELAXATION TIME

!real(kind=8), save, allocatable:: GH(:,:),GHW(:,:)

integer:: IT

contains

!------------------------------------------------------------------------------

  subroutine Init_Kinetic_Variables(Deg_Max,TotNumDomain)
    implicit none
    ! declare subroutine arguments
    integer:: Deg_Max(2)                ! Pass in the Nxi_max
    integer:: TotNumDomain
   

    integer:: ND1,ND2 
    integer:: ierr
   
!    IT=0
    
    ND1=Deg_Max(1); ND2=Deg_Max(2);

    allocate(  SR(0:ND1,0:ND2,1:TotNumDomain), &
               SUx(0:ND1,0:ND2,1:TotNumDomain), &
               SUy(0:ND1,0:ND2,1:TotNumDomain), &
               SE(0:ND1,0:ND2,1:TotNumDomain), &
               R_loc(0:ND1,0:ND2,1:TotNumDomain), &
               Ux(0:ND1,0:ND2,1:TotNumDomain), &
               Uy(0:ND1,0:ND2,1:TotNumDomain), &
               ET(0:ND1,0:ND2,1:TotNumDomain), &
               P(0:ND1,0:ND2,1:TotNumDomain), &
               T(0:ND1,0:ND2,1:TotNumDomain), &
               Z(0:ND1,0:ND2,1:TotNumDomain), &
              stat=ierr)

    if ( ierr .ne. 0 ) then
       write(*,*)'Cannot allocate memory for state variable q'
       write(*,*)'Abort!'
       stop
    endif

    return 
    
  end subroutine Init_Kinetic_Variables

!-------------------------------------------------------------------------------  
  
end module Kinetic_Var
