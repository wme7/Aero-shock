 subroutine Data_Output(dtt,nx,pp)
  use State_Var
  use Legendre
  use RK_Var
  use universal_const
!  character(len=2) :: NumCase
  implicit none
  integer      :: i,j,ii,jj,L 
  integer      :: lid
  integer      :: nx,pp
  integer      :: p
  integer:: ierr

 real(kind=8):: dtt
  
 character(len=20) :: file_tec
 character(len=5) :: NQ,NumDM
! real(kind=8) :: xo
! real(kind=8),allocatable :: SR_o(:,:),SU_o(:,:),SUy_o(:,:),SE_o(:,:)
! real(kind=8),allocatable :: R_o(:,:),Ux_o(:,:),Uy_o(:,:),ET_o(:,:)
! real(kind=8),allocatable :: P_o(:,:),T_o(:,:),Z_o(:,:),F_out(:,:,:,:)
! real(kind=8),allocatable :: Leg_Grid_xi1_o(:,:,:),Leg_Grid_xi2_o(:,:,:)
! real(kind=8),allocatable :: leg_tb_o(:,:),PN(:),PD(:)
! real(kind=8),allocatable :: LG_grids_o(:), ValuesOfPolyNatGrids_o(:)
! real(kind=8),allocatable :: xo1(:)

 p=pp-1
    
write(NQ,1000) p
write(NumDM,2000) nx

1000  format(I2)
2000  format(I4)

  file_tec =  './N'//trim(adjustl(NumDM))//'Q'//trim(adjustl(NQ))//'.dat'

  lid = 200
write(6,*) file_tec
  open(lid,file=file_tec)

!        write(lid,*) '"CFL =',CCFL,' dt=',dtt,'"'
!        write(lid,*) '"METHOD=','DG"'
!        write(lid,*) '"TotNum_DM =',TotNum_DM,'"'
!        write(lid,*) '"STATISTICS (0=Maxwellian, 1=Fermion, -1=Boson) =',IT,'"'

!            write(lid,*) '"NAVIER-STOKES, ','Tau =',Tau,'"'
!            write(lid,*) '"Configuration = ',I_case,'"'
!            write(lid,*) '"Polynomial Degree = ',PDeg1,'"'


!    write(lid,888)
! 888 format (1X,'VARIABLES = "x","y","n","p","z"')
do i=1,nx
do j=1,pp
   write (lid,*) x(i,j), Rs(i,j),Us(i,j) !,Ts(i,j)
end do
enddo

  close(lid)




! 1000 format(2e23.15)

 end subroutine Data_Output
