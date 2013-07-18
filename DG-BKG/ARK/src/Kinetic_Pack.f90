
subroutine BGKimexS_nsoli(xcur, fcur, iupar, rupar)
  use Legendre, ONLY: F_alt, Leg_Grid_xi1, Leg_Grid_xi2, Vx, Vy, IGH, GHW, &
        LGLWeights_Grid_xi1,LGLWeights_Grid_xi2
  use MD2D_Grid !, ONLY: Jacobian

  use Metric_Var
  use NorVec_Var
  use Kinetic_Var, ONLY: GHW, Tau, IT
  use RK_Var
  use universal_const
  implicit none
  real(kind=8), intent(in):: xcur(1:IGH*IGH*(PDeg1+1)*(PDeg2+1))
  real(kind=8), intent(out) :: fcur(1:IGH*IGH*(PDeg1+1)*(PDeg2+1))
  integer :: iupar(2)
  real(kind=8):: rupar(2)

  integer:: i,j,ii,jj,mm
  real(kind=8) :: sum_loc,gamma_ark
  real(kind=8):: vin(1:IGH,1:IGH,0:PDeg1,0:PDeg2)
  real(kind=8):: ft_a(1:IGH,1:IGH,0:PDeg1,0:PDeg2)
  real(kind=8):: SR_loc(0:PND1,0:PND2),SE_loc(0:PND1,0:PND2)
  real(kind=8):: SUx_loc(0:PND1,0:PND2),SUy_loc(0:PND1,0:PND2)
  real(kind=8):: R_loc(0:PND1,0:PND2),ET_loc(0:PND1,0:PND2)
  real(kind=8):: Ux_loc(0:PND1,0:PND2),Uy_loc(0:PND1,0:PND2)
  real(kind=8):: Z_loc(0:PND1,0:PND2),T_loc(0:PND1,0:PND2)
  real(kind=8):: F0_loc(1:IGH,1:IGH,0:PND1,0:PND2)
  real(kind=8):: F_loc(1:IGH,1:IGH,0:PND1,0:PND2)
  real(kind=8):: FS_loc(1:IGH,1:IGH,0:PND1,0:PND2)
  real(kind=8) :: FS_tmp(0:PND1,0:PND2)

  real(kind=8) :: ZA,ZB,ZC,GA1,GB1,GC1,GA2,GB2,GC2,PSIA,PSIB,PSIC
  real(kind=8) :: PP
  integer:: DDK_l

vin=reshape(xcur,(/IGH,IGH,PDeg1+1,PDeg2+1 /) )
gamma_ark=const_a_I(2,1)

      DDK_l=iupar(1)

F_loc =0d0
      do ii=1,IGH
      do jj=1,IGH

      do i=0,PDeg1
      do j=0,PDeg2
!        F_loc(ii,jj,:,:,DDK_l)=F_loc(ii,jj,:,:,DDK_l)+F_alt(ii,jj,i,j,DDK_l)* &
!          Leg_Grid_xi1(:,:,i)*Leg_Grid_xi2(:,:,j)
        F_loc(ii,jj,:,:)=F_loc(ii,jj,:,:)+vin(ii,jj,i,j)* &
          Leg_Grid_xi1(:,:,i)*Leg_Grid_xi2(:,:,j)
      end do
      end do

      end do
      end do
SR_loc = 0d0
SE_loc = 0d0
SUx_loc = 0d0
SUy_loc = 0d0

    do ii=1,IGH
    do jj=1,IGH
      SR_loc(:,:) =  SR_loc(:,:)+GHW(ii)*GHW(jj)*F_loc(ii,jj,:,:)
     SUx_loc(:,:) = SUx_loc(:,:)+GHW(ii)*GHW(jj)*F_loc(ii,jj,:,:)*Vx(ii)
     SUy_loc(:,:) = SUy_loc(:,:)+GHW(ii)*GHW(jj)*F_loc(ii,jj,:,:)*Vy(jj)
      SE_loc(:,:) =  SE_loc(:,:)+GHW(ii)*GHW(jj)*F_loc(ii,jj,:,:)* &
        (Vx(ii)**2+Vy(jj)**2)/2d0
    end do
    end do

     R_loc    = SR_loc
     Ux_loc   = SUx_loc/SR_loc
     Uy_loc   = SUy_loc/SR_loc
     ET_loc   = SE_loc
         do ii=0, PND1
         do jj=0, PND2

          T_loc(ii,jj) = (2d0* ET_loc(ii,jj)/R_loc(ii,jj))-&
            (Ux_loc(ii,jj)**2+Uy_loc(ii,jj)**2)
          Z_loc(ii,jj) = R_loc(ii,jj)/(pi*T_loc(ii,jj))
!          P_loc(ii,jj) = R_loc(ii,jj)*T_loc(ii,jj)/2d0

if ( T_loc(ii,jj) .lt. 0) then
write(6,*) "T",T_loc(ii,jj),ii,jj,DDK_l
stop
endif
         end do
         end do

    do i=1,IGH
    do j=1,IGH

     do ii=0, PND1
     do jj=0, PND2

       PP = ( (Vx(i)-Ux_loc(ii,jj))**2 + (Vy(j)-Uy_loc(ii,jj))**2 ) / T_loc(ii,jj)
       F0_loc(i,j,ii,jj) = 1d0/((EXP(PP)/Z_loc(ii,jj))+dble(IT))

    if (F0_loc(i,j,ii,jj) .gt. 1d20) then
    write(6,*) "Infty",i,j,ii,jj,DDK_l
    write(6,*) T_loc(ii,jj),Z_loc(ii,jj),R_loc(ii,jj),PP
    stop
    endif
     END DO
     END DO
  END DO
  END DO

!F_tmp=F_loc-F_0
FS_loc=0d0

!  do DDK=1,TotNum_DM
      DDK_l=iupar(1)

      do ii=1,IGH
      do jj=1,IGH

      do i=0,PDeg1
      do j=0,PDeg2

        FS_tmp=F_loc(ii,jj,:,:)-F0_loc(ii,jj,:,:)
        FS_tmp=FS_tmp*Leg_Grid_xi1(:,:,i)*Leg_Grid_xi2(:,:,j)*&
         LGLWeights_Grid_xi1(:,:)*LGLWeights_Grid_xi2(:,:)*Jacobian(:,:,DDK_l)

        sum_loc=0d0
        do mm=0,PND1
         sum_loc=sum_loc+sum(FS_tmp(:,mm))
        enddo

         FS_loc(ii,jj,i,j)=sum_loc/Tau
        end do
        end do

      end do
      end do

! enddo
         ft_a= vin - F_alt(:,:,:,:,DDK_l)-gamma_ark*dt*FS_loc

fcur=reshape(ft_a,(/ IGH*IGH*(PDeg1+1)*(PDeg2+1) /) )

return
end subroutine BGKimexS_nsoli

subroutine Compute_Source
  use Legendre
  use MD2D_Grid
  use Metric_Var
  use NorVec_Var
  use Kinetic_Var

  implicit none
  integer:: i,j,ii,jj,mm
  real(kind=8) :: FS_tmp(0:PND1,0:PND2)
  real(kind=8) :: sum_loc
  real(kind=8),allocatable :: temp(:,:)
  real(kind=8),allocatable :: temp1(:,:),tempsum1(:,:)
  real(kind=8),allocatable :: temp3(:,:),tempsum3(:,:)

!F_tmp=F_loc-F_0
FS=0d0

  do DDK=1,TotNum_DM

      do ii=1,IGH
      do jj=1,IGH

      do i=0,PDeg1
      do j=0,PDeg2

        FS_tmp=F_loc(ii,jj,:,:,DDK)-F_0(ii,jj,:,:,DDK)
        FS_tmp=FS_tmp*Leg_Grid_xi1(:,:,i)*Leg_Grid_xi2(:,:,j)*&
         LGLWeights_Grid_xi1(:,:)*LGLWeights_Grid_xi2(:,:)*Jacobian(:,:,DDK)

        sum_loc=0d0
        do mm=0,PND1
         sum_loc=sum_loc+sum(FS_tmp(:,mm))
        enddo

         FS(ii,jj,i,j,DDK)=sum_loc/Tau
        end do
        end do

      end do
      end do

 enddo
!write(6,*) "F_loc",F_loc(1:3,1,0:2,0:2,3)
!write(6,*) "F_0",F_0(1:3,1,0:2,0:2,3)
!write(6,*) "FS",FS(1:3,1,0:2,0:2,3)

end subroutine Compute_Source


SUBROUTINE EQUILIBRIUM !(NX,NY,NV,C,V,F,R,UX,UY,ET)
use Legendre
use MD2D_Grid
use Kinetic_Var
use universal_const

IMPLICIT NONE
INTEGER i,j,ii,jj
real(kind=8) :: PP

   do DDK=1,TotNum_DM

    do i=1,IGH
    do j=1,IGH

     do ii=0, PND1
     do jj=0, PND2

                PP = ( (Vx(i)-Ux(ii,jj,DDK))**2 + (Vy(j)-Uy(ii,jj,DDK))**2 ) / T(ii,jj,DDK)
                F_0(i,j,ii,jj,DDK) = 1d0/((EXP(PP)/Z(ii,jj,DDK))+dble(IT))
    if (F_0(i,j,ii,jj,DDK) .gt. 1d20) then
    write(6,*) "Infty",i,j,ii,jj,DDK
    write(6,*) T(ii,jj,DDK),Z(ii,jj,DDK),R_loc(ii,jj,DDK),PP
    stop
    endif
     END DO
     END DO
  END DO
  END DO
END DO


END SUBROUTINE EQUILIBRIUM

SUBROUTINE CALDOM !(NX,NY,NV,C,V,F,R,UX,UY,ET)
use Legendre
use MD2D_Grid
use Kinetic_Var
use universal_const

IMPLICIT NONE

INTEGER i,j,ii,jj,L
logical isnan

!REAL(8) SR,SUx,SUy,SE
!REAL, DIMENSION (IGH) :: C,V
!REAL, DIMENSION (IMAXGRIDX,IMAXGRIDY) :: R,UX,UY,ET
!REAL, DIMENSION (IGH,IGH,IMAXGRIDX,IMAXGRIDY) :: F

 real(kind=8) :: ZA,ZB,ZC,GA1,GB1,GC1,GA2,GB2,GC2,PSIA,PSIB,PSIC

! Compute SU, SE, SAV, SR

SR=0d0
SUx=0d0
SUy=0d0
SE=0d0

F_loc =0d0

    do DDK=1,TotNum_DM

      do ii=1,IGH
      do jj=1,IGH

      do i=0,PDeg1
      do j=0,PDeg2
        F_loc(ii,jj,:,:,DDK)=F_loc(ii,jj,:,:,DDK)+F_alt(ii,jj,i,j,DDK)*Leg_Grid_xi1(:,:,i)*Leg_Grid_xi2(:,:,j)
      end do
      end do

      end do
      end do

    do ii=1,IGH
    do jj=1,IGH
      SR(:,:,DDK) =  SR(:,:,DDK)+GHW(ii)*GHW(jj)*F_loc(ii,jj,:,:,DDK)
     SUx(:,:,DDK) = SUx(:,:,DDK)+GHW(ii)*GHW(jj)*F_loc(ii,jj,:,:,DDK)*Vx(ii)
     SUy(:,:,DDK) = SUy(:,:,DDK)+GHW(ii)*GHW(jj)*F_loc(ii,jj,:,:,DDK)*Vy(jj)
      SE(:,:,DDK) =  SE(:,:,DDK)+GHW(ii)*GHW(jj)*F_loc(ii,jj,:,:,DDK)*(Vx(ii)**2+Vy(jj)**2)/2d0
    end do
    end do

   end do
 !!!!!!!!!!!!!!!!!!!!!!!!!!

          R_loc    = SR
          Ux   = SUx/SR
          Uy   = SUy/SR
          ET   = SE
!write(6,*) "F_0",F_0(1:3,1,0:2,0:2,1) 
!write(6,*) "F_alt",F_alt(1:3,1,0:2,0:2,1)
!write(6,*) "F_loc",F_loc(1:3,1,0:2,0:2,1)

if (isnan(SR(2,2,1))) then
write(6,*) "NaN SR"
!write(6,*) "F_loc",F_loc(:,:,:,:,1)
stop
endif


  select case(IT)

    case (0) ! MAXWELLIAN
        do DDK=1,TotNum_DM
         do ii=0, PND1
         do jj=0, PND2

          T(ii,jj,DDK) = (2d0* ET(ii,jj,DDK)/R_loc(ii,jj,DDK))-(Ux(ii,jj,DDK)**2+Uy(ii,jj,DDK)**2)
          Z(ii,jj,DDK) = R_loc(ii,jj,DDK)/(pi*T(ii,jj,DDK))
          P(ii,jj,DDK) = R_loc(ii,jj,DDK)*T(ii,jj,DDK)/2d0

if ( T(ii,jj,DDK) .lt. 0) then
write(6,*) "T",T(ii,jj,DDK),ii,jj,DDK
write(6,*) R_loc(ii,jj,DDK),Ux(ii,jj,DDK),Uy(ii,jj,DDK),ET(ii,jj,DDK)
write(6,*) "F_loc",F_loc(:,:,ii,jj,DDK)

stop
endif
         end do
         end do
        end do

       case(-1)
! maxwellian = 0., fermion = 1., boson = -1
        do DDK=1,TotNum_DM
         do ii=0, PND1
         do jj=0, PND2

          ZA = 0.0001d0
          ZB = 0.99d0
  DO WHILE (ABS(ZA-ZB) .GT. 1D-13)
    GA1 = 0d0
    GB1 = 0d0
    GA2 = 0d0
    GB2 = 0d0
       DO L = 1, 100
            GA1 = GA1 + (ZA**L) /dble(L)
            GB1 = GB1 + (ZB**L) /dble(L)
            GA2 = GA2 + (ZA**L) /dble(L**2)
            GB2 = GB2 + (ZB**L) /dble(L**2)
        END DO
    PSIA = 2d0*ET(ii,jj,DDK) - (GA2*(R_loc(ii,jj,DDK)/GA1)**2)/pi - &
 R_loc(ii,jj,DDK)*(Ux(ii,jj,DDK)**2+Uy(ii,jj,DDK)**2)
    PSIB = 2d0*ET(ii,jj,DDK) - (GB2*(R_loc(ii,jj,DDK)/GB1)**2)/pi - &
 R_loc(ii,jj,DDK)*(Ux(ii,jj,DDK)**2+Uy(ii,jj,DDK)**2)
    ZC = (ZA+ZB)/2d0
    GC1 = 0d0
    GC2 = 0d0
        DO L = 1, 100
            GC1 = GC1 + (ZC**L)/dble(L)
            GC2 = GC2 + (ZC**L)/dble(L**2)
        END DO
    PSIC = 2d0*ET(ii,jj,DDK) - (GC2*(R_loc(ii,jj,DDK)/GC1)**2)/pi - &
R_loc(ii,jj,DDK)*(Ux(ii,jj,DDK)**2+Uy(ii,jj,DDK)**2)

    IF ((PSIA*PSIC) .LT. 0) THEN
        ZB = ZC
    ELSE
        ZA = ZC
    END IF
  END DO
        Z(ii,jj,DDK) = ZC;
!        T(ii,jj,DDK) = R(ii,jj,DDK)**2 / (pi*GC1);
        T(ii,jj,DDK) = R_loc(ii,jj,DDK) / (pi*GC1);
        P(ii,jj,DDK) = ET(ii,jj,DDK) - 0.5d0* R_loc(ii,jj,DDK) *( Ux(ii,jj,DDK)**2+Uy(ii,jj,DDK)**2)
         end do
     end do
  enddo
        case(1)
!           %(IT == 1)
!maxwellian = 0., fermion = 1., boson = -1
            do DDK=1,TotNum_DM
             do ii=0, PND1
             do jj=0, PND2

          ZA = 0.0001d0
          ZB = 0.99d0
  DO WHILE (ABS(ZA-ZB) .GT. 1D-13)
    GA1 = 0d0
    GB1 = 0d0
    GA2 = 0d0
    GB2 = 0d0
       DO L = 1, 100
            GA1 = GA1 + (ZA**L) * (-1)**(L-1)/dble(L)
            GB1 = GB1 + (ZB**L) * (-1)**(L-1)/dble(L)
            GA2 = GA2 + (ZA**L) * (-1)**(L-1)/dble(L**2)
            GB2 = GB2 + (ZB**L) * (-1)**(L-1)/dble(L**2)
        END DO
    PSIA = 2d0*ET(ii,jj,DDK) - (GA2*(R_loc(ii,jj,DDK)/GA1)**2)/pi - &
 R_loc(ii,jj,DDK)*(Ux(ii,jj,DDK)**2+Uy(ii,jj,DDK)**2)
    PSIB = 2d0*ET(ii,jj,DDK) - (GB2*(R_loc(ii,jj,DDK)/GB1)**2)/pi - &
 R_loc(ii,jj,DDK)*(Ux(ii,jj,DDK)**2+Uy(ii,jj,DDK)**2)
    ZC = (ZA+ZB)/2d0
    GC1 = 0d0
    GC2 = 0d0
        DO L = 1, 100
            GC1 = GC1 + (ZC**L) * (-1)**(L-1)/dble(L)
            GC2 = GC2 + (ZC**L) * (-1)**(L-1)/dble(L**2)
        END DO
    PSIC = 2d0*ET(ii,jj,DDK) - (GC2*(R_loc(ii,jj,DDK)/GC1)**2)/pi - &
R_loc(ii,jj,DDK)*(Ux(ii,jj,DDK)**2+Uy(ii,jj,DDK)**2)

    IF ((PSIA*PSIC) .LT. 0) THEN
        ZB = ZC
    ELSE
        ZA = ZC
    END IF
  END DO
        Z(ii,jj,DDK) = ZC;
!        T(ii,jj,DDK) = R(ii,jj,DDK)**2 / (pi*GC1);
        T(ii,jj,DDK) = R_loc(ii,jj,DDK) / (pi*GC1);

        P(ii,jj,DDK) = ET(ii,jj,DDK) - 0.5d0* R_loc(ii,jj,DDK) *( Ux(ii,jj,DDK)**2+Uy(ii,jj,DDK)**2)
         end do
     end do
 enddo

            end select !if IT

!  write(6,*) "IT",IT
!  write(6,*) "Ux",Ux(0:2,0:2,1) !:TotNum_DM)
!  write(6,*) "Uy",Uy(0:2,0:2,1) ! :TotNum_DM)
!  write(6,*) "Z",Z(0:2,0:2,1) !:TotNum_DM)
!  write(6,*) "T",T(0:2,0:2,1) !:TotNum_DM)
!  write(6,*) "SR",SR(0:2,0:2,1) !:TotNum_DM)
!  write(6,*) "R",R_loc(0:2,0:2,1) !:TotNum_DM)
!  write(6,*) "ET",ET(0:2,0:2,1) !:TotNum_DM)


if (isnan(Z(2,2,1))) then
write(6,*) "NaN Z"
write(6,*) "R",R_loc(0:2,0:2,1)
 write(6,*) "Ux",Ux(0:2,0:2,1)
 write(6,*) "Uy",Uy(0:2,0:2,1)
 write(6,*) "Z",Z(0:2,0:2,1)
 write(6,*) "T",T(0:2,0:2,1)
write(6,*) "ET",ET(0:2,0:2,1)
stop
endif

end SUBROUTINE CALDOM

subroutine DG2D_Kinetic_Initial(Time_final,dt_eq)
  use Metric_Var
  use NorVec_Var
  use MD2D_Grid
  use Legendre
  use universal_const
  use State_Var
  use Kinetic_Var

  implicit none

  integer:: i,j,k,ierr,ii,jj
!  real(kind=8) :: a_vec(1:2),
  real(kind=8) ::Time_final
  real(kind=8) :: UX1,UY1,Z1,Ti1
  real(kind=8) :: UX2,UY2,Z2,Ti2
  real(kind=8) :: UX3,UY3,Z3,Ti3
  real(kind=8) :: UX4,UY4,Z4,Ti4
  real(kind=8) :: PP, MF,eigMF,dt_eq,Cx,Cy

  real(kind=8) :: wp(20)
  real(kind=8),allocatable :: b_ini(:,:),Pi_Pm(:,:),Pi_Pj(:,:)
  real(kind=8),allocatable :: temp(:,:),tempsum(:),func(:,:,:)
  real(kind=8),allocatable :: LGLWeights_2D(:,:)
  real(kind=8),allocatable :: ValOfPolyNatGrids(:) ! ,cs_1d(:),w_1d(:)
  allocate( b_ini(0:PDeg1,0:PDeg2),Pi_Pm(0:PDeg1,0:PDeg2), &
            Pi_Pj(0:PND1,0:PND2),func(0:PND1,0:PND2,1:TotNum_DM), &
            temp(0:PND1,0:PND2),tempsum(0:PND1), &
            LGLWeights_2D(0:PND1,0:PND1) ,stat=ierr )

  LGLWeights_2D = LGLWeights_Grid_xi1*LGLWeights_Grid_xi2

  if (ierr .ne. 0) then
     write(*,*)"Can't allocate u"
     stop
  endif

if (I_prob .ge. 1) then
allocate( ValOfPolyNatGrids(1:IGH))
write(6,*) IGH

GH=0d0

call ZEHEGA(IGH,GH,ValOfPolyNatGrids)
call WEHEGA(IGH,GH,GHW)
write(6,*) "GH",GH
!write(6,*) GHW
GHW=GHW*exp(GH**2)
write(6,*) "GHW",GHW

wp  =(/0.898591961453,0.704332961176,0.62227869619,0.575262442852, &
    0.544851742366,0.524080350949,0.509679027117,0.499920871336, &
    0.493843385272,0.490921500667,0.490921500667,0.493843385272, & 
    0.499920871336,0.509679027117,0.524080350949,0.544851742366, &
    0.575262442852,0.62227869619,0.704332961176,0.898591961453/)

write(6,*) "wp",wp

Vx=GH
Vy=GH

!GH(1)=a_vec(1)
!GH(2)=a_vec(2)
!GH(2:IGH)=a_vec(2)

write(6,*) GH
deallocate( ValOfPolyNatGrids)

endif
select case(I_Case)
  case(5)
NumCase='05'
    Z1   = 0.142d0
    UX1  = -0.75d0
    UY1  = -0.5d0
    Ti1   = 2.078d0

    Z2   = 0.4253d0
    UX2  = -0.75d0
    UY2  = 0.5d0
    Ti2   = 1.1494d0

    Z3   = 0.142d0
    UX3  = 0.75d0
    UY3  = 0.5d0
    Ti3   = 2.078d0

    Z4   = 0.6635d0
    UX4  = 0.75d0
    UY4  = -0.5d0
    Ti4   = 0.87685d0
  case(17)! init17
NumCase='17'

Z2=0.37d0
UY2=-0.3d0
UX2=0d0
Ti2=1.25d0

Z1=0.14d0
UY1=-0.4d0 
UX1=0.1d0
Ti1=2.08d0

Z3=0.38d0
UY3=0.12d0
UX3=0d0
Ti3=0.77d0

Z4=0.14d0
UY4=-0.98d0 
UX4=0d0
Ti4=1.31d0

case(13)! init13
NumCase='13'
Z2=0.4253d0
UY2=0.3d0
UX2=0d0
Ti2=1.1494d0

Z1=0.142d0
UY1=-0.3d0
UX1=0.d0
Ti1=2.0782d0

Z3=0.4448d0
UY3=0.697d0
UX3=0d0
Ti3=0.7083d0

Z4=0.151d0
UY4=0.26254d0
UX4=0d0
Ti4=1.273d0
case default
write(6,*) "NO SUCH Case",I_Case
stop
end select
!setup initial condition
if (DM_Range .gt. 0) then
Cx=0.5d0
Cy=0.5d0
else
Cx=0d0
Cy=0d0
endif

if (I_prob .lt. 1) then
       Z=Z1
       T=Ti1
       Ux=UX1
       Uy=UY1


elseif (I_prob .le. 3) then


   do DDK=1,TotNum_DM

!     do ii=0, PND1
!     do jj=0, PND2
     if (x1(2,2,DDK) .le. Cx) then
      if (x2(2,2,DDK) .le. Cy) then
       Z(0:PND1,0:PND2,DDK)=Z3
       T(0:PND1,0:PND2,DDK)=Ti3
       Ux(0:PND1,0:PND2,DDK)=UX3
       Uy(0:PND1,0:PND2,DDK)=UY3
      else
       Z(0:PND1,0:PND2,DDK)=Z2
       T(0:PND1,0:PND2,DDK)=Ti2
       Ux(0:PND1,0:PND2,DDK)=UX2
       Uy(0:PND1,0:PND2,DDK)=UY2
      endif
     else 
      if (x2(2,2,DDK) .le. Cy) then
       Z(0:PND1,0:PND2,DDK)=Z4
       T(0:PND1,0:PND2,DDK)=Ti4
       Ux(0:PND1,0:PND2,DDK)=UX4
       Uy(0:PND1,0:PND2,DDK)=UY4
      else
       Z(0:PND1,0:PND2,DDK)=Z1
       T(0:PND1,0:PND2,DDK)=Ti1
       Ux(0:PND1,0:PND2,DDK)=UX1
       Uy(0:PND1,0:PND2,DDK)=UY1
      endif
     endif


!     END DO
!     END DO
END DO

elseif (I_prob .gt. 4) then

  do DDK = 1,TotNum_DM
     do i = 0,PND1
        do j = 0,PND2
           func(i,j,DDK) = dcos(2*pi*(x1(i,j,DDK)+x2(i,j,DDK)))
        enddo
     enddo
  enddo
endif
 
call EQUILIBRIUM


!   do i,j=1,IGH
!     do ii,jj=0, PND1
!  F_0(i,j,ii,jj,DDK)

  do i = 0,PDeg1
     do j = 0,PDeg2
        Pi_Pm(i,j) = sum(Leg_Grid_xi1(:,0,i)*Leg_Grid_xi1(:,0,j)*LGLWeights_Grid_xi1(:,0))
     enddo
  enddo

  do DDK = 1,TotNum_DM
!     B(1,DDK) = a_vec(1)*dx2_dxi2(0,0,DDK)-a_vec(2)*dx1_dxi2(0,0,DDK)
!     B(2,DDK) = -a_vec(1)*dx2_dxi1(0,0,DDK)+a_vec(2)*dx1_dxi1(0,0,DDK)

     Bq(2,2,DDK) = dx2_dxi2(0,0,DDK)
     Bq(2,1,DDK) = dx2_dxi1(0,0,DDK)
     Bq(1,2,DDK) = dx1_dxi2(0,0,DDK)
     Bq(1,1,DDK) = dx1_dxi1(0,0,DDK)
!     B(1,DDK) = a_vec(1)*Bq(2,2,DDK)-a_vec(2)*Bq(1,2,DDK)
!     B(2,DDK) = -a_vec(1)*Bq(2,1,DDK)+a_vec(2)*Bq(1,1,DDK)

  do i = 1,IGH
     do j = 1,IGH
     B_gen(1,i,j) = GH(i)*Bq(2,2,DDK)-GH(j)*Bq(1,2,DDK)
     B_gen(2,i,j) = -GH(i)*Bq(2,1,DDK)+GH(j)*Bq(1,1,DDK)
     enddo
  enddo
do ii=1,IGH
do jj=1,IGH
     b_ini = 0.0d0

     do i = 0,PDeg1
     do j = 0,PDeg2

           Pi_Pj = Leg_Grid_xi1(:,:,i)*Leg_Grid_xi2(:,:,j)
           temp(:,:) = F_0(ii,jj,:,:,DDK)*Pi_Pj*Jacobian(:,:,DDK)*LGLWeights_2D

           do k = 0,PND1
              tempsum(k) = sum(temp(k,:))
           enddo
           b_ini(i,j) = sum(tempsum)
           A1(i,j,DDK) = Pi_Pm(i,i)*Pi_Pm(j,j)*Jacobian(0,0,DDK)
           A1(i,j,DDK) = 1d0/A1(i,j,DDK)
        enddo
     enddo
     F_new(ii,jj,:,:,DDK) = b_ini*A1(:,:,DDK)
  enddo
enddo
enddo
  do i = 0,PDeg1
     do j = 0,PDeg2
        tempsum = Leg_Grid_xi1(:,0,i)*DLeg_Grid_xi1(:,0,j)*LGLWeights_Grid_xi1(:,0)
        B_tal1x(i,j) = sum(tempsum)

        tempsum = Leg_Grid_xi2(0,:,j)*Leg_Grid_xi2(0,:,j)*LGLWeights_Grid_xi2(0,:)
        B_tal1y(i,j) = sum(tempsum)

        tempsum = Leg_Grid_xi1(:,0,i)*Leg_Grid_xi1(:,0,i)*LGLWeights_Grid_xi1(:,0)
        B_tal2x(i,j) = sum(tempsum)

     enddo
  enddo
!  write(6,*) "F_new",F_new(12,10,0:PDeg1,0:PDeg2,1)
!  write(6,*) "F_new",F_new(12,11,0:PDeg1,0:PDeg2,1)

  F_alt=F_new
  deallocate( b_ini,Pi_Pj,Pi_Pm,temp,tempsum,func,LGLWeights_2D )
  call CALDOM
  write(6,*) "Theta",IT,"I_prob",I_prob
  write(6,*) "Ux1",UX1,Ux(1:2,1:2,1:4)
  write(6,*) "Uy1",UY1,Uy(1:2,1:2,1:4)
  write(6,*) "Z1",Z1,Z(1:2,1:2,1:4)
  write(6,*) "T1",Ti1,T(1:2,1:2,1:4)
  write(6,*) "R_loc",R_loc(1:2,1:2,1:4)
ii=1
jj=1

do DDK = 1,TotNum_DM
    do i=1,IGH
    do j=1,IGH

PP = ( (Vx(i)-Ux(ii,jj,DDK))**2 + (Vy(j)-Uy(ii,jj,DDK))**2 ) / T(ii,jj,DDK)

MF=(exp( PP )*GHW(i)*GHW(j)* (2d0*R_loc(ii,jj,DDK)*T(ii,jj,DDK) + &
 pi* ( PP*T(ii,jj,DDK) )**2 *Z(ii,jj,DDK) ))/&
( pi*R_loc(ii,jj,DDK)*T(ii,jj,DDK)**2*( exp( PP ) + IT*Z(ii,jj,DDK) )**2 )
eigMF=max(abs(MF),eigMF)
enddo
enddo
enddo
write(6,*) "Tau",Tau
dt_eq=Tau/(1d0+eigMF*IGH**2*(PDeg1+1)**2)

write(6,*) "End of DG2D Initial"
end subroutine DG2D_Kinetic_Initial

logical function isnan(x)
real x
isnan = x .ne. x
end
