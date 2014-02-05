 subroutine Tecplt_Output(dtt)
  use State_Var
  use Kinetic_Var
  use MD2D_Grid
  use Legendre
  use RK_Var
  use universal_const
!  character(len=2) :: NumCase
  implicit none
  integer      :: i,j,ii,jj,L 
  integer      :: lid
  integer      :: NDo1,NDo2
    integer:: ierr

  real(kind=8):: dtt
  
  character(len=20) :: file_tec
character(len=5) :: NQ
 real(kind=8) :: ZA,ZB,ZC,GA1,GB1,GC1,GA2,GB2,GC2,PSIA,PSIB,PSIC,xo
 real(kind=8),allocatable :: SR_o(:,:),SUx_o(:,:),SUy_o(:,:),SE_o(:,:)
 real(kind=8),allocatable :: R_o(:,:),Ux_o(:,:),Uy_o(:,:),ET_o(:,:)
 real(kind=8),allocatable :: P_o(:,:),T_o(:,:),Z_o(:,:),F_out(:,:,:,:)
 real(kind=8),allocatable :: Leg_Grid_xi1_o(:,:,:),Leg_Grid_xi2_o(:,:,:)
 real(kind=8),allocatable :: leg_tb_o(:,:),PN(:),PD(:)
 real(kind=8),allocatable :: LG_grids_o(:), ValuesOfPolyNatGrids_o(:)
 real(kind=8),allocatable :: xo1(:),xo2(:)


    NDo1=nint(400/sqrt(dble(TotNum_DM)))
    NDo1=max(NDo1,PND1)
    NDo2=NDo1
    
!NQ=adjustl(PDeg1)
!write(6,*) NQ
write(NQ,1000) PDeg1
1000  format(I1)
select case(IT)
case(0)
  file_tec =  './'//NumCase//'M'//NumDM//'Q'//trim(NQ)//'.plt'
case(1)
  file_tec =  './'//NumCase//'F'//NumDM//'Q'//trim(NQ)//'.plt'
case(-1)
  file_tec =  './'//NumCase//'B'//NumDM//'Q'//trim(NQ)//'.plt'
end select
  lid = 200
write(6,*) file_tec
  open(lid,file=file_tec)

        write(lid,*) '"CFL =',CCFL,' dt=',dtt,'"'
        write(lid,*) '"METHOD=','DG"'
        write(lid,*) '"TotNum_DM =',TotNum_DM,'"'
        write(lid,*) '"STATISTICS (0=Maxwellian, 1=Fermion, -1=Boson) =',IT,'"'
!          if (isolver .eq. 1) then

            write(lid,*) '"NAVIER-STOKES, ','Tau =',Tau,'"'
            write(lid,*) '"Configuration = ',I_case,'"'
            write(lid,*) '"Polynomial Degree = ',PDeg1,'"'

!          else
!            write(lid,*) '"EULER'
!          end if

    write(lid,888)
 888 format (1X,'VARIABLES = "x","y","n","p","z"')

if (NDo1 .gt. PND1) then

   do DDK=1,TotNum_DM

    write(lid,*) 'ZONE T="Dm',DDK,'",i=',PND2+1,',j=',PND2+1

         do ii=0, PND1
         do jj=0, PND2

   write (lid,*) x1(ii,jj,DDK),x2(ii,jj,DDK), R_loc(ii,jj,DDK), T(ii,jj,DDK), Z(ii,jj,DDK)

end do
end do
end do

else
  allocate( Leg_Grid_xi1_o(0:NDo1,0:NDo2,0:PDeg1),Leg_Grid_xi2_o(0:NDo1,0:NDo2,0:PDeg1), stat=ierr)
  allocate( leg_tb_o(0:NDo1,0:PDeg1),PN(0:PDeg1),PD(0:PDeg1))
  allocate( LG_grids_o(0:NDo1), ValuesOfPolyNatGrids_o(0:NDo1))
  allocate( xo1(0:NDo1), xo2(0:NDo1))

      call ZELEGL(NDo1, LG_grids_o, ValuesOfPolyNatGrids_o)
write(6,*) "LG_grids_o",LG_grids_o
      do i=0,NDo1
       xo=LG_grids_o(i)
       call LPN(PDeg1,xo,PN,PD)
       leg_tb_o(i,0:PDeg1) = PN(0:PDeg1)
    enddo



  do j = 0,PDeg1
     do i = 0,NDo1
        Leg_Grid_xi1_o(i,0:NDo2,j) = leg_tb_o(i,j)
     enddo
     Leg_Grid_xi2_o(0:NDo1,0:NDo1,j) = transpose(Leg_Grid_xi1_o(0:NDo1,0:NDo1,j))
  enddo




    allocate(  SR_o(0:NDo1,0:NDo2), &
               SUx_o(0:NDo1,0:NDo2), &
               SUy_o(0:NDo1,0:NDo2), &
               SE_o(0:NDo1,0:NDo2), &
               R_o(0:NDo1,0:NDo2), &
               Ux_o(0:NDo1,0:NDo2), &
               Uy_o(0:NDo1,0:NDo2), &
               ET_o(0:NDo1,0:NDo2), &
               P_o(0:NDo1,0:NDo2), &
               T_o(0:NDo1,0:NDo2), &
               Z_o(0:NDo1,0:NDo2), &
               F_out(1:IGH,1:IGH,0:NDo1,0:NDo2), &
              stat=ierr)




! Compute SU, SE, SAV, SR

    do DDK=1,TotNum_DM
 do ii=0, NDo1
 xo1(ii)= x1(0,1,DDK)+(x1(PND1,1,DDK)-x1(0,1,DDK))*(LG_grids_o(ii)-LG_grids_o(0))/2d0
enddo
 do jj=0, NDo1
 xo2(jj)= x2(1,0,DDK)+(x2(1,PND1,DDK)-x2(1,0,DDK))*(LG_grids_o(jj)-LG_grids_o(0))/2d0

enddo
F_out =0d0
SR_o=0d0
SUx_o=0d0
SUy_o=0d0
SE_o=0d0


      do ii=1,IGH
      do jj=1,IGH

      do i=0,PDeg1
      do j=0,PDeg2
        F_out(ii,jj,:,:)=F_out(ii,jj,:,:)+F_alt(ii,jj,i,j,DDK)*Leg_Grid_xi1_o(:,:,i)*&
Leg_Grid_xi2_o(:,:,j)
      end do
      end do

      end do
      end do

    do ii=1,IGH
    do jj=1,IGH
      SR_o(:,:) =  SR_o(:,:)+GHW(ii)*GHW(jj)*F_out(ii,jj,:,:)
     SUx_o(:,:) = SUx_o(:,:)+GHW(ii)*GHW(jj)*F_out(ii,jj,:,:)*Vx(ii)
     SUy_o(:,:) = SUy_o(:,:)+GHW(ii)*GHW(jj)*F_out(ii,jj,:,:)*Vy(jj)
      SE_o(:,:) =  SE_o(:,:)+GHW(ii)*GHW(jj)*F_out(ii,jj,:,:)*(Vx(ii)**2+Vy(jj)**2)/2d0
    end do
    end do

 !!!!!!!!!!!!!!!!!!!!!!!!!!

          R_o    = SR_o
          Ux_o   = SUx_o/SR_o
          Uy_o   = SUy_o/SR_o
          ET_o   = SE_o

  select case(IT)

    case (0) ! MAXWELLIAN
         do ii=0, NDo1
         do jj=0, NDo1

          T_o(ii,jj) = (2d0* ET_o(ii,jj)/R_o(ii,jj))-(Ux_o(ii,jj)**2+Uy_o(ii,jj)**2)
          Z_o(ii,jj) = R_o(ii,jj)/(pi*T_o(ii,jj))
          P_o(ii,jj) = R_o(ii,jj)*T_o(ii,jj)/2d0

         end do
         end do

       case(-1)
! maxwellian = 0., fermion = 1., boson = -1
         do ii=0, NDo1
         do jj=0, NDo1

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
    PSIA = 2d0*ET_o(ii,jj) - (GA2*(R_o(ii,jj)/GA1)**2)/pi - &
 R_o(ii,jj)*(Ux_o(ii,jj)**2+Uy_o(ii,jj)**2)
    PSIB = 2d0*ET_o(ii,jj) - (GB2*(R_o(ii,jj)/GB1)**2)/pi - &
 R_o(ii,jj)*(Ux_o(ii,jj)**2+Uy_o(ii,jj)**2)
    ZC = (ZA+ZB)/2d0
    GC1 = 0d0
    GC2 = 0d0
        DO L = 1, 100
            GC1 = GC1 + (ZC**L)/dble(L)
            GC2 = GC2 + (ZC**L)/dble(L**2)
        END DO
    PSIC = 2d0*ET_o(ii,jj) - (GC2*(R_o(ii,jj)/GC1)**2)/pi - &
R_o(ii,jj)*(Ux_o(ii,jj)**2+Uy_o(ii,jj)**2)

    IF ((PSIA*PSIC) .LT. 0) THEN
        ZB = ZC
    ELSE
        ZA = ZC
    END IF
  END DO
        Z_o(ii,jj) = ZC;
!        T_o(ii,jj) = R(ii,jj)**2 / (pi*GC1);
        T_o(ii,jj) = R_o(ii,jj) / (pi*GC1);
        P_o(ii,jj) = ET_o(ii,jj) - 0.5d0* R_o(ii,jj) *( Ux_o(ii,jj)**2+Uy_o(ii,jj)**2)
         end do
     end do
        case(1)
!           %(IT == 1)
!maxwellian = 0., fermion = 1., boson = -1
             do ii=0, NDo1
             do jj=0, NDo1

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
    PSIA = 2d0*ET_o(ii,jj) - (GA2*(R_o(ii,jj)/GA1)**2)/pi - &
 R_o(ii,jj)*(Ux_o(ii,jj)**2+Uy_o(ii,jj)**2)
    PSIB = 2d0*ET_o(ii,jj) - (GB2*(R_o(ii,jj)/GB1)**2)/pi - &
 R_o(ii,jj)*(Ux_o(ii,jj)**2+Uy_o(ii,jj)**2)
    ZC = (ZA+ZB)/2d0
    GC1 = 0d0
    GC2 = 0d0
        DO L = 1, 100
            GC1 = GC1 + (ZC**L) * (-1)**(L-1)/dble(L)
            GC2 = GC2 + (ZC**L) * (-1)**(L-1)/dble(L**2)
        END DO
    PSIC = 2d0*ET_o(ii,jj) - (GC2*(R_o(ii,jj)/GC1)**2)/pi - &
R_o(ii,jj)*(Ux_o(ii,jj)**2+Uy_o(ii,jj)**2)

    IF ((PSIA*PSIC) .LT. 0) THEN
        ZB = ZC
    ELSE
        ZA = ZC
    END IF
  END DO
        Z_o(ii,jj) = ZC;
!        T_o(ii,jj) = R(ii,jj)**2 / (pi*GC1);
        T_o(ii,jj) = R_o(ii,jj) / (pi*GC1);

        P_o(ii,jj) = ET_o(ii,jj) - 0.5d0* R_o(ii,jj) *( Ux_o(ii,jj)**2+Uy_o(ii,jj)**2)
         end do
     end do

            end select !if IT

    write(lid,*) 'ZONE T="Dm',DDK,'",i=',NDo1+1,',j=',NDo1+1

         do ii=0, NDo1
         do jj=0, NDo1

   write (lid,*) xo1(ii),xo2(jj), R_o(ii,jj), T_o(ii,jj), Z_o(ii,jj)

end do
end do
 enddo

 deallocate(SR_o,SUx_o,SUy_o,SE_o)
 deallocate(R_o,Ux_o,Uy_o,ET_o)
 deallocate(P_o,T_o,Z_o,F_out)
 deallocate(Leg_Grid_xi1_o,Leg_Grid_xi2_o)
 deallocate(leg_tb_o,PN,PD)
 deallocate(LG_grids_o, ValuesOfPolyNatGrids_o)
 deallocate(xo1,xo2)
end if
  close(lid)




! 1000 format(2e23.15)

 end subroutine Tecplt_Output
