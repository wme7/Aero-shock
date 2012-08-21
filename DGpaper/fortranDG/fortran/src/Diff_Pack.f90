subroutine diff_x1(N_max,N,A,DA,Domain_Num)
  use Legendre     ! Legendre.f90
  use Metric_Var   ! Metric_Var.f90
  implicit none
  ! Declare subroutine arguments
    integer:: N_max(1:2)
    integer:: N(1:2)
    integer:: Domain_Num
    real(kind=8)::  A(0:N_max(1),0:N_max(2))
    real(kind=8):: DA(0:N_max(1),0:N_max(2))
  !
  ! Declare local arguments
    integer:: i,j
    integer:: N1,N2
  !
    N1=N(1); N2=N(2);

    DA(0:N1,0:N2) = dxi1_dx1(0:N1,0:N2,Domain_Num) * &
         matmul(Diff_xi1(0:N1,0:N1,N1), A(0:N1,0:N2)) &
                  + dxi2_dx1(0:N1,0:N2,Domain_Num) * &
         matmul(A(0:N1,0:N2), Diff_xi2(0:N2,0:N2,N2)) 

  return

end subroutine diff_x1

!========================================================================

subroutine diff_x2(N_max,N,A,DA,Domain_Num)
  use Legendre     ! Legendre.f90
  use Metric_Var   ! Metric_Var.f90
  implicit none

  ! Declare subroutine arguments
    integer:: N_max(1:2)
    integer:: N(1:2)
    integer:: Domain_Num
    real(kind=8)::  A(0:N_max(1),0:N_max(2))
    real(kind=8):: DA(0:N_max(1),0:N_max(2))

  ! Declare local arguments
    integer:: i,j
    integer:: N1,N2

  ! subroutine begin
  !
    N1=N(1); N2=N(2);

    DA(0:N1,0:N2) = dxi1_dx2(0:N1,0:N2,Domain_Num) * &
         matmul(Diff_xi1(0:N1,0:N1,N1), A(0:N1,0:N2)) &
                  + dxi2_dx2(0:N1,0:N2,Domain_Num) * &
         matmul(A(0:N1,0:N2), Diff_xi2(0:N2,0:N2,N2)) 

  return

end subroutine diff_x2

!========================================================================
subroutine differentiate_xi1(N_max,N,A,DA)
  use Legendre     ! Legendre.f90
  implicit none
  ! Declare subroutine arguments
    integer:: N_max(1:2)
    integer:: N(1:2)
    real(kind=8)::  A(0:N_max(1),0:N_max(2))
    real(kind=8):: DA(0:N_max(1),0:N_max(2))
  !
  ! Declare local arguments
    integer:: i
    integer:: N1,N2
  !
  ! subroutine begin
  !
    N1=N(1); N2=N(2);

    DA(0:N1,0:N2)=matmul(Diff_xi1(0:N1,0:N1,N1),A(0:N1,0:N2))

  return

end subroutine differentiate_xi1
!====================================================================!

subroutine differentiate_xi2(N_max,N,A,DA)
  use Legendre     ! Legendre.f90
  implicit none
  ! Declare subroutine arguments
    integer:: N_max(1:2)
    integer:: N(1:2)
    real(kind=8)::  A(0:N_max(1),0:N_max(2))
    real(kind=8):: DA(0:N_max(1),0:N_max(2))
  !
  ! Declare local arguments
    integer:: i
    integer:: N1,N2
  !
  ! subroutine begin
  !
    N1=N(1); N2=N(2);
    DA(0:N1,0:N2)=matmul(A(0:N1,0:N2),Diff_xi2(0:N2,0:N2,N2))

  return

end subroutine differentiate_xi2


subroutine Inner_flux1_loop!(DK_Connect)
!  do inx_edg=1,in1_EG
!      EG_Connect = EG_in_Connect(1:5,inx_edg)

  use Legendre
  use MD2D_Grid
  use Metric_Var
  use NorVec_Var
  implicit none
  integer::DK_Connect(1:5),DK,i,j,ii,ierr,inx_edg
  real(kind=8),allocatable :: temp(:,:)
  real(kind=8),allocatable :: temp1(:,:),tempsum1(:,:)
  real(kind=8),allocatable :: temp3(:,:),tempsum3(:,:)
  allocate( temp(0:PDeg1,0:PDeg2), &
temp1(0:PDeg1,0:PDeg2),tempsum1(0:PDeg2,0:0),&
temp3(0:PDeg1,0:PDeg2),tempsum3(0:PDeg2,0:0), stat=ierr  )

do inx_edg=1,in1_EG
DK_connect = EG_in_Connect(1:5,inx_edg)


! edg 1
  DDK = DK_connect(2)
  DK = DK_connect(4)
do i=1,IGH
do j=1,IGH
      temp = transpose(F_alt(i,j,:,:,DDK) )
      temp1 = E_gen(2,1,i,j)*temp*C2x_Edge1

      temp = transpose(F_alt(i,j,:,:,DK) )
      temp1 = temp1 + E_gen(1,1,i,j)*temp*C1x_Edge1

      temp = transpose(F_alt(i,j,:,:,DK)*C2x_Edge3)
      temp3 = E_gen(2,3,i,j)*temp
      temp = transpose(F_alt(i,j,:,:,DDK )*C1x_Edge3)
      temp3 = temp3 + E_gen(1,3,i,j)*temp

      do ii=0,PDeg2
         tempsum1(ii,0) = sum ( temp1(:,ii) )
         tempsum3(ii,0) = sum ( temp3(:,ii) )
      enddo

      temp = matmul(tempsum1,Cy_Edge1)
      flx_F(i,j,:,:,DDK) = flx_F(i,j,:,:,DDK) + temp
      temp = matmul(tempsum3,Cy_Edge3)
      flx_F(i,j,:,:,DK) = flx_F(i,j,:,:,DK) + temp

enddo
enddo
enddo

  deallocate( temp,temp1,tempsum1,temp3,tempsum3)
end subroutine Inner_flux1_loop

subroutine Inner_flux2_loop
  use Legendre
  use MD2D_Grid
  use Metric_Var
  use NorVec_Var
  implicit none
  integer::DK_Connect(1:5),DK,i,j,ii,ierr,inx_edg
  real(kind=8),allocatable :: temp(:,:)
  real(kind=8),allocatable :: temp2(:,:),tempsum2(:,:)
  real(kind=8),allocatable :: temp4(:,:),tempsum4(:,:)
  logical isnan 

  allocate( temp(0:PDeg1,0:PDeg2), &
temp2(0:PDeg1,0:PDeg2),tempsum2(0:0,0:PDeg2),&
temp4(0:PDeg1,0:PDeg2),tempsum4(0:0,0:PDeg2), stat=ierr  )

!do inx_edg=1,in1_EG
do inx_edg=(1+in1_EG),(in1_EG+in2_EG)

DK_connect = EG_in_Connect(1:5,inx_edg)

! edg 2
  DDK = DK_connect(2)
  DK = DK_connect(4)

do i=1,IGH
do j=1,IGH

      temp2 = E_gen(2,2,i,j)*F_alt(i,j,:,:,DDK)*C1x_Edge1
      temp2 = temp2 + E_gen(1,2,i,j)*F_alt(i,j,:,:,DK )*C2x_Edge1

      temp4 = E_gen(2,4,i,j)*F_alt(i,j,:,:,DK)*C2x_Edge1
      temp4 = temp4 + E_gen(1,4,i,j)*F_alt(i,j,:,:,DDK )*C1x_Edge1


      do ii=0,PDeg2
         tempsum2(0,ii) = sum ( temp2(:,ii) )
         tempsum4(0,ii) = sum ( temp4(:,ii) )
      enddo

      temp = matmul(Cy_Edge2,tempsum2)
      flx_F(i,j,:,:,DDK) = flx_F(i,j,:,:,DDK) + temp
      temp = matmul(Cy_Edge4,tempsum4)
      flx_F(i,j,:,:,DK) = flx_F(i,j,:,:,DK) + temp

if (isnan(flx_F(i,j,1,1,DK))) then
write(6,*) "NaN flx in 2",i,j,DK,DDK
write(6,*) "DK=",DK,"F_alt",F_alt(i,j,0:3,0:3,DK)
write(6,*) "DDK=",DDK,"F_alt",F_alt(i,j,0:3,0:3,DDK)
write(6,*) tempsum2(0,0:PDeg2)
write(6,*) tempsum4(0,0:PDeg2)
stop
endif

enddo
enddo
enddo
  deallocate( temp,temp2,tempsum2,temp4,tempsum4)

end subroutine Inner_flux2_loop


subroutine Outedg1_flux_loop
  use Legendre
  use MD2D_Grid
  use Metric_Var
  use NorVec_Var
  implicit none
  integer::DK_Connect(1:5),DK,NDK,i,j,ii,ierr,inx_edg
  real(kind=8),allocatable :: temp(:,:)
  real(kind=8),allocatable :: temp1(:,:),tempsum1(:,:)
  real(kind=8),allocatable :: temp3(:,:),tempsum3(:,:)
  allocate( temp(0:PDeg1,0:PDeg2), &
temp1(0:PDeg1,0:PDeg2),tempsum1(0:PDeg2,0:0),&
temp3(0:PDeg1,0:PDeg2),tempsum3(0:PDeg2,0:0), stat=ierr  )

do inx_edg=1,out1_EG

DK_connect = EG_out_Connect(1:5,inx_edg)
! edg 1
  DDK = DK_connect(2)
  DK = DK_connect(4)

do i=1,IGH
do j=1,IGH

if (DK .gt. 0) then
write(6,*) "wrong BD_type",DDK,DK

else if (DK .lt. 0) then
!Periodic BC
NDK=-DK
      temp = transpose(F_alt(i,j,:,:,DDK) )
      temp1 = E_gen(2,1,i,j)*temp*C2x_Edge1

      temp = transpose(F_alt(i,j,:,:,NDK) )
      temp1 = temp1 + E_gen(1,1,i,j)*temp*C1x_Edge1

else
!Neumann BC

      temp = transpose(F_alt(i,j,:,:,DDK) )
      temp1 = (E_gen(1,1,i,j)+E_gen(2,1,i,j))*temp*C2x_Edge1

endif


      do ii=0,PDeg2
         tempsum1(ii,0) = sum ( temp1(:,ii) )
      enddo

      temp = matmul(tempsum1,Cy_Edge1)
      flx_F(i,j,:,:,DDK) = flx_F(i,j,:,:,DDK) + temp

enddo
enddo
enddo

  deallocate( temp,temp1,tempsum1,temp3,tempsum3)
end subroutine Outedg1_flux_loop

subroutine Outedg2_flux_loop
  use Legendre
  use MD2D_Grid
  use Metric_Var
  use NorVec_Var
  implicit none
  integer::DK_Connect(1:5),DK,NDK,ii,i,j,ierr,inx_edg
  real(kind=8),allocatable :: temp(:,:)
  real(kind=8),allocatable :: temp2(:,:),tempsum2(:,:)
  real(kind=8),allocatable :: temp4(:,:),tempsum4(:,:)

  allocate( temp(0:PDeg1,0:PDeg2), &
temp2(0:PDeg1,0:PDeg2),tempsum2(0:0,0:PDeg2),&
temp4(0:PDeg1,0:PDeg2),tempsum4(0:0,0:PDeg2), stat=ierr  )

  !do inx_edg=1,in1_EG

  do inx_edg=(1+out1_EG+out3_EG),(out1_EG+out2_EG+out3_EG)

DK_connect = EG_out_Connect(1:5,inx_edg)

! edg 2
  DDK = DK_connect(2)
  DK = DK_connect(4)

do i=1,IGH
do j=1,IGH

if (DK .gt. 0) then
write(6,*) "wrong BD_type",DDK,DK

else if (DK .lt. 0) then
!Periodic BC
NDK=-DK
      temp2 = E_gen(2,2,i,j)*F_alt(i,j,:,:,DDK)*C1x_Edge1
      temp2 = temp2 + E_gen(1,2,i,j)*F_alt(i,j,:,:,NDK )*C2x_Edge1

else
!Neumann BC
      temp2 = (E_gen(1,2,i,j)+E_gen(2,2,i,j))*F_alt(i,j,:,:,DDK)*C1x_Edge1


end if


      do ii=0,PDeg2
         tempsum2(0,ii) = sum ( temp2(:,ii) )
      enddo

      temp = matmul(Cy_Edge2,tempsum2)
      flx_F(i,j,:,:,DDK) = flx_F(i,j,:,:,DDK) + temp

enddo
enddo
enddo

  deallocate( temp,temp2,tempsum2,temp4,tempsum4)

end subroutine Outedg2_flux_loop

subroutine Outedg3_flux_loop!(DK_Connect)
  use Legendre
  use MD2D_Grid
  use Metric_Var
  use NorVec_Var
  implicit none
  integer::DK_Connect(1:5),DK,NDDK,i,j,ii,ierr,inx_edg
  real(kind=8),allocatable :: temp(:,:)
  real(kind=8),allocatable :: temp1(:,:),tempsum1(:,:)
  real(kind=8),allocatable :: temp3(:,:),tempsum3(:,:)
  allocate( temp(0:PDeg1,0:PDeg2), &
temp1(0:PDeg1,0:PDeg2),tempsum1(0:PDeg2,0:0),&
temp3(0:PDeg1,0:PDeg2),tempsum3(0:PDeg2,0:0), stat=ierr  )

  do inx_edg=(1+out1_EG),(out1_EG+out3_EG)

DK_connect = EG_out_Connect(1:5,inx_edg)
! edg 3
  DDK = DK_connect(4)
  DK = DK_connect(2)

do i=1,IGH
do j=1,IGH

if (DDK .gt. 0) then
write(6,*) "wrong BD_type",DDK

else if (DDK .lt. 0) then
!Periodic BC
NDDK=-DDK

      temp = transpose(F_alt(i,j,:,:,DK)*C2x_Edge3)
      temp3 = E_gen(2,3,i,j)*temp

      temp = transpose(F_alt(i,j,:,:,NDDK )*C1x_Edge3)
      temp3 = temp3 + E_gen(1,3,i,j)*temp

else
!Neumann BC
      temp = transpose(F_alt(i,j,:,:,DK)*C2x_Edge3)
      temp3 = (E_gen(1,3,i,j)+E_gen(2,3,i,j))*temp

end if

      do ii=0,PDeg2
         tempsum3(ii,0) = sum ( temp3(:,ii) )
      enddo

      temp = matmul(tempsum3,Cy_Edge3)
      flx_F(i,j,:,:,DK) = flx_F(i,j,:,:,DK) + temp

enddo
enddo
enddo

  deallocate( temp,temp1,tempsum1,temp3,tempsum3)
end subroutine Outedg3_flux_loop

subroutine Outedg4_flux_loop
  use Legendre
  use MD2D_Grid
  use Metric_Var
  use NorVec_Var
  implicit none
  integer::DK_Connect(1:5),DK,NDDK,i,j,ii,ierr,inx_edg
  real(kind=8),allocatable :: temp(:,:)
  real(kind=8),allocatable :: temp2(:,:),tempsum2(:,:)
  real(kind=8),allocatable :: temp4(:,:),tempsum4(:,:)

  allocate( temp(0:PDeg1,0:PDeg2), &
temp2(0:PDeg1,0:PDeg2),tempsum2(0:0,0:PDeg2),&
temp4(0:PDeg1,0:PDeg2),tempsum4(0:0,0:PDeg2), stat=ierr  )

  do inx_edg=(1+out1_EG+out2_EG+out3_EG),(out1_EG+out2_EG+out3_EG+out4_EG)


DK_connect = EG_out_Connect(1:5,inx_edg)

! edg 4
  DDK = DK_connect(4)
  DK = DK_connect(2)

do i=1,IGH
do j=1,IGH

if (DDK .gt. 0) then
write(6,*) "wrong BD_type"
else if (DDK .lt. 0) then
!Periodic BC
NDDK=-DDK
      temp4 = E_gen(2,4,i,j)*F_alt(i,j,:,:,DK)*C2x_Edge1
      temp4 = temp4 + E_gen(1,4,i,j)*F_alt(i,j,:,:,NDDK )*C1x_Edge1

else
!Neumann BC
      temp4 = (E_gen(1,4,i,j)+E_gen(2,4,i,j))*F_alt(i,j,:,:,DK)*C2x_Edge1

!      temp4 = temp4 + E(1,4,DK)*u_alt(:,:,DDK )*C1x_Edge1

end if


      do ii=0,PDeg2
         tempsum4(0,ii) = sum ( temp4(:,ii) )
      enddo

      temp = matmul(Cy_Edge4,tempsum4)
      flx_F(i,j,:,:,DK) = flx_F(i,j,:,:,DK) + temp

enddo
enddo
enddo

  deallocate( temp,temp2,tempsum2,temp4,tempsum4)

end subroutine Outedg4_flux_loop


!====================================================================!
subroutine DG2D_flux
  use Legendre
  use MD2D_Grid
  use Metric_Var
  use NorVec_Var
  use State_Var
  implicit none
  integer :: i,j,ii,ierr,DK_Connect(1:4),EG_Connect(1:5),inx_edg
  real(kind=8),allocatable :: edge_flux(:,:,:),temp(:,:)
  real(kind=8),allocatable :: temp1(:,:),tempsum1(:,:),temp2(:,:),tempsum2(:,:)
  real(kind=8),allocatable :: temp3(:,:),tempsum3(:,:),temp4(:,:),tempsum4(:,:)
  real(kind=8),allocatable :: tempb1(:,:),tempb2(:,:)
  integer:: e_flag

flx_F=0d0

  allocate( tempb1(0:PDeg1,0:PDeg2),tempb2(0:PDeg1,0:PDeg2), &
            stat=ierr  )

!---------------------------------------------
! Interior Edges
!---------------------------------------------
!write(6,*) "F_alt",F_alt(1:3,1,0:2,0:2,3)

      call Compute_Source
      call Inner_flux1_loop!(EG_Connect)
!write(6,*) "Inner 1 flx_F",flx_F(1:3,1,0:2,0:2,3)

      call Inner_flux2_loop!(EG_Connect)
!write(6,*) "Inner 2 flx_F",flx_F(1:3,1,0:2,0:2,3)

!---------------------------------------------
! Outer Edges
!---------------------------------------------

if (TotNum_out_EG .gt. 0) then
call Outedg1_flux_loop
!write(6,*) "Out 1 flx_F",flx_F(1:3,1,0:2,0:2,3)

call Outedg3_flux_loop
call Outedg2_flux_loop
call Outedg4_flux_loop
endif


!---------------------------------------------
! End of Edge Flux
!---------------------------------------------

  do DDK=1,TotNum_DM
  do i=1,IGH
   do j=1, IGH

      tempb1 = matmul(transpose(B_tal1x),F_alt(i,j,:,:,DDK))
      tempb2 = matmul(F_alt(i,j,:,:,DDK),B_tal1x) 

!      u_t(:,:,DDK) = B(1,DDK)*(tempb1*B_tal1y) + B(2,DDK)*(tempb2*B_tal2x)
!      F_t(i,j,:,:,DDK) = (a_vec(1)*Bq(2,2,DDK)-a_vec(2)*Bq(1,2,DDK))*(tempb1*B_tal1y)+&
!                    (-a_vec(1)*Bq(2,1,DDK)+a_vec(2)*Bq(1,1,DDK))*(tempb2*B_tal2x)

      F_t(i,j,:,:,DDK) = (GH(i)*Bq(2,2,DDK)-GH(j)*Bq(1,2,DDK))*(tempb1*B_tal1y)+&
                    (-GH(i)*Bq(2,1,DDK)+GH(j)*Bq(1,1,DDK))*(tempb2*B_tal2x)

      F_t(i,j,:,:,DDK) = F_t(i,j,:,:,DDK) - flx_F(i,j,:,:,DDK)-FS(i,j,:,:,DDK)
      F_t(i,j,:,:,DDK) = F_t(i,j,:,:,DDK) *A1(:,:,DDK)
   enddo
  enddo
  enddo
  deallocate( tempb1,tempb2 )
endsubroutine DG2D_flux

!====================================================================
subroutine DG_compute_error(error_L2,error_max,time)!,a_vec)
  use Legendre
  use MD2D_Grid
  use Metric_Var
  use NorVec_Var
  use universal_const
  use State_Var
  implicit none
  integer:: i,j,ierr,ii,jj
  real(kind=8):: error_L2,error_max,tempmax,time!,a_vec(2)
  
  real(kind=8),allocatable :: uhhij(:,:,:,:,:),tempsum(:,:,:),temp(:,:,:,:)
  real(kind=8),allocatable :: error_L2ij(:,:),error_maxij(:,:)

  allocate( uhhij(1:IGH,1:IGH,0:PND1,0:PND2,1:TotNum_DM),&
     tempsum(1:IGH,1:IGH,0:PND2),temp(1:IGH,1:IGH,0:PND1,0:PND2),&
     error_L2ij(1:IGH,1:IGH),error_maxij(1:IGH,1:IGH),&
     stat=ierr  )

  if (I_prob .ge. 2) then
  do DDK=1,TotNum_DM
   do ii=1,IGH
   do jj=1,IGH
     do i=0,PND1
        do j=0,PND2
           exa_solij(ii,jj,i,j,DDK) = dcos(2*pi*(x1(i,j,DDK)+x2(i,j,DDK)-(GH(ii)+GH(jj))*time))
        enddo
     enddo
        enddo
     enddo
  enddo
  else
  do DDK=1,TotNum_DM
   do ii=1,IGH
   do jj=1,IGH
     do i=0,PND1
        do j=0,PND2
           exa_solij(ii,jj,i,j,DDK) = dsin(2*pi*(x1(i,j,DDK)+x2(i,j,DDK)-(GH(ii)+GH(jj))*time))
        enddo
     enddo
        enddo
     enddo
  enddo
  endif
  uhhij = 0.d0  

  do DDK=1,TotNum_DM
   do ii=1,IGH
   do jj=1,IGH

     do i=0,PDeg1
        do j=0,PDeg2
      uhhij(ii,jj,:,:,DDK)=uhhij(ii,jj,:,:,DDK)+F_new(ii,jj,i,j,DDK)*(Leg_Grid_xi1(:,:,i)*Leg_Grid_xi2(:,:,j))    
        enddo 
     enddo
        enddo
     enddo
  enddo

  error_L2ij = 0.d0
  error_maxij = 0.d0

  do DDK=1,TotNum_DM
   do ii=1,IGH
   do jj=1,IGH
     temp(ii,jj,:,:) = (exa_solij(ii,jj,:,:,DDK)-uhhij(ii,jj,:,:,DDK))**2
     temp(ii,jj,:,:) = temp(ii,jj,:,:)*Jacobian(:,:,DDK)*LGLWeights_Grid_xi1(:,:)*&
       LGLWeights_Grid_xi2(:,:)
     do i=0,PND2  
        tempsum(ii,jj,i)=sum(temp(ii,jj,:,i))
     enddo
     error_L2ij(ii,jj) = error_L2ij(ii,jj) + sum(tempsum(ii,jj,:))

     temp(ii,jj,:,:) = abs(exa_solij(ii,jj,:,:,DDK)-uhhij(ii,jj,:,:,DDK))
     do i=0,PND1
        tempsum(ii,jj,i)=maxval(temp(ii,jj,:,i))
     enddo
     tempmax = maxval(tempsum(ii,jj,:))
     error_maxij(ii,jj) = max(error_maxij(ii,jj),tempmax) 

   enddo
enddo

  enddo

  error_L2 = sqrt(error_L2ij(2,1)) 
  error_max =  error_maxij(2,1)
  deallocate( uhhij,tempsum,temp )

endsubroutine DG_compute_error
!====================================================================

