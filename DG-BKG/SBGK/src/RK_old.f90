
subroutine BGKimexL_nsoli(xcur, fcur)
!(n, xcur, fcur, rpar, ipar, itrmf)
!     call DGimexL(NV*nnx*pp, Fi, rwork, rpar, ipar, itrmf)
use State_Var, ONLY: F, NV, IT, VIS, V, VW,dx_g,dt_g,nx_g, pp_g
use Legendre
use RK_Var

implicit none

real(kind=8), intent(inout):: xcur(1:NV*nx_g*pp_g)
real(kind=8),intent(out) :: fcur(1:NV*nx_g*pp_g)

integer :: pp !=ipar(2)
integer :: nnx !=ipar(1)
real(kind=8):: dx 
real(kind=8):: dt 
real(kind=8):: gamma_ark

integer :: p,i,j
integer :: K,m
real(kind=8):: FC(1:pp_g), FB(1:pp_g)
real(kind=8):: FS_loc(1:NV,1:nx_g,1:pp_g)
real(kind=8):: vin(1:NV,1:nx_g,1:pp_g)
real(kind=8):: ft_a(1:NV,1:nx_g,1:pp_g)
real(kind=8):: FEQ_l(1:NV,1:nx_g,1:pp_g)
real(kind=8):: Mtemp(NV,pp_g)
real(kind=8):: F_loc(1:NV,1:pp_g)

real(kind=8):: SR(1:nx_g,1:pp_g), SU
real(kind=8):: SE,SAV
real(kind=8):: Rloc(1:nx_g,1:pp_g), Uloc(1:nx_g,1:pp_g)
real(kind=8):: Zloc(1:nx_g,1:pp_g), Tloc(1:nx_g,1:pp_g)
real(kind=8):: ETloc(1:nx_g,1:pp_g)
logical isnan

!NV=ipar(3)
 pp=pp_g
nnx=nx_g
dx = dx_g
dt = dt_g

gamma_ark=const_a_I(2,1)
!gamma_ark=rpar(3)

write(6,*) "BGK",gamma_ark
p=pp-1
vin=reshape(xcur,(/ NV,nnx,pp /) )
write(6,*) "RU, size of vin",size(vin),shape(vin)
!call Comp_RU_F(nnx,pp,vin)
do i=1,nnx
    Mtemp=0d0
    do K=1,NV
        Mtemp(K,:)=vin(K,i,:)
    enddo

    F_loc(:,:)=matmul(Mtemp,Pleg)
        SR(i,:) = matmul(VW,F_loc)

    do m=1,pp
!        SR(i,:) = matmul(VW,F_loc)
        SU = sum( VW*F_loc(:,m)* V)
        SE = sum( VW*F_loc(:,m)* V**2)/2d0
        SAV= sum( VW*F_loc(:,m)* abs(V))

        Rloc(i,m)    = SR(i,m)
        Uloc(i,m)    = SU/SR(i,m)
        ETloc(i,m)   = SE
!        AV(i,m)   = SAV
    enddo
enddo
!write(6,*) "ZTP b()",b
!call Comp_ZTP_F(nnx,pp)
    do i=1,nnx
        do m=1,pp
            Tloc(i,m)    = 4d0*ETloc(i,m)/Rloc(i,m) - 2d0*Uloc(i,m)**2
            Zloc(i,m)    = Rloc(i,m) / sqrt(pi* Tloc(i,m))
            !Ps(i,m) = ET(i,m) - 0.5d0 * Rs(i,m) * Us(i,m)**2
        enddo
    enddo
!write(6,*) "FEQ"
!call EQUILIBRIUM_F(nnx,pp)
   do i=1,nnx
    do K=1,NV
     do m=1,pp
    FEQ_l(K,i,m)   = 1d0/((exp( (V(K)-Uloc(i,m))**2 /Tloc(i,m))/Zloc(i,m)) + dble(IT) )
    if (FEQ_l(K,i,m) .gt. 1d20) then
    write(6,*) "Infty",i,K,m
    stop
    endif
     END DO
  END DO
END DO

!write(6,*) "Source", FS_loc
!call Compute_Source_F(nnx,pp,dx,vin,FS_loc)
      do i=1,nnx
        do K = 1,NV

         FC(:)=vin(K,i,:)
         FB(:)=FEQ_l(K,i,:)
         FC=matmul(FC,Pleg)-FB
          do j=0,p
           FS_loc(K,i,j+1)=sum(FC*Pleg(j+1,:)*wl)*dx/2d0/VIS(i,j+1)
          enddo
        enddo
      enddo

write(6,*) "Comp NV,nx,pp", NV, nnx, pp,maxval(abs(reshape(FS_loc,(/ NV*nnx*pp /) )))

do K=1,NV
  do j=1,pp
      FS_loc(K,:,j)=-FS_loc(K,:,j)*b(j)
!       do i=1,nnx
!       if (isnan(FS_loc(K,i,j))) then
!        write(6,*) "FS NaN K=",K,"Ij",i,j
!      stop
!      endif
!       enddo
  enddo
enddo

write(6,*) "FS",maxval(abs(reshape(FS_loc,(/ NV*nnx*pp /) )))

ft_a= vin - F - gamma_ark*dt*FS_loc

fcur=reshape(ft_a, (/ NV*nnx*pp /) )
write(6,*) "fcur", maxval(abs(fcur))

return
end subroutine

subroutine BGKimexL(n, xcur, fcur, rpar, ipar, itrmf)
!     call DGimexL(NV*nnx*pp, Fi, rwork, rpar, ipar, itrmf)
use State_Var, ONLY: F, NV, IT, VIS, V, VW
use Legendre
use RK_Var

implicit none
integer,intent(in):: n
real(kind=8),intent(in):: rpar(3)
integer,intent(in) :: ipar(1:2)

integer,intent(out) :: itrmf
real(kind=8), intent(inout):: xcur(1:n)
real(kind=8),intent(out) :: fcur(1:n)

integer :: pp !=ipar(2)
integer :: nnx !=ipar(1)
real(kind=8):: dx != rpar(1)
real(kind=8):: dt != rpar(2)
real(kind=8):: gamma_ark

integer :: p,i,j
integer :: K,m
real(kind=8):: FC(1:ipar(2)), FB(1:ipar(2))
real(kind=8):: FS_loc(1:NV,1:ipar(1),1:ipar(2))
real(kind=8):: vin(1:NV,1:ipar(1),1:ipar(2))
real(kind=8):: ft_a(1:NV,1:ipar(1),1:ipar(2))
real(kind=8):: FEQ_l(1:NV,1:ipar(1),1:ipar(2))
real(kind=8):: Mtemp(NV,ipar(2))
real(kind=8):: F_loc(1:NV,1:ipar(2))

real(kind=8):: SR(1:ipar(1),1:ipar(2)), SU
real(kind=8):: SE,SAV
real(kind=8):: Rloc(1:ipar(1),1:ipar(2)), Uloc(1:ipar(1),1:ipar(2))
real(kind=8):: Zloc(1:ipar(1),1:ipar(2)), Tloc(1:ipar(1),1:ipar(2))
real(kind=8):: ETloc(1:ipar(1),1:ipar(2))
logical isnan

!NV=ipar(3)
 pp=ipar(2)
nnx=ipar(1)
dx = rpar(1)
dt = rpar(2)

gamma_ark=const_a_I(2,1)
!gamma_ark=rpar(3)

write(6,*) "BGK",ipar,rpar,gamma_ark
p=pp-1
vin=reshape(xcur,(/ NV,nnx,pp /) )
write(6,*) "RU, size of vin",size(vin),shape(vin)
!call Comp_RU_F(nnx,pp,vin)
do i=1,nnx
    Mtemp=0d0
    do K=1,NV
        Mtemp(K,:)=vin(K,i,:)
    enddo

    F_loc(:,:)=matmul(Mtemp,Pleg)
        SR(i,:) = matmul(VW,F_loc)

    do m=1,pp
!        SR(i,:) = matmul(VW,F_loc)
        SU = sum( VW*F_loc(:,m)* V)
        SE = sum( VW*F_loc(:,m)* V**2)/2d0
        SAV= sum( VW*F_loc(:,m)* abs(V))

        Rloc(i,m)    = SR(i,m)
        Uloc(i,m)    = SU/SR(i,m)
        ETloc(i,m)   = SE
!        AV(i,m)   = SAV
    enddo
enddo
!write(6,*) "ZTP b()",b
!call Comp_ZTP_F(nnx,pp)
    do i=1,nnx
        do m=1,pp
            Tloc(i,m)    = 4d0*ETloc(i,m)/Rloc(i,m) - 2d0*Uloc(i,m)**2
            Zloc(i,m)    = Rloc(i,m) / sqrt(pi* Tloc(i,m))
            !Ps(i,m) = ET(i,m) - 0.5d0 * Rs(i,m) * Us(i,m)**2
        enddo
    enddo
!write(6,*) "FEQ"
!call EQUILIBRIUM_F(nnx,pp)
   do i=1,nnx
    do K=1,NV
     do m=1,pp
    FEQ_l(K,i,m)   = 1d0/((exp( (V(K)-Uloc(i,m))**2 /Tloc(i,m))/Zloc(i,m)) + dble(IT) )
    if (FEQ_l(K,i,m) .gt. 1d20) then
    write(6,*) "Infty",i,K,m
    stop
    endif
     END DO
  END DO
END DO

!write(6,*) "Source", FS_loc
!call Compute_Source_F(nnx,pp,dx,vin,FS_loc)
      do i=1,nnx
        do K = 1,NV

         FC(:)=vin(K,i,:)
         FB(:)=FEQ_l(K,i,:)
         FC=matmul(FC,Pleg)-FB
          do j=0,p
           FS_loc(K,i,j+1)=sum(FC*Pleg(j+1,:)*wl)*dx/2d0/VIS(i,j+1)
          enddo
        enddo
      enddo

write(6,*) "Comp NV,nx,pp", NV, nnx, pp,maxval(abs(reshape(FS_loc,(/ NV*nnx*pp /) )))

do K=1,NV
  do j=1,pp
      FS_loc(K,:,j)=-FS_loc(K,:,j)*b(j)
!       do i=1,nnx
!       if (isnan(FS_loc(K,i,j))) then
!        write(6,*) "FS NaN K=",K,"Ij",i,j
!      stop
!      endif
!       enddo
  enddo
enddo

write(6,*) "FS",maxval(abs(reshape(FS_loc,(/ NV*nnx*pp /) )))

ft_a= vin - F - gamma_ark*dt*FS_loc

fcur=reshape(ft_a, (/ NV*nnx*pp /) )
itrmf = 0
write(6,*) "fcur", maxval(abs(fcur))

return
end subroutine

subroutine jacvDG(n,xcur,fcur,ijob,v,z,rpar,ipar,itrmjv)
! --------------------------------------------------------------------
! This subroutine evaluates Jacobian-vector products and preconditioner
! solves for the Bratu test problem.
! --------------------------------------------------------------------
      implicit none
      integer n, ijob, ipar(*), itrmjv
      double precision xcur(n), fcur(n), v(n), z(n), rpar(*)
      double precision cl, cr, dexp, h2l
      integer i, j, j1, j2, nx, ny

write(6,*) "ijob", ijob

      itrmjv = 0

return
end subroutine

! --------------------------------------------------------------------


! --------------------------------------------------------------------
  subroutine ARK_marching_1D(nnx,pp,l,dt) 
    use State_Var
    use RK_Var
    implicit none    
    integer:: nnx ,pp ,l,K 
    real(kind=8):: dt 
! Local Variables
    integer:: i,j

    if (l .LT. RK_Stage) then

        F_new=F_new+ dt*const_b(l)*(F_s(:,:,:,l)+F_ns(:,:,:,l))
        F=Fold

         do j=1,l 
          F = F + dt*(const_a_I(l+1,j)*F_s(:,:,:,j) + const_a_E(l+1,j)*F_ns(:,:,:,j))
         enddo

    else
          F_new=F_new+ dt*const_b(l)*(F_s(:,:,:,l)+F_ns(:,:,:,l))
    endif

     call Comp_RU(nnx,pp)
     call Comp_ZTP(nnx,pp)
  end subroutine ARK_marching_1D

  !------------------------------------------------------------------------------

  subroutine DG_flux_ark(nnx,pp,l,dx) !,c,b,A,F,u_t)
    use State_Var, ONLY: NV,F, F_s, F_ns, V, FS, BC_type
    use RK_Var
    use Legendre
    use Nitsol_Var
    implicit none

  integer:: nnx ,pp ,l
  real(kind=8):: dx

  integer::NVh
  integer:: i,j,K
  real(kind=8):: tempAu(1:pp), FC(1:pp), FR(1:pp),FU(1:pp)

  real(kind=8):: Fi(NV*nnx*pp),fout(NV*nnx*pp)
  real(kind=8):: fnrm, ftol,gamma_ark,dtt
  real(kind=8)::rpari, tol(2),it_hist(NV*nnx*pp,3)
  integer:: ipari(4)

  integer     :: itrmf, iterm
  integer:: info(6)

  double precision ddot, dnrm2
  external ddot, dnrm2
  external BGKIMEXL,JACVDG, BGKIMEXL_nsoli

gamma_ark=const_a_I(2,1)
dtt=rpar(2)

NVh=NV/2
tol=dx*1d-6

ipari=(/ 40, 40, 2, 20 /)
rpari=-0.1d0


!  Compute the source term 
call EQUILIBRIUM(nnx,pp)
call Compute_Source(nnx,pp,dx)

if (l .eq. 1) then

do K=1,NVh
! Right-going part
  tempAu=matmul(F(K,1,1:pp),A)
 if (BC_type .eq. 0) then
  !BC no-flux
        FC(:)=F(K,1,:)
        FU=FC
        FR(:)=FS(K,1,:)
!     F_s(K,1,:,1)=( (-FR)' .* b)
!     F_ns(K,1,:,1)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FU'.* c) * c')' .* b)
!     F_ns(K,1,:,1)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FU'.* c) * c')' .* b)+( (-FR)' .* b)

      F_s(K,1,:,l)=-FR
      F_ns(K,1,:,l)=V(K)*(tempAu(:)-sum(FC(1:pp))+sum( FU*c )*c)
  
      elseif (BC_type .eq.  -1) then
   !BC reflecting
        FC(:)=F(K,1,:)
        FU(:)=F(NV-K+1,1,:)
        FR(:)=FS(K,1,:)
!        F_s(K,1,:,1)=( (-FR)' .* b)
!        F_ns(K,1,:,1)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FU'.* c) * c')' .* b)

      F_s(K,1,:,l)=-FR
      F_ns(K,1,:,l)=V(K)*(tempAu(:)-sum(FC(1:pp))+sum( FU*c )*c)
      else
     write(6,*) "BC_type not supported",BC_type
     stop
   endif
    
  do i=2,nnx
  tempAu=matmul(F(K,i,1:pp),A)
      FU(:)=F(K,i-1,:)
      FC(:)=F(K,i,:)
      FR(:)=FS(K,i,:)

      F_s(K,i,:,l)=-FR
      F_ns(K,i,:,l)=V(K)*(tempAu(:)-sum(FC(1:pp))+sum(FU)*c)
  enddo
! Left-going part
    
  do i=1,nnx-1
  tempAu=matmul(F(NVh+K,i,1:pp),A)
      FU(:)=F(NVh+K,i+1,:)
      FC(:)=F(NVh+K,i,:)
      FR(:)=FS(NVh+K,i,:)
      F_s(NVh+K,i,:,l)=-FR
      F_ns(NVh+K,i,:,l)=V(NVh+K)*(tempAu(:)-sum(FU(1:pp)*c)+sum( FC*c )*c)
  enddo

  tempAu=matmul(F(NVh+K,nnx,1:pp),A)
 if (BC_type .eq. 0) then
  !BC no-flux
        FC(:)=F(NVh+K,nnx,:)
        FU=FC
        FR(:)=FS(NVh+K,nnx,:)

      F_s(NVh+K,nnx,:,l)=-FR
      F_ns(NVh+K,nnx,:,l)=V(NVh+K)*(tempAu(:)-sum(FU(1:pp))+sum( FC*c )*c)
  
      elseif (BC_type .eq.  -1) then
   !BC reflecting
        FC(:)=F(NVh+K,nnx,:)
        FU(:)=F(NV-K+1,nnx,:)
        FR(:)=FS(NVh+K,nnx,:)

      F_s(NVh+K,nnx,:,l)=-FR
      F_ns(NVh+K,nnx,:,l)=V(NVh+K)*(tempAu(:)-sum(FU(1:pp))+sum( FC*c )*c)
      else
     write(6,*) "BC_type not supported",BC_type
     stop
   endif

enddo !K 

     do j=1,pp
        F_s(:,:,j,l)=F_s(:,:,j,l)*b(j)
        F_ns(:,:,j,l)=F_ns(:,:,j,l)*b(j)
     enddo

else
! Stages 2-6

Fi=reshape(F,(/ NV*nnx*pp /) )
write(6,*) "size Fi", size(Fi)
     call BGKimexL(NV*nnx*pp, Fi, fout, rpar, ipar, itrmf)

     fnrm = dnrm2(NV*nnx*pp, fout, 1)
     ftol = rlftol*1d0 !fnrm
call NSOLI(Fi,BGKimexL_nsoli,tol, ipari, rpari, Fi, it_hist, iterm, NV*nnx*pp)

!call nitsol(NV*nnx*pp, Fi, BGKimexL, jacvDG, ftol, stptol,  &
!    input, info, rwork, rpar, ipar, iterm, ddot, dnrm2)

if (iterm .ne. 0) then
    write(6,*) "Iterm", iterm,ftol,maxval(abs(rwork)),info
stop
endif

    F=reshape(Fi,(/ NV,nnx,pp /))


do K=1,NVh
! Right-going part
  tempAu=matmul(F(K,1,1:pp),A)
 if (BC_type .eq. 0) then
  !BC no-flux
        FC(:)=F(K,1,:)
        FU=FC
        FR(:)=FS(K,1,:)

      F_s(K,1,:,l)=-FR
      F_ns(K,1,:,l)=V(K)*(tempAu(:)-sum(FC(1:pp))+sum( FU*c )*c)
  
      elseif (BC_type .eq.  -1) then
   !BC reflecting
        FC(:)=F(K,1,:)
        FU(:)=F(NV-K+1,1,:)
        FR(:)=FS(K,1,:)
      F_s(K,1,:,l)=-FR
      F_ns(K,1,:,l)=V(K)*(tempAu(:)-sum(FC(1:pp))+sum( FU*c )*c)
      else
     write(6,*) "BC_type not supported",BC_type
     stop
   endif
    
  do i=2,nnx
  tempAu=matmul(F(K,i,1:pp),A)
      FU(:)=F(K,i-1,:)
      FC(:)=F(K,i,:)
      FR(:)=FS(K,i,:)
      F_s(K,i,:,l)=-FR
      F_ns(K,i,:,l)=V(K)*(tempAu(:)-sum(FC(1:pp))+sum(FU)*c)
  enddo
! Left-going part
    
  do i=1,nnx-1
  tempAu=matmul(F(NVh+K,i,1:pp),A)
      FU(:)=F(NVh+K,i+1,:)
      FC(:)=F(NVh+K,i,:)
      FR(:)=FS(NVh+K,i,:)
      F_s(NVh+K,i,:,l)=-FR
      F_ns(NVh+K,i,:,l)=V(NVh+K)*(tempAu(:)-sum(FU(1:pp)*c)+sum( FC*c )*c)
  enddo

  tempAu=matmul(F(NVh+K,nnx,1:pp),A)
 if (BC_type .eq. 0) then
  !BC no-flux
        FC(:)=F(NVh+K,nnx,:)
        FU=FC
        FR(:)=FS(NVh+K,nnx,:)
      F_s(NVh+K,nnx,:,l)=-FR
      F_ns(NVh+K,nnx,:,l)=V(NVh+K)*(tempAu(:)-sum(FU(1:pp))+sum( FC*c )*c)
  
      elseif (BC_type .eq.  -1) then
   !BC reflecting
        FC(:)=F(NVh+K,nnx,:)
        FU(:)=F(NV-K+1,nnx,:)
        FR(:)=FS(NVh+K,nnx,:)
      F_s(NVh+K,nnx,:,l)=-FR
      F_ns(NVh+K,nnx,:,l)=V(NVh+K)*(tempAu(:)-sum(FU(1:pp))+sum( FC*c )*c)
      else
     write(6,*) "BC_type not supported",BC_type
     stop
   endif

enddo !K 
     do j=1,pp
        F_s(:,:,j,l)=F_s(:,:,j,l)*b(j)
        F_ns(:,:,j,l)=F_ns(:,:,j,l)*b(j)
     enddo

endif !l (RK)

  end subroutine DG_flux_ark


  !------------------------------------------------------------------------------
  subroutine ERK_marching_1D(nnx,pp,l,dt) 
    use State_Var
    use RK_Var
    implicit none    
    integer:: nnx ,pp ,l,K 
    real(kind=8):: dt 
! Local Variables
    integer:: i,j

    if (l .LT. RK_Stage) then

        F_new=F_new+ dt*const_b(l)*(F_s(:,:,:,l)+F_ns(:,:,:,l))
        F=Fold

         do j=1,l 
          F = F + dt*(const_a_I(l+1,j)*F_s(:,:,:,j) + const_a_E(l+1,j)*F_ns(:,:,:,j))
         enddo

    else
          F_new=F_new+ dt*const_b(l)*(F_s(:,:,:,l)+F_ns(:,:,:,l))
    endif

     call Comp_RU(nnx,pp)
     call Comp_ZTP(nnx,pp)
  end subroutine ERK_marching_1D

  !------------------------------------------------------------------------------

  subroutine DG_flux_and_bdd(nnx,pp,l,dx) !,c,b,A,F,u_t)
    use State_Var
    use RK_Var
    use Legendre
    implicit none

  integer:: nnx ,pp ,l
  real(kind=8):: dx

  integer::NVh
  integer:: i,j,K
  real(kind=8):: tempAu(1:pp), FC(1:pp), FR(1:pp),FU(1:pp)
 !, FS(1:nnx,1:pp) 

NVh=NV/2

!  Compute the source term 
call EQUILIBRIUM(nnx,pp)
call Compute_Source(nnx,pp,dx)

!do i=1,nnx
!      FC(:)=matmul(F(i,:),Pleg)
  !      FC(:)=FC*FC
  !      FC=(FC'*Pleg).^2
! do j=1,pp
!   FS(i,j)=sum (FC(:)*Pleg(j,:)*wl(:))*dx/2d0
! end do            
!end do

do K=1,NVh
! Right-going part
  tempAu=matmul(F(K,1,1:pp),A)
 if (BC_type .eq. 0) then
  !BC no-flux
        FC(:)=F(K,1,:)
        FU=FC
        FR(:)=FS(K,1,:)
!     F_s(K,1,:,1)=( (-FR)' .* b)
!     F_ns(K,1,:,1)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FU'.* c) * c')' .* b)
!     F_ns(K,1,:,1)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FU'.* c) * c')' .* b)+( (-FR)' .* b)

      F_ns(K,1,:,l)=V(K)*(tempAu(:)-sum(FC(1:pp))+sum( FU*c )*c)-FR
  
      elseif (BC_type .eq.  -1) then
   !BC reflecting
        FC(:)=F(K,1,:)
        FU(:)=F(NV-K+1,1,:)
        FR(:)=FS(K,1,:)
!        F_s(K,1,:,1)=( (-FR)' .* b)
!        F_ns(K,1,:,1)=( (V(K)*A'*FC -V(K)* sum(FC)+ V(K)* sum(FU'.* c) * c')' .* b)
      F_ns(K,1,:,l)=V(K)*(tempAu(:)-sum(FC(1:pp))+sum( FU*c )*c)-FR
      else
     write(6,*) "BC_type not supported",BC_type
     stop
   endif
    
  do i=2,nnx
  tempAu=matmul(F(K,i,1:pp),A)
      FU(:)=F(K,i-1,:)
      FC(:)=F(K,i,:)
      FR(:)=FS(K,i,:)
      F_ns(K,i,:,l)=V(K)*(tempAu(:)-sum(FC(1:pp))+sum(FU)*c)-FR

!     do j=1,pp
!        F_ns(i,j,l)=tempAu(j)-sum(F(i,:))+sum(F(i-1,:))*c(j)+FS(i,j)
!     enddo
!     do j=1,pp
!        F_ns(i,j,l)=F_ns(i,j,l)*b(j)
!     enddo
  enddo
! Left-going part
    
  do i=1,nnx-1
  tempAu=matmul(F(NVh+K,i,1:pp),A)
      FU(:)=F(NVh+K,i+1,:)
      FC(:)=F(NVh+K,i,:)
      FR(:)=FS(NVh+K,i,:)
      F_ns(NVh+K,i,:,l)=V(NVh+K)*(tempAu(:)-sum(FU(1:pp)*c)+sum( FC*c )*c)-FR
  enddo

  tempAu=matmul(F(NVh+K,nnx,1:pp),A)
 if (BC_type .eq. 0) then
  !BC no-flux
        FC(:)=F(NVh+K,nnx,:)
        FU=FC
        FR(:)=FS(NVh+K,nnx,:)
      F_ns(NVh+K,nnx,:,l)=V(NVh+K)*(tempAu(:)-sum(FU(1:pp))+sum( FC*c )*c)-FR
  
      elseif (BC_type .eq.  -1) then
   !BC reflecting
        FC(:)=F(NVh+K,nnx,:)
        FU(:)=F(NV-K+1,nnx,:)
        FR(:)=FS(NVh+K,nnx,:)
      F_ns(NVh+K,nnx,:,l)=V(NVh+K)*(tempAu(:)-sum(FU(1:pp))+sum( FC*c )*c)-FR
      else
     write(6,*) "BC_type not supported",BC_type
     stop
   endif

enddo !K 


     do j=1,pp
        F_ns(:,:,j,l)=F_ns(:,:,j,l)*b(j)
     enddo


  end subroutine DG_flux_and_bdd
  !------------------------------------------------------------------------------


  real(kind=8) function func0(x)
  implicit none
  real(kind=8)::x

 if ((x>0.4d0) .and.(x<0.6d0)) then
   func0=(x-0.4d0)**10*(x-0.6d0)**10*10d0**20
 else
   func0=0d0
 end if

  return
  end function

logical function isnan(x)
real x !(kind=8) x
isnan = x .ne. x
end
