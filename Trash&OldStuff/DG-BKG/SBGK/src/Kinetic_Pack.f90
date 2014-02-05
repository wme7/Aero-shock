subroutine Compute_Source(nx,pp,dx)
use Legendre
use State_Var
use universal_const

IMPLICIT NONE
integer:: nx,pp,p
real(kind=8)::dx

integer:: i,K,j
real(kind=8)::FC(1:pp),FB(1:pp)

p=pp-1

      do i=1,nx
        do K = 1,NV

         FC(:)=F(K,i,:)
         FB(:)=FEQ(K,i,:)
         FC=matmul(FC,Pleg)-FB
          do j=0,p
           FS(K,i,j+1)=sum(FC*Pleg(j+1,:)*wl)*dx/2d0/VIS(i,j+1)
          enddo
        enddo               
      enddo


end subroutine Compute_Source

SUBROUTINE EQUILIBRIUM(nx,pp) 
use Legendre
use State_Var
use universal_const

IMPLICIT NONE
INTEGER:: i,K,m
integer:: nx,pp

   do i=1,nx
    do K=1,NV
     do m=1,pp
    FEQ(K,i,m)   = 1d0/((exp( (V(K)-Us(i,m))**2 /Ts(i,m))/Zs(i,m)) + dble(IT) )
    if (FEQ(K,i,m) .gt. 1d20) then
    write(6,*) "Infty",i,K,m
    stop
    endif
     END DO
  END DO
END DO
end SUBROUTINE EQUILIBRIUM

subroutine  Comp_RU(nnx,pp)
use State_Var
use Legendre
implicit none
integer:: nnx,pp
integer:: i,m,K
real(kind=8):: Mtemp(NV,pp)
real(kind=8):: F_loc(1:NV,1:pp)
real(kind=8):: SR(1:nnx,1:pp),SU(1:nnx,1:pp),SE(1:nnx,1:pp),SAV(1:nnx,1:pp)


do i=1,nnx

    Mtemp=0d0
    do K=1,NV
        Mtemp(K,:)=F(K,i,:)
    enddo

    F_loc(:,:)=matmul(Mtemp,Pleg)

    do m=1,pp
        SR(i,:) = matmul(VW,F_loc)
        SU(i,m) = sum( VW*F_loc(:,m)* V)
        SE(i,m) = sum( VW*F_loc(:,m)* V**2)/2d0
        SAV(i,m)= sum( VW*F_loc(:,m)* abs(V))
        
        Rs(i,m)    = SR(i,m)
        Us(i,m)    = SU(i,m)/SR(i,m)
        ET(i,m)   = SE(i,m)
        AV(i,m)   = SAV(i,m)
    enddo 
enddo

end subroutine Comp_RU

subroutine Comp_ZTP(nnx,pp)
use Legendre
use State_Var
use universal_const
IMPLICIT NONE
integer,intent(in):: nnx,pp
integer:: i,m,K,L
real(kind=8):: ZA, ZB, ZC, GA12, GB12, GA32, GB32, GC12, GC32, GC52
real(kind=8):: PSIA, PSIB, PSIC

if (IT .eq. 0) then
    do i=1,nnx
        do m=1,pp
            Ts(i,m)    = 4d0*ET(i,m)/Rs(i,m) - 2d0*Us(i,m)**2
            Zs(i,m)    = Rs(i,m) / sqrt(pi* Ts(i,m))
            Ps(i,m) = ET(i,m) - 0.5d0 * Rs(i,m) * Us(i,m)**2
    if (Ts(i,m) <0) then
       write(6,*) "T is Negative", i, m
       stop
    endif
        enddo
    enddo
else
    do i=1,nnx
        do m=1,pp

            ZA = 0.0001d0
            ZB = 0.99d0
            do while (abs(ZA-ZB) > 1d-10)
                GA12 = 0d0
                GB12 = 0d0
                GA32 = 0d0
                GB32 = 0d0
                do L = 1,50
                    if (IT .eq. 1) then
                        GA12 = GA12 + (ZA**L)*(-1)**(L-1)/dsqrt(dble(L))
                        GB12 = GB12 + (ZB**L)*(-1)**(L-1)/dsqrt(dble(L))
                        GA32 = GA32 + (ZA**L)*(-1)**(L-1)/(dble(L)*dsqrt(dble(L)))
                        GB32 = GB32 + (ZB**L)*(-1)**(L-1)/(dble(L)*dsqrt(dble(L)))
                    else
                        GA12 = GA12 + (ZA**L)/dsqrt(dble(L))
                        GB12 = GB12 + (ZB**L)/dsqrt(dble(L))
                        GA32 = GA32 + (ZA**L)/(dble(L)*dsqrt(dble(L)))
                        GB32 = GB32 + (ZB**L)/(dble(L)*dsqrt(dble(L)))
                    endif
                enddo
                PSIA = 2d0*ET(i,m) - GA32*(Rs(i,m)/GA12)**3/(2d0*pi) - Rs(i,m)*Us(i,m)**2
                PSIB = 2d0*ET(i,m) - GB32*(Rs(i,m)/GB12)**3/(2d0*pi) - Rs(i,m)*Us(i,m)**2
                ZC = (ZA + ZB)/2d0
                GC12 = 0d0
                GC32 = 0d0
                GC52 = 0d0
                do L = 1,50
                    if  (IT .eq. 1) then
                        GC12 = GC12 + (ZC**L)*(-1)**(L-1)/dsqrt(dble(L))
                        GC32 = GC32 + (ZC**L)*(-1)**(L-1)/(dble(L)*dsqrt(dble(L)))
                        GC52 = GC52 + (ZC**L)*(-1)**(L-1)/(dble(L**2)*dsqrt(dble(L)))
                    else
                        GC12 = GC12 + (ZC**L)/dsqrt(dble(L))
                        GC32 = GC32 + (ZC**L)/(dble(L)*dsqrt(dble(L)))
                        GC52 = GC52 + (ZC**L)/(dble(L**2)*dsqrt(dble(L)))
                    endif
                enddo
                PSIC = 2d0*ET(i,m) - GC32*(Rs(i,m)/GC12)**3/(2d0*pi) - Rs(i,m)*Us(i,m)**2
                
                if ((PSIA*PSIC) < 0d0) then
                    ZB = ZC
                else
                    ZA = ZC
                endif
            enddo
            Zs(i,m) = ZC
            Ts(i,m) = Rs(i,m)**2 / (pi*GC12**2 )
            Ps(i,m) = ET(i,m) - 0.5d0 * Rs(i,m) * Us(i,m)**2
            
        enddo
      enddo
endif !if IT
end subroutine Comp_ZTP

subroutine Comp_RU_F(nnx,pp,Fv)
use State_Var, ONLY: Rs, Us, ET, AV,NV, VW, V
use Legendre, ONLY: Pleg
implicit none
integer,intent(in):: nnx,pp
integer:: i,m,K
real(kind=8),intent(in):: Fv(NV,nnx,pp)
real(kind=8):: Mtemp(NV,pp)
real(kind=8):: F_loc(1:NV,1:pp)
real(kind=8):: SR(1:nnx,1:pp),SU(1:nnx,1:pp),SE(1:nnx,1:pp),SAV(1:nnx,1:pp)
do i=1,nnx

    Mtemp=0d0
    do K=1,NV
        Mtemp(K,:)=Fv(K,i,:)
    enddo

    F_loc(:,:)=matmul(Mtemp,Pleg)

    do m=1,pp
        SR(i,:) = matmul(VW,F_loc)
        SU(i,m) = sum( VW*F_loc(:,m)* V)
        SE(i,m) = sum( VW*F_loc(:,m)* V**2)/2d0
        SAV(i,m)= sum( VW*F_loc(:,m)* abs(V))
        
        Rs(i,m)    = SR(i,m)
        Us(i,m)    = SU(i,m)/SR(i,m)
        ET(i,m)   = SE(i,m)
        AV(i,m)   = SAV(i,m)
    enddo 
enddo

end subroutine Comp_RU_F

subroutine Comp_ZTP_F(nnx,pp, IT, Rloc, Uloc, ETloc, Zloc, Tloc)

!use Legendre
!use State_Var, ONLY: Ts, Zs, Ps, ET, Rs, Us, IT
use universal_const
IMPLICIT NONE
integer, intent(in):: nnx,pp, IT
integer:: i,m,K,L
real(kind=8):: ZA, ZB, ZC, GA12, GB12, GA32, GB32, GC12, GC32, GC52
real(kind=8):: PSIA, PSIB, PSIC
real(kind=8), intent(in) :: Rloc(1:nnx,1:pp), Uloc(1:nnx,1:pp)
real(kind=8), intent(out):: Zloc(1:nnx,1:pp), Tloc(1:nnx,1:pp)
real(kind=8), intent(in) :: ETloc(1:nnx,1:pp)

if (IT .eq. 0) then
    do i=1,nnx
        do m=1,pp
            Tloc(i,m)    = 4d0*ETloc(i,m)/Rloc(i,m) - 2d0*Uloc(i,m)**2
            Zloc(i,m)    = Rloc(i,m) / sqrt(pi* Tloc(i,m))
!            Ps(i,m) = ETloc(i,m) - 0.5d0 * Rloc(i,m) * Uloc(i,m)**2
        enddo
    enddo
else
    do i=1,nnx
        do m=1,pp

            ZA = 0.0001d0
            ZB = 0.99d0
            do while (abs(ZA-ZB) > 1d-10)
                GA12 = 0d0
                GB12 = 0d0
                GA32 = 0d0
                GB32 = 0d0
                do L = 1,50
                    if (IT .eq. 1) then
                        GA12 = GA12 + (ZA**L)*(-1)**(L-1)/dsqrt(dble(L))
                        GB12 = GB12 + (ZB**L)*(-1)**(L-1)/dsqrt(dble(L))
                        GA32 = GA32 + (ZA**L)*(-1)**(L-1)/(dble(L)*dsqrt(dble(L)))
                        GB32 = GB32 + (ZB**L)*(-1)**(L-1)/(dble(L)*dsqrt(dble(L)))
                    else
                        GA12 = GA12 + (ZA**L)/dsqrt(dble(L))
                        GB12 = GB12 + (ZB**L)/dsqrt(dble(L))
                        GA32 = GA32 + (ZA**L)/(dble(L)*dsqrt(dble(L)))
                        GB32 = GB32 + (ZB**L)/(dble(L)*dsqrt(dble(L)))
                    endif
                enddo
                PSIA = 2d0*ETloc(i,m) - GA32*(Rloc(i,m)/GA12)**3/(2d0*pi) - Rloc(i,m)*Uloc(i,m)**2
                PSIB = 2d0*ETloc(i,m) - GB32*(Rloc(i,m)/GB12)**3/(2d0*pi) - Rloc(i,m)*Uloc(i,m)**2
                ZC = (ZA + ZB)/2d0
                GC12 = 0d0
                GC32 = 0d0
                GC52 = 0d0
                do L = 1,50
                    if  (IT .eq. 1) then
                        GC12 = GC12 + (ZC**L)*(-1)**(L-1)/dsqrt(dble(L))
                        GC32 = GC32 + (ZC**L)*(-1)**(L-1)/(dble(L)*dsqrt(dble(L)))
                        GC52 = GC52 + (ZC**L)*(-1)**(L-1)/(dble(L**2)*dsqrt(dble(L)))
                    else
                        GC12 = GC12 + (ZC**L)/dsqrt(dble(L))
                        GC32 = GC32 + (ZC**L)/(dble(L)*dsqrt(dble(L)))
                        GC52 = GC52 + (ZC**L)/(dble(L**2)*dsqrt(dble(L)))
                    endif
                enddo
                PSIC = 2d0*ETloc(i,m) - GC32*(Rloc(i,m)/GC12)**3/(2d0*pi) - Rloc(i,m)*Uloc(i,m)**2
                
                if ((PSIA*PSIC) < 0d0) then
                    ZB = ZC
                else
                    ZA = ZC
                endif
            enddo
            Zloc(i,m) = ZC
            Tloc(i,m) = Rloc(i,m)**2 / (pi*GC12**2 )
!            Ps(i,m) = ETloc(i,m) - 0.5d0 * Rloc(i,m) * Uloc(i,m)**2
            
        enddo
      enddo
endif !if IT
end subroutine Comp_ZTP_F

subroutine Compute_Source_F(nx,pp,dx,Fv,Fso)
use Legendre, ONLY: Pleg, wl
use State_Var, ONLY: FEQ, VIS, NV
use universal_const

IMPLICIT NONE
integer, intent(in):: nx,pp

real(kind=8), intent(in):: Fv(NV,nx,pp)
real(kind=8), intent(out):: Fso(NV,nx,pp)
real(kind=8), intent(in) ::dx

integer:: i,K,j,p
real(kind=8)::FC(1:pp),FB(1:pp)

p=pp-1

      do i=1,nx
        do K = 1,NV

         FC(:)=Fv(K,i,:)
         FB(:)=FEQ(K,i,:)
         FC=matmul(FC,Pleg)-FB
          do j=0,p
           Fso(K,i,j+1)=sum(FC*Pleg(j+1,:)*wl)*dx/2d0/VIS(i,j+1)
          enddo
        enddo               
      enddo


end subroutine Compute_Source_F

SUBROUTINE EQUILIBRIUM_F(nx,pp) 
!use Legendre
use State_Var, ONLY: FEQ,V,Us,Ts,Zs,IT,NV
use universal_const

IMPLICIT NONE
integer, intent(in):: nx,pp

INTEGER:: i,K,m

   do i=1,nx
    do K=1,NV
     do m=1,pp
    FEQ(K,i,m)   = 1d0/((exp( (V(K)-Us(i,m))**2 /Ts(i,m))/Zs(i,m)) + dble(IT) )
    if (FEQ(K,i,m) .gt. 1d20) then
    write(6,*) "Infty",i,K,m
    stop
    endif
     END DO
  END DO
END DO
end SUBROUTINE EQUILIBRIUM_F

subroutine Comp_ZTP_S(nnx,pp, IT, Rloc, Uloc, ETloc, Zloc, Tloc)

!use Legendre
!use State_Var, ONLY: Ts, Zs, Ps, ET, Rs, Us, IT
use universal_const
IMPLICIT NONE
integer, intent(in):: nnx,pp, IT
integer:: m,K,L
real(kind=8):: ZA, ZB, ZC, GA12, GB12, GA32, GB32, GC12, GC32, GC52
real(kind=8):: PSIA, PSIB, PSIC
real(kind=8), intent(in) :: Rloc(1:pp), Uloc(1:pp)
real(kind=8), intent(out):: Zloc(1:pp), Tloc(1:pp)
real(kind=8), intent(in) :: ETloc(1:pp)

if (IT .eq. 0) then
        do m=1,pp
            Tloc(m)    = 4d0*ETloc(m)/Rloc(m) - 2d0*Uloc(m)**2
            Zloc(m)    = Rloc(m) / sqrt(pi* Tloc(m))
!            Ps(m) = ETloc(m) - 0.5d0 * Rloc(m) * Uloc(m)**2
        enddo
else
        do m=1,pp

            ZA = 0.0001d0
            ZB = 0.99d0
            do while (abs(ZA-ZB) > 1d-10)
                GA12 = 0d0
                GB12 = 0d0
                GA32 = 0d0
                GB32 = 0d0
                do L = 1,50
                    if (IT .eq. 1) then
                        GA12 = GA12 + (ZA**L)*(-1)**(L-1)/dsqrt(dble(L))
                        GB12 = GB12 + (ZB**L)*(-1)**(L-1)/dsqrt(dble(L))
                        GA32 = GA32 + (ZA**L)*(-1)**(L-1)/(dble(L)*dsqrt(dble(L)))
                        GB32 = GB32 + (ZB**L)*(-1)**(L-1)/(dble(L)*dsqrt(dble(L)))
                    else
                        GA12 = GA12 + (ZA**L)/dsqrt(dble(L))
                        GB12 = GB12 + (ZB**L)/dsqrt(dble(L))
                        GA32 = GA32 + (ZA**L)/(dble(L)*dsqrt(dble(L)))
                        GB32 = GB32 + (ZB**L)/(dble(L)*dsqrt(dble(L)))
                    endif
                enddo
                PSIA = 2d0*ETloc(m) - GA32*(Rloc(m)/GA12)**3/(2d0*pi) - Rloc(m)*Uloc(m)**2
                PSIB = 2d0*ETloc(m) - GB32*(Rloc(m)/GB12)**3/(2d0*pi) - Rloc(m)*Uloc(m)**2
                ZC = (ZA + ZB)/2d0
                GC12 = 0d0
                GC32 = 0d0
                GC52 = 0d0
                do L = 1,50
                    if  (IT .eq. 1) then
                        GC12 = GC12 + (ZC**L)*(-1)**(L-1)/dsqrt(dble(L))
                        GC32 = GC32 + (ZC**L)*(-1)**(L-1)/(dble(L)*dsqrt(dble(L)))
                        GC52 = GC52 + (ZC**L)*(-1)**(L-1)/(dble(L**2)*dsqrt(dble(L)))
                    else
                        GC12 = GC12 + (ZC**L)/dsqrt(dble(L))
                        GC32 = GC32 + (ZC**L)/(dble(L)*dsqrt(dble(L)))
                        GC52 = GC52 + (ZC**L)/(dble(L**2)*dsqrt(dble(L)))
                    endif
                enddo
                PSIC = 2d0*ETloc(m) - GC32*(Rloc(m)/GC12)**3/(2d0*pi) - Rloc(m)*Uloc(m)**2
                
                if ((PSIA*PSIC) < 0d0) then
                    ZB = ZC
                else
                    ZA = ZC
                endif
            enddo
            Zloc(m) = ZC
            Tloc(m) = Rloc(m)**2 / (pi*GC12**2 )
!            Ps(m) = ETloc(m) - 0.5d0 * Rloc(m) * Uloc(m)**2
            
        enddo
endif !if IT
end subroutine Comp_ZTP_S
