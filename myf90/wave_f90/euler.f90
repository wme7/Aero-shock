      subroutine euler(tn,j)
!     **********************
!       compute u(.,.,j+1) from u(.,.,j)
!       compute v(.,.,j+1) from v(.,.,j)
      use mod_wave
      implicit none
      double precision,intent(in)::tn
      integer,intent(in)::j
      double precision ffunc(kt)
!      double precision:: AA(kt+1)
      double precision::bc,aux,fu,fv
      integer::m,jj,i,l

!      write(6,*) "euler"
!          do m=1,kt
!            ffunc(m)=4.d0*(tn-dt*(xt(m)-1.d0)/2.d0)**3
!          enddo
!          do jj=0,kt-1
!            AA(jj+1)=0.d00
!           do m=1,kt
!            AA(jj+1)= AA(jj+1)+ffunc(m)*plt(m,jj+1)*wt(m)
!           enddo
!            AA(jj+1)= AA(jj+1)*(2.d00*jj+1.d00)/2.d00
!          end do
      bc=0.d0   
!          do m=0,j
!           do jj=0,k
!           bc=bc+cbin(j+1,m+1)*(-2.d0)**m*AA(jj+1)*dl(jj+1,m+1)
!           end do
!          end do

           bc=0.d0


!      --- evaluating the fluxes
      do 11 i=1,nx-1
        call intflux(i,fu,fv,u,v,j)
        fluxu(i)=fu
        fluxv(i)=fv
!       write(6,*) nproblem
   11 continue
        if ((tn .le.1) .or. (nproblem .le. 2)) then
        call bdryflux(1,fu,fv,tn,j)
        fluxu(0)=fu
        fluxv(0)=fv
        call bdryflux(2,fu,fv,tn,j)
        fluxu(nx)=fu
        fluxv(nx)=fv
        else
!        call pbdryflux(fu,fv,u,v,j)
!        fluxu(nx)=fu
!        fluxv(nx)=fv
!        fluxu(0)=fu
!        fluxv(0)=fv
!        write(6,*) fu,fv
        end if
!     
!      --- computing uo at t=tn1
!      --- the interior values
      do i=1,nx
         do l=0,k
!         computing u
         aux=0.d00
         do m=0,l-1
            aux=aux+v(i,m+1,j)*(1-(-1)**(m+l))
         enddo
         u(i,l+1,j+1)=u(i,l+1,j)+c(i)*c(i)*dt/dx*dble(2*l+1)*  &
         (-aux - fluxu(i)+fluxu(i-1)*(-1)**l)
!    
!         computing v
         aux=0.d00
         do m=0,l-1
            aux=aux+u(i,m+1,j)*(1-(-1)**(m+l))
         enddo
         v(i,l+1,j+1)=v(i,l+1,j)+dt/dx*dble(2*l+1)*  &
         (-aux - fluxv(i)+fluxv(i-1)*(-1)**l)
         enddo
         u(i,1,j+1)=u(i,1,j+1)+dt*bc
         v(i,1,j+1)=v(i,1,j+1)-dt*bc
        enddo
!             
      end subroutine euler
