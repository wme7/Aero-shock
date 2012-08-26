      subroutine bdryflux(i,fu,fv,tn,j)
      use mod_wave
      implicit none
!
      double precision, external:: func
      integer, intent(in)::i,j
      double precision, intent(in)::tn
      double precision,intent(out)::fu,fv
      double precision:: ffunc(kt)
!      double precision::AA(kt)
      double precision:: s,bc,aux
      double precision::utrplus,vtrplus,utrmins,vtrmins
      integer::l,m,jj,tst 
! 
!     flux at x=-a 
!     ------------ 
         tst=1 
      if (i.eq.1) then 
 
! 
!      evaluation of the traces 
!      ------------------------ 
      utrplus=0.d00 
      vtrplus=0.d00 
! 
      do l=0,k 
       utrplus=utrplus+u(1,l+1,j)*(-1)**l 
       vtrplus=vtrplus+v(1,l+1,j)*(-1)**l 
      enddo 
!      evaluation of the fluxes 
!      ------------------------ 
! 
!********************************************************* 
!       --- the boundary condition 
!           s=-a-tn-j*dt 
!       bc=func(s) 
       s=-a-tn 
       if (nbc .gt. 0) then 
        if (nbc .eq. 1) then 
!         write(6,*) "scheme 1" 
 
!      computing the moments of the initial data in unc' 
!      and storing them in hi' 
! 
!        evaluating the function unc' at the quadrature points 
         if (tst .eq. 0) then 
          do m=1,kt 
            ffunc(m)=func(-a-(tn-dt*(xt(m)-1.d0)/2.d0)) 
!            ffunc(m)=ffunc(m)+(tn-dt*(xt(m)-1.d0)/2.d0)**4 
          enddo 
          do jj=0,kt-1 
            AA(jj+1)=0.d00 
           do m=1,kt 
            AA(jj+1)= AA(jj+1)+ffunc(m)*plt(m,jj+1)*wt(m) 
           enddo 
            AA(jj+1)= AA(jj+1)*(2.d00*jj+1.d00)/2.d00 
          end do 
         end if 
 
      bc=0.d0 
!        write(6,*) AA 
!        write(6,*) dl 
          do m=0,j 
          do jj=0,kt-1 
           bc=bc+cbin(j+1,m+1)*(-2.d0)**m*AA(jj+1)*dl(jj+1,m+1) 
        end do 
          end do  
        endif 
 
 
        if (nbc .eq. 2) then 
!         write(6,*) "scheme 2" 
 
          bc=f(1) 
 
!		  if (f(1) .gt. .95) write(6,*) f 
          do m=1,j 
		  if (j .eq. 0) write(6,*) "alert j=",j 
          aux=AA(m) 
          do jj=m+1,kt-1 
            aux=aux+AA(jj)*dt**(jj-m)*dl(jj-1,jj-m) 
          end do 
          bc=bc+cbin(j+1,m+1)*aux*(dt**m)*dl(m,m) 
 
          end do  
        endif 
 
 
       endif 
!             
 
      fu=-vtrplus 
      fv=-2*bc-vtrplus 
!       write(6,*) fu,fv 
!      fu=bc 
!      fv=-bc  
!      else 
!      fu=-vtrplus 
!      fv=-utrplus 
!      endif 
 
!      
      else 
! 
!     flux at x=+a 
!     ------------ 
! 
!      evaluation of the traces 
!      ------------------------ 
      utrmins=0.d00 
      vtrmins=0.d00 
! 
      do l=0,k 
       utrmins=utrmins+u(nx,l+1,j) 
       vtrmins=vtrmins+v(nx,l+1,j) 
      enddo 
! 
!      evaluation of the fluxes 
!      ------------------------ 
      fu=-vtrmins 
      fv=-fu 
!      fv=-utrmins 
!   
      endif 
! 
!      deallocate(ffunc) 
!      write(6,*) cbin 
      return 
   
      end subroutine bdryflux 
