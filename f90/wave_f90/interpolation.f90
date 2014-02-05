      subroutine interpolation(tn)
!     ***************
!     interpolation           
!     ***************
      use mod_wave
      implicit none
      double precision, intent(in)::tn
      double precision,external::func
      double precision::s
      double precision,dimension(kt)::ffunc
       integer, dimension(ks)::IA
      double precision,dimension(kt,2)::BB
      integer::i,j,m,il,l,ll,lo,in,jj


!       --- the boundary condition
       s=-a-tn
       if (nbc .gt. 0) then
        if (nbc .eq. 1) then
!         write(6,*) "scheme 1"

!      computing the moments of the initial data in unc'
!      and storing them in hi'
!
!        evaluating the function unc' at the quadrature points
          do m=1,kt
            ffunc(m)=func(-a-(tn-dt*(xt(m)-1.d0)/2.d0))
          enddo
!
          do jj=0,kt-1
            AA(jj+1)=0.d00
           do m=1,kt
            AA(jj+1)= AA(jj+1)+ffunc(m)*plt(m,jj+1)*wt(m)
           enddo
            AA(jj+1)= AA(jj+1)*(2.d00*jj+1.d00)/2.d00
          end do
        endif
       endif
        if (nbc .eq. 2) then

! *************************************
! ENO       
!**************************************
!          do m=1,kt
!		  f(m)=m**(kt-1)
!		  enddo


!		  dt=1.d0
          AA(1)=(f(2)-f(1))/dt
           do l=1,kt-1
            BB(l,1)=-(f(l+1)-f(l))/dt
           enddo
           AA(1)=BB(1,1)
           
           lo=1
           ll=2
           do m=1,kt-2
           do l=1,kt-1-m
           BB(l,ll)=-(BB(l+1,lo)-BB(l,lo))/dt/(dble(m+1))
		  
           enddo
!		   write(6,*) (BB(jj,ll),jj=1,kt-1)
           AA(m+1)=BB(1,ll)
           l=lo
           lo=ll
           ll=l
           enddo

        endif
!       write(6,*) f
!       write(6,*) AA
!	   pause
      return
      end subroutine interpolation
