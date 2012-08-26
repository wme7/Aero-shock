      subroutine pordr
!     ***************
      use mod_wave
      implicit none
      double precision::epl1,epli,esl1,esli
      double precision::epl2,esl2
      double precision::ds,dp,alfam,aux,axxm
      integer::i
!      --- control
      if (ne.eq.1) return
            pordl2(1)=0.d0
            pordl1(1)=0.d0
            pordli(1)=0.d0
!                       
!            --- initialization of the iterations on ne
            esl2=perrl2(1)
            esl1=perrl1(1)
            esli=perrli(1)
!
            ds=ddx(1)
!      
!            --- iterations on ne
            do i=2,ne
               epl2=esl2
               epl1=esl1
               epli=esli
!          
               dp=ds
!          
               esl2=perrl2(i)
               esl1=perrl1(i)
               esli=perrli(i)
!
               ds=ddx(i)
!

!               --- computing the L1-order
               if (epl1.gt.0..and.esl1.gt.0..and.ds.ne.dp) then
                  pordl1(i)=dlog(epl1/esl1)/dlog(dp/ds)
               else
                  pordl1(i)=0.d0
               endif
!               --- computing the L2-order
               if (epl2.gt.0..and.esl2.gt.0..and.ds.ne.dp) then
                  pordl2(i)=dlog(epl2/esl2)/dlog(dp/ds)
               else
                  pordl2(i)=0.d0
               endif
!               write(6,*) esl1,pordl1(i)
!               --- computing the Li-order
               if (epli.gt.0..and.esli.gt.0..and.ds.ne.dp) then
                  pordli(i)=dlog(epli/esli)/dlog(dp/ds)
               else
                  pordli(i)=0.d0
               endif
             enddo
!   11          continue
! 
!         --- computation of the constants
!            --- iterations on ne
            do  i=1,ne
               ds=ddx(i) 
!
!               --- computing the L1-constant
               alfam=pordl1(i)
               aux=alfam*dlog(ds)
               if (dabs(aux).gt.10) then
                  axxm=1.
               else
                  axxm=dexp(aux)
               endif
!
               constl1(i)=perrl1(i)/axxm   
!
!               --- computing the Li-constant
               alfam=pordli(i)
               aux=alfam*dlog(ds)
               if (dabs(aux).gt.10) then
                  axxm=1.
               else
                  axxm=dexp(aux)
               endif
!
               constli(i)=perrli(i)/axxm   
!
             enddo
!   12          continue
!   10    continue
!
      return
      end subroutine pordr
