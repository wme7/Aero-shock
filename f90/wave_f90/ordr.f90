      subroutine ordr
!     ***************
      use mod_wave
      implicit none
      double precision::epl1,epli,esl1,esli,epl2,esl2
      double precision::ds,dp,alfam,aux,axxm
      integer::i
!      --- control
      if (ne.eq.1) return
!      --- if indi=0 set indinierr=1
!      if (indi.eq.0) then
!         indinierr=1
!      else
!         indinierr=nierr
!         endif
!      --- loop on the 'error intervals'
!      do 10 ierr=1,indinierr 
!         --- computation of the orders
            ordl2(1)=0.d0
            ordl1(1)=0.d0
            ordli(1)=0.d0
!                       
!            --- initialization of the iterations on ne
            esl2=errl2(1)
            esl1=errl1(1)
            esli=errli(1)
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
               esl2=errl2(i)
               esl1=errl1(i)
               esli=errli(i)
!
               ds=ddx(i)
!
!               --- computing the L1-order
               if (epl1.gt.0..and.esl1.gt.0..and.ds.ne.dp) then
                  ordl1(i)=dlog(epl1/esl1)/dlog(dp/ds)
               else
                  ordl1(i)=0.d0
               endif
!               --- computing the L2-order
               if (epl2.gt.0..and.esl2.gt.0..and.ds.ne.dp) then
                  ordl2(i)=dlog(epl2/esl2)/dlog(dp/ds)
               else
                  ordl2(i)=0.d0
               endif
!               --- computing the Li-order
               if (epli.gt.0..and.esli.gt.0..and.ds.ne.dp) then
                  ordli(i)=dlog(epli/esli)/dlog(dp/ds)
               else
                  ordli(i)=0.d0
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
               alfam=ordl1(i)
               aux=alfam*dlog(ds)
               if (dabs(aux).gt.10) then
                  axxm=1.
               else
                  axxm=dexp(aux)
               endif
!
               constl1(i)=errl1(i)/axxm   
!
!               --- computing the Li-constant
               alfam=ordli(i)
               aux=alfam*dlog(ds)
               if (dabs(aux).gt.10) then
                  axxm=1.
               else
                  axxm=dexp(aux)
               endif
!
               constli(i)=errli(i)/axxm   
!
             enddo
!   12          continue
!   10    continue
!
      return
      end subroutine ordr
