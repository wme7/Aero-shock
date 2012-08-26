      subroutine ordr
c     ***************
c
      include 'com'
c
c      --- control
      if (ne.eq.1) return
c
c         --- computation of the orders
            ordl1(1)=0.
            ordli(1)=0.
c                       
c            --- initialization of the iterations on ne
            esl1=errl1(1)
            esli=errli(1)
c
            ds=ddx(1)
c      
c            --- iterations on ne
            do 11 i=2,ne
               epl1=esl1
               epli=esli
c          
               dp=ds
c          
               esl1=errl1(i)
               esli=errli(i)
c
               ds=ddx(i)
c
c               --- computing the L1-order
               if (epl1.gt.0..and.esl1.gt.0..and.ds.ne.dp) then
                  ordl1(i)=dlog(epl1/esl1)/dlog(dp/ds)
               else
                  ordl1(i)=0.
               endif
c
c               --- computing the Li-order
               if (epli.gt.0..and.esli.gt.0..and.ds.ne.dp) then
                  ordli(i)=dlog(epli/esli)/dlog(dp/ds)
               else
                  ordli(i)=0.
               endif
c
   11          continue
c 
c         --- computation of the constants
c
c            --- iterations on ne
            do 12 i=1,ne
               ds=ddx(i) 
c
c               --- computing the L1-constant
               alfam=ordl1(i)
               aux=alfam*dlog(ds)
               if (dabs(aux).gt.10) then
                  axxm=1.
               else
                  axxm=dexp(aux)
               endif
c
               constl1(i)=errl1(i)/axxm   
c
c               --- computing the Li-constant
               alfam=ordli(i)
               aux=alfam*dlog(ds)
               if (dabs(aux).gt.10) then
                  axxm=1.
               else
                  axxm=dexp(aux)
               endif
c
               constli(i)=errli(i)/axxm   
c
   12       continue
c
      return
      end
