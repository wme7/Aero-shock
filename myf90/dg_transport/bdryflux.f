      subroutine bdryflux(i,fu,tn,j)
c     ******************************
c
      include 'com'
c
c     flux at x=0
c     ------------
      if (i.eq.1) then
c
c      evaluation of the traces
c      ------------------------
      utrmins=0.d00
c
      do l=0,k
       utrmins=utrmins+u(nx,l+1,j)
      enddo
c            
      fu=utrmins
c     
      else
c
c     flux at x=rl
c     ------------
c
c      evaluation of the traces
c      ------------------------
      utrmins=0.d00
c
      do l=0,k
       utrmins=utrmins+u(nx,l+1,j)
      enddo
c
c      evaluation of the fluxes
c      ------------------------
      fu=utrmins
c  
      endif
c
      return
      end
