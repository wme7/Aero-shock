      subroutine intflux(i,fu,uu,j)
c     *****************************
c
      include 'com'
c
      double precision uu(NINT,KK+1,0:KK+1)
c
c      evaluation of the traces
c      ------------------------
      utrplus=0.d00
      utrmins=0.d00
c
      do l=0,k
       utrmins=utrmins+uu(i,l+1,j)
      enddo
c
c      evaluation of the fluxes
c      ------------------------
      fu=utrmins
c
      return
      end
