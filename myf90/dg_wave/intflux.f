      subroutine intflux(i,fu,fv,uu,vv,j)
c     ***********************************
c
      include 'com'
c
      double precision uu(NINT,KK+1,0:KK+1),vv(NINT,KK+1,0:KK+1)
c
c      evaluation of the traces
c      ------------------------
      utrplus=0.d00
      utrmins=0.d00
c
      vtrplus=0.d00
      vtrmins=0.d00
c
      do l=0,k
       utrplus=utrplus+(-1)**l*uu(i+1,l+1,j)
       utrmins=utrmins+uu(i,l+1,j)
c
       vtrplus=vtrplus+(-1)**l*vv(i+1,l+1,j)
       vtrmins=vtrmins+vv(i,l+1,j)
      enddo
c
c      evaluation of the fluxes
c      ------------------------
      uave=(utrplus+utrmins)/2.0d00
      vave=(vtrplus+vtrmins)/2.0d00
c
      ujmp=utrplus-utrmins
      vjmp=vtrplus-vtrmins
c
      vhat=vave+c11(i)*ujmp-c12(i)*vjmp
      uhat=uave+c12(i)*ujmp+c22(i)*vjmp
ccccc      vhat=vave
ccccc      uhat=uave
c
      fu=-vhat
      fv=-uhat
c
      return
      end
