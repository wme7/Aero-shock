      subroutine bdryflux(i,fu,fv,uu,vv,tn,j)
c     ***************************************
c
      include 'com'
c
      double precision uu(NINT,KK+1,0:KK+1),vv(NINT,KK+1,0:KK+1)
c
c     flux at x=-a
c     ------------
      if (i.eq.1) then
c
c      evaluation of the traces
c      ------------------------
      utrplus=0.d00
      vtrplus=0.d00
c
      do l=0,k
       utrplus=utrplus+u(1,l+1,j)*(-1)**l
       vtrplus=vtrplus+v(1,l+1,j)*(-1)**l
      enddo
c      evaluation of the fluxes
c      ------------------------
c
c       --- the boundary condition
       bc=func(-a-tn-j*dt)
c
       fu=-vtrplus
       fv=-2*bc-vtrplus
c     
      else
c
c     flux at x=+a
c     ------------
c
c      evaluation of the traces
c      ------------------------
      utrmins=0.d00
      vtrmins=0.d00
c
      do l=0,k
       utrmins=utrmins+u(nx,l+1,j)
       vtrmins=vtrmins+v(nx,l+1,j)
      enddo
c
c      evaluation of the fluxes
c      ------------------------
      fu=-vtrmins
      fv=-fu
c  
      endif
c
      return
      end
