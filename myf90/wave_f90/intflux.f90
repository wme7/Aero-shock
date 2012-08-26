      subroutine intflux(i,fu,fv,uu,vv,j)
!     ***********************************
!
      use mod_wave
      implicit none
    
!
      integer, intent(in)::i,j
      double precision,dimension(nx,k+1,0:kt),intent(in)::uu,vv
      double precision,intent(out)::fu,fv
      double precision::utrplus,vtrplus,utrmins,vtrmins
      double precision::uave,vave,ujmp,vjmp,uhat,vhat
      integer::l
!
!      evaluation of the traces
!      ------------------------
      utrplus=0.d00
      utrmins=0.d00
!
      vtrplus=0.d00
      vtrmins=0.d00
!
      do l=0,k
       utrplus=utrplus+(-1)**l*uu(i+1,l+1,j)
       utrmins=utrmins+uu(i,l+1,j)
!
       vtrplus=vtrplus+(-1)**l*vv(i+1,l+1,j)
       vtrmins=vtrmins+vv(i,l+1,j)
      enddo
!
!      evaluation of the fluxes
!      ------------------------
      uave=(utrplus+utrmins)/2.0d00
      vave=(vtrplus+vtrmins)/2.0d00
!
      ujmp=utrplus-utrmins
      vjmp=vtrplus-vtrmins
!
      vhat=vave+c11(i)*ujmp-c12(i)*vjmp
      uhat=uave+c12(i)*ujmp+c22(i)*vjmp
!      vhat=vave
!      uhat=uave
!
      fu=-vhat
      fv=-uhat
!
      return
      end subroutine intflux
