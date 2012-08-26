      subroutine cdt
c     **************
c
      include 'com'
c
      dimension ccfl(0:KK,KK+1)
c
      data ccfl(0,1)/1.d0/
c
      data ccfl(1,2)/.333d0/
      data ccfl(1,3)/.409d0/
c
      data ccfl(2,3)/.209d0/
      data ccfl(2,4)/.235d0/
      data ccfl(2,5)/.271d0/
c
      data ccfl(3,3)/.130d0/
      data ccfl(3,4)/.145d0/
      data ccfl(3,5)/.167d0/
      data ccfl(3,6)/.185d0/
      data ccfl(3,7)/.206d0/
c
      data ccfl(4,3)/.089d0/
      data ccfl(4,4)/.100d0/
      data ccfl(4,5)/.115d0/
      data ccfl(4,6)/.127d0/
      data ccfl(4,7)/.142d0/
      data ccfl(4,8)/.154d0/
      data ccfl(4,9)/.168d0/
c
      data ccfl(5,3)/.066d0/
      data ccfl(5,4)/.073d0/
      data ccfl(5,5)/.085d0/
      data ccfl(5,6)/.093d0/
      data ccfl(5,7)/.104d0/
      data ccfl(5,8)/.114d0/
      data ccfl(5,9)/.124d0/
      data ccfl(5,10)/.134d0/
      data ccfl(5,11)/.144d0/
c
c      --- cmax = maximum speed
      if (nproblem.eq.1.or.nproblem.eq.2) then
         cmax=dabs(cc)
         if (cmax.lt.1.d00) then
            cmax=1.d00
         endif
      endif
c
c      --- cdt = dt
c
      if (k.gt.6.or.k.lt.0.or.
     *    kt.gt.12.or.kt.lt.1.or.
     *    kt.gt.2*k+1) then
       write(06,111)
 111   format(t5,'This scheme might be unstable!')
      else
c       --- ensuring the stability of the scheme
       if (cfl.gt.1.d00*ccfl(k,kt)) cfl=1.d00*ccfl(k,kt)
      endif
c
      dt=cfl*dx/cmax
c
      return
      end
