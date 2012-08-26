      subroutine cdt
!     **************
!
      use mod_wave
!
      implicit none 
      double precision,dimension(0:8,12)::ccfl
      double precision::cmax      
!
      data ccfl(0,1)/1.d0/
      data ccfl(0,2)/1.d0/
      data ccfl(0,3)/1.256d0/
      data ccfl(0,4)/1.396d0/
      data ccfl(0,5)/1.608d0/
      data ccfl(0,6)/1.776d0/
      data ccfl(0,7)/1.997d0/
      data ccfl(0,8)/2.156d0/
      data ccfl(0,9)/2.350d0/
      data ccfl(0,10)/2.534d0/
      data ccfl(0,11)/2.725d0/
      data ccfl(0,12)/2.911d0/

!
      data ccfl(1,2)/.333d0/
      data ccfl(1,3)/.409d0/
      data ccfl(1,4)/0.464d0/
      data ccfl(1,5)/0.534d0/
      data ccfl(1,6)/0.592d0/
      data ccfl(1,7)/0.659d0/
      data ccfl(1,8)/0.718d0/
      data ccfl(1,9)/0.783d0/
      data ccfl(1,10)/0.844d0/
      data ccfl(1,11)/0.908d0/
      data ccfl(1,12)/0.970d0/

!
      data ccfl(2,3)/.209d0/
      data ccfl(2,4)/.235d0/
      data ccfl(2,5)/.271d0/
      data ccfl(2,6)/.300d0/
      data ccfl(2,7)/.333d0/
      data ccfl(2,8)/.364d0/
      data ccfl(2,9)/.396d0/
      data ccfl(2,10)/.428d0/
      data ccfl(2,11)/.460d0/
      data ccfl(2,12)/.491d0/

!
      data ccfl(3,3)/.130d0/
      data ccfl(3,4)/.145d0/
      data ccfl(3,5)/.167d0/
      data ccfl(3,6)/.185d0/
      data ccfl(3,7)/.206d0/
      data ccfl(3,8)/.225d0/
      data ccfl(3,9)/.245d0/
      data ccfl(3,10)/.264d0/
      data ccfl(3,11)/.284d0/
      data ccfl(3,12)/.303d0/

!
      data ccfl(4,3)/.089d0/
      data ccfl(4,4)/.100d0/
      data ccfl(4,5)/.115d0/
      data ccfl(4,6)/.127d0/
      data ccfl(4,7)/.142d0/
      data ccfl(4,8)/.154d0/
      data ccfl(4,9)/.168d0/
      data ccfl(4,10)/.182d0/
      data ccfl(4,11)/.195d0/
      data ccfl(4,12)/.209d0/
!
      data ccfl(5,3)/.066d0/
      data ccfl(5,4)/.073d0/
      data ccfl(5,5)/.085d0/
      data ccfl(5,6)/.093d0/
      data ccfl(5,7)/.104d0/
      data ccfl(5,8)/.114d0/
      data ccfl(5,9)/.124d0/
      data ccfl(5,10)/.134d0/
      data ccfl(5,11)/.144d0/
      data ccfl(5,12)/.153d0/
!
      data ccfl(6,3)/.051d0/
      data ccfl(6,4)/.056d0/
      data ccfl(6,5)/.065d0/
      data ccfl(6,6)/.072d0/
      data ccfl(6,7)/.080d0/
      data ccfl(6,8)/.087d0/
      data ccfl(6,9)/.095d0/
      data ccfl(6,10)/.103d0/
      data ccfl(6,11)/.111d0/
      data ccfl(6,12)/.118d0/
!
      data ccfl(7,3)/.040d0/
      data ccfl(7,4)/.045d0/
      data ccfl(7,5)/.052d0/
      data ccfl(7,6)/.057d0/
      data ccfl(7,7)/.064d0/
      data ccfl(7,8)/.070d0/
      data ccfl(7,9)/.076d0/
      data ccfl(7,10)/.082d0/
      data ccfl(7,11)/.088d0/
      data ccfl(7,12)/.094d0/
!
      data ccfl(8,3)/.033d0/ 
      data ccfl(8,4)/.037d0/
      data ccfl(8,5)/.042d0/
      data ccfl(8,6)/.047d0/
      data ccfl(8,7)/.052d0/
      data ccfl(8,8)/.057d0/
      data ccfl(8,9)/.062d0/
      data ccfl(8,10)/.067d0/
      data ccfl(8,11)/.072d0/
      data ccfl(8,12)/.077d0/
!
!      --- cmax = maximum speed
!      if (nproblem.eq.1.or.nproblem.eq.2) then
         cmax=dabs(cc)
         if (cmax.lt.1.d00) then
            cmax=1.d00
         endif
!      endif

!      --- cdt = dt

      if (k.gt.8.or.k.lt.0.or.  &
         kt.gt.12.or.kt.lt.1.or.  &
        kt.gt.2*k+1) then
       write(06,111)
 111   format(t5,'This scheme might be unstable!')
      else
!       --- ensuring the stability of the scheme
       if (cfl.gt.1.d00*ccfl(k,kt)) cfl=1.d00*ccfl(k,kt)
      endif
       if (cfl.gt.1.d00*ccfl(k,kt)) cfl=1.d00*ccfl(k,kt)
!
      dt=cfl*dx/cmax
      dt0=dt
!      write(6,*) "dt=",dt
!
      return
      end subroutine cdt
