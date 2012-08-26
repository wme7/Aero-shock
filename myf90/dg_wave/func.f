      double precision function func(y)
c     ---------------------------------
c
      include 'com'
      double precision y
c
      integer p
      double precision h,ww
      parameter (p=10,h=1.d00,ww=5.d-01)
c
c      --- bump
      if (nproblem.eq.1.or.nproblem.eq.2) then 
c     
       xxx=(y+a+ww)/ww
       if (dabs(xxx).lt..5d00) then
        func=h*(2.0d00*xxx-1)**p*(2.0d00*xxx+1)**p
       else
        func=0.d00
       endif
c
      endif
c
      return
      end
