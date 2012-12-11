      double precision function func(y,tt)
c     ------------------------------------
c
      include 'com'
      double precision y,tt
c
c      --- sinusoidal
      if (nproblem.eq.1) then 
c     
       twopi=2.d00*dacos(-1.d00)
       func=sin(twopi*(y-tt)) 
      return
c
      endif
c
c      --- step
      if (nproblem.eq.2) then 
c
       if (y-tt.ge.0.d00) then
         yy=y-tt
       else     
         yy=(y-tt)-int(y-tt)+1.d00
       endif
c     
       if (yy.lt..25d00.or.yy.gt..75d00) then
        func=0.d00
       else
        func=1.d00
       endif
       return
      endif
c
c      --- bump
      if (nproblem.eq.3) then 
c
       if (y-tt.ge.0.d00) then
         yy=y-tt
       else     
         yy=(y-tt)-int(y-tt)+1.d00
       endif
c
       if (abs(yy-.5).lt..25d00) then
        func=64*(yy-.25d00)**2*(yy-.75d00)**2
       else
        func=0.d00
       endif
       return
c
      endif
c
      return
      end
