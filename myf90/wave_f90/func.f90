      double precision function func(y)
!     ---------------------------------
!
      use mod_wave
      implicit none
      double precision,intent(in):: y
      double precision,external::phifun1
!
      integer:: p
      double precision:: h,ww,xxx
      parameter (p=10,h=1.d00,ww=5.d-01)
!
!      --- bump
!      if (nproblem.eq.1.or.nproblem.eq.2) then 
!     
       xxx=(y+a+ww)/ww
       func=phifun1(xxx)
!       if (dabs(xxx).lt..5d00) then
!        func=h*(2.0d00*xxx-1)**p*(2.0d00*xxx+1)**p
!       else
!        func=0.d00
!       endif
!
!      endif
!
      return
      end function func

