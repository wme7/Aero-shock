      double precision function phifun1(xxx)
!     ---------------------------------
!
      implicit none
      double precision,intent(in):: xxx
!
      integer:: p
      double precision:: h
      parameter (p=10,h=1.d00)
!
!      --- bump
!      if (nproblem.eq.1.or.nproblem.eq.2) then 
!     
       if (dabs(xxx).lt..5d00) then
        phifun1=h*(2.0d00*xxx-1)**p*(2.0d00*xxx+1)**p
       else
        phifun1=0.d00
       endif
!
!
      return
      end function phifun1

      double precision function phifun2(xxx)
!     ---------------------------------
!
      implicit none
      double precision,intent(in):: xxx
!
      integer:: p
      double precision:: h
      parameter (p=10,h=1.d00)
!
!      --- bump
!      if (nproblem.eq.1.or.nproblem.eq.2) then 
!     
       if (dabs(xxx).lt..5d00) then
        phifun2=h*(0.5d00-dabs(xxx))*2.d0
       else
        phifun2=0.d00
       endif
!
!
      return
      end function phifun2
