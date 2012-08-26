      subroutine exac(ua,va,y)
!     ---------------------------------
!
      use mod_wave
      implicit none
      double precision,intent(in):: y
      double precision,intent(out)::ua,va
      double precision,external::func
      double precision::z,uo,t2
!
      integer:: p,nnt
      double precision:: h,ww
      parameter (p=10,h=1.d00,ww=5.d-01)
!
!      --- bump
      z=y-t 
      if (nproblem.eq.1) then
        ua=func(z)
        va=-ua
      endif
      if (nproblem.eq.3) then
       nnt=int(t/4)
        z=y-(t-4.d0*nnt)
        ua=func(z)
        va=-ua
      endif

      if (nproblem.eq.2) then
        if (t .le. 1.25d0) then  
         ua=func(z)
         va=-ua
        elseif (t .le. 1.75d0) then
         if (y .gt. -rl) then
         t2=(y+rl)/2.d0
         z=-rl-(t-t2)
!         z=y-2.d0*(t-1.25d0)-1.25d0
         va=-2.d0*func(z)/3.d0
         ua=-2.d0*va
         else
         t2=-(y+rl)
          z=-rl-(t-t2)
!         z=y+t-1.25d0-1.25d0
         ua=func(z)/3.d0
         va=ua
         z=y-t
         uo=func(z)
         ua=ua+uo
         va=va-uo
         endif
        else

         if (y .gt. -rl) then
         t2=(y+rl)/2.d0
         z=-rl-(t-t2)  
         va=-2.d0*func(z)/3.d0
         ua=-2.d0*va
         else
         t2=-(y+rl)
          z=-rl-(t-t2)
         ua=func(z)/3.d0
         va=ua
         endif         

        endif
      endif
!
      return
      end subroutine exac

