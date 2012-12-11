      subroutine input
c     ****************
c
      include 'com'
c
c      --- input
      if (ne.eq.1) then
         write(06,210)
  210    format(t5,'problem number :'/
     *          t20,' 1 :'
     *         /t28,' phi(x) = sin(2 pi x)'/
     *          t20,' 2 :'
     *         /t28,' phi(x) = 1  for x in (.25, .75)'/
     *         /t34,      '  = 0  otherwise '/
     *          t20,' 3 :'
     *         /t28,' phi(x) = 64*(x-.25)**2*(x-.75)**2',
     *              ' for x in (.25, .75)'/
     *         /t34,      '  = 0  otherwise ')
c
         write(06,31)
   31    format(t5,'problem number :')
         read(05,*) nproblem 
c
c         --- computing nexist
         if (nproblem.eq.1) nexist=1
         if (nproblem.eq.2) nexist=1
         if (nproblem.eq.3) nexist=1
c
         rl=1.d00
c
         write(06,310)
310      format(t5,'Polynomial degree k:')
         read(05,*) k
c
         write(06,311)
 311     format(t5,'Order of the RK method:')
         read(05,*) kt
c
c         write(06,22)
c   22    format(t5,'scheme number :'/
c    *          t23,'  1: Upwinding')
c
c         write(06,32)
c   32    format(t5,'scheme number :')
c         read(05,*) nsch                       

         nsch =1 
c
         write(06,39)
   39    format(t5,'Final time T :')
         read(05,*) t
c
         cfl=.99d00/(2*k+1)
         write(06,400)
 400     format(t5,'Do we take cfl=.99/(2k+1) ? yes(1), no(0)')
         read(05,*) iread
c
         if (iread.ne.1) then
            write(06,40)
   40       format(t5,'the cfl number:')
            read(05,*) cfl
         endif
      endif 
c
      call cdt          
c
      write(06,41)
   41 format(t5,'number of intervals in [0,1] :')
      read(05,*) nx
c
      dx=rl/nx
      ddx(ne)=dx
c
c      computing dt
      dt=cfl*dx
      nt=t/dt
c
      return
c
      end
