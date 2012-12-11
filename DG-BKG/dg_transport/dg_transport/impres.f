      subroutine impres
c     *****************
c      
      include 'com'
c
c      --- title
      write(20,110) t,rl
c
      if (nproblem.eq.1) write(20,201)
      if (nproblem.eq.2) write(20,202)
      if (nproblem.eq.3) write(20,203)
c
c      --- method
      write(20,300) k,kt
c
c      --- parameters
      write(20,611) (1./ddx(i),cfl,i=1,ne)
c
c      --- L2-errors
      write(20,640) (errl1(i),constl1(i),ordl1(i),i=1,ne)
c
c      --- L-infinity-errors
      write(20,740) (errli(i),constli(i),ordli(i),i=1,ne)
c
      return
c
c      --- formats 
 110  format(//t5,' du/dt + du/dx = 0 in (0,T=',f5.2,')x(0,L=',f4.2,
     *     ')'/t5,' u(t=0)=uo         on (0,L)'
     *        /t5,' periodic boundary conditions'/)
c
 201  format(t5,' uo(x) = sin(2 pi x)')
 202  format(t5,' uo(x) = 1  for x in (.25, .75)'/
     *      /t11,      '  = 0  otherwise ')
 203  format(t5,' uo(x) = 64*(x-.25)**2*(x-.75)**2',
     *          ' for x in (.25, .75)'/
     *      /t11,      ' = 0  otherwise ')
c
 300  format(/t6,'Polynomial degree:', i5,
     *       /t6,'Order of RK method:', i4)
c
  611 format( /t5,'     1/dx      cfl'
     *        /t5,'  -------------------'/
     *       (t5,f9.2,f11.6))
c
  640 format(//t9,'L2-error',t26,'constant',t42,'order'/t6,44('-')
     *        /(t5,3e15.8))
c 
  740 format(//t9,'Li-error',t26,'constant',t42,'order'/t6,44('-')
     *        /(t5,3e15.8))
c
      end
