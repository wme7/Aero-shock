c**************************
c*     --- COMMONS        *
c**************************
      include 'com'
c
c*******************************
c*     --- INSTRUCTIONS        *
c*******************************
c
c      --- Self-introduction
c      --- *****************
      write(06,999)
  999 format(//t5,' This code solves numerically the'
     *        /t5,' initial boundary value problem for the'
     *        /t5,' transport equation:'
     *       //t5,' du/dt + du/dx = 0 in (0,T)x(0,1)'
     *        /t5,' u(x,t=0)=phi(x)   for x in (0,1)'
     *        /t5,' with periodic boundary conditions.'///)
c
c      --- Entering Mode of Operation:
c      --- ***************************
         iread=05
c
c      --- PREPROCESSING.
c      --- **************
      ne=0
   99 ne=ne+1
c
c         --- Checking the current number of different meshes
         if (ne.gt.100) then
            write(06,98)
   98       format(t5,' Only 100 different meshes are allowed!')
            stop
         endif
c
c         --- data input.
         call input
c
c         --- initialization
         call init
c
c      --- PROCESSING.
c      --- ***********
c
c         --- computation of the approximate solution uh.
         open(unit=13,file='apprx')
         open(unit=12,file='exact')
c
         call calc
c
c         --- computation of the error [Ph(u)-uh].
c         --- Ph(u) is the interpolant of u used to
c         --- define the initial condition.
c         --- The error is computed only if the exact solution
c         --- has been programmed, i.e., if nexist = 1.
         if (nexist.eq.1) call errors
c
c      --- FLOW-CONTROL.
c      --- *************
c
c      --- ind=1 => we are going to refine the mesh,
c      --- all the remaining data is left unchanged;
c      --- ind=0 => no refinement => go to POST-PROCESSING
      write(06,02)
 02   format(t5,'another discretization ?: yes(1), no(0)')
      read (05,*) ind  
      if (ind.ne.0) go to 99
c
c      --- POST-PROCESSING.
c      --- ****************
c
c         --- computing the rates of convergence and constants,
c         --- and opening the results file "res" and printing 
c         --- the results in it.
         if (nexist.eq.1) then
            call ordr
c
            open(unit=20,file='res')
            call impres
c
         endif     
c
c****************
c*    END       *
c****************
      end
