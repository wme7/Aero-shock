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
     *        /t5,' Wave equation:'
     *       //t5,' 1/c^2 du/dt - dv/dx = 0, in (0,T)x(-a,a)'
     *       //t5,'       dv/dt - du/dx = 0, in (0,T)x(-a,a)'
     *        /t5,' u(x,t=0)=phi(x), v(x,t=0)=-phi(x),  on (-a,a)'
     *        /t5,' u+cv=0 at x=a (c=1 at x=a,-a),'
     *        /t5,' u-cv=2 phi(x-t) at x=-a for t in (0,T)'///)
c
c      --- Entering Mode of Operation:
c      --- ***************************
         iread=05
         write(06,01)
   01    format(t5,' enter the operation mode: ( 1 -> interactive,'/
     *          t33,' 0 -> read data in file DATA)')
         read(05,*) inter
         if (inter .eq.0) then
            iread=10
            open(unit=10,file='data')
         endif
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
         call input(inter)
c
c         --- initialization
         call init
c
c      --- PROCESSING.
c      --- ***********
c
c         --- computation of the approximate solution uh.
         open(unit=12,file='sol.m')
         call calc
c
c         --- computation of the error [Ph(u)-uh].
c         --- Ph(u) is the interpolant of u used to
c         --- define the initial condition.
c         --- The error is computed only if the exact solution
c         --- has been programmed, i.e., if nexist = 1.
ccccccc         if (nexist.eq.1) call errors
c
c      --- FLOW-CONTROL.
c      --- *************
c
c      --- ind=1 => we are going to refine the mesh,
c      --- all the remaining data is left unchanged;
c      --- ind=0 => no refinement => go to POST-PROCESSING
      if (inter.gt.0) write(06,02)
   02 format(t5,'another discretization ?')
      read (iread,*) ind  
      if (ind.ne.0) go to 99
c
c      --- POST-PROCESSING.
c      --- ****************
c
c         --- opening the approximate solution file "sol",
c         --- and saving uh in it.
ccccc         if (imp.ne.0) then
ccccc            open(unit=12,file='sol.m')
            call impvar
ccccc         endif         
c
c         --- computing the rates of convergence and constants,
c         --- and opening the results file "res" and printing 
c         --- the results in it.
ccccc         if (nexist.eq.1) then
ccccc            call ordr
cccccc
ccccc            open(unit=20,file='res')
ccccc            call impres
cccccc
ccccc         endif     
c
c****************
c*    END       *
c****************
      end
