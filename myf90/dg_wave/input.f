      subroutine input(inter)
c     ****************
c
      include 'com'
c
c      --- input
      if (inter.le.0) then
c
c         --- reading everything from the file 'data'
       else
c
c         --- reading everything from the screen
         if (ne.eq.1) then
            write(06,30)
   30       format(t5,' detailed questions? ( 1 -> yes, 0 -> no )')
            read(05,*) idetail
c            
            if (idetail.eq.1) then
               write(06,210)
  210          format(t5,'problem number :'/
     *                t20,' 1 :'
     *               /t28,' sound speed c=1 on (-a,a)'
     *               /t28,' transparent boundaries at a=2'
     *               /t28,' phi(t) = bump')
c
               write(06,211)
 211           format(t5,'problem number :'/
     *                t20,' 2 :'
     *               /t28,' sound speed c=2 on (-L,L), L=1'
     *               /t28,' transparent boundaries at a=2'
     *               /t28,' phi(t) = bump')
            endif
c
            write(06,31)
   31       format(t5,'problem number :')
            read(05,*) nproblem 
c
            if (nproblem.eq.1) then
             cc=1.d00
             rl=1.d00
             a =2.d00
            endif
c
            if (nproblem.eq.2) then
             cc=2.d00
             rl=1.d00
             a =2.d00
            endif           
c
            write(06,310)
310         format(t5,'Polynomial degree k:')
            read(05,*) k
c
            write(06,311)
 311        format(t5,'Order of the RK method:')
            read(05,*) kt
c
c            --- computing nexist
c            --- ***************************
            if (nproblem.eq.1) nexist=1
            if (nproblem.eq.2) nexist=1
c
c            if (idetail.eq.1) then
c               write(06,22)
c   22          format(t5,'scheme number :'/
c     *                t23,'  1: Upwinding')
c            endif
c
c            write(06,32)
c   32       format(t5,'scheme number :')
c            read(05,*) nsch                         

             nsch=1
c
c            write(06,38)
c   38       format(t30,'1->save exact and approx solutions at the end ')
c            read(05,*) imp
c
            write(06,39)
   39       format(t5,'Final time T :')
            read(05,*) t
c
            write(06,40)
   40       format(t5,'cfl=|c|dt/dx :')
            read(05,*) cfl
            endif
c
         write(06,41)
   41    format(t5,'number of intervals in [-L,L] :')
         read(05,*) nnxx
c
         dx=2.d00*rl/nnxx
         nx=int(a/rl)*nnxx
c
c         computing dt
         call cdt(cfl,k,kt)
         nt=t/dt
c
         return
c
         endif
c
      end
