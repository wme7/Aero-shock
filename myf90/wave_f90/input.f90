      subroutine input(inter)
!     ****************
!
      use mod_wave
      implicit none
      integer, intent(in)::inter
      integer::idetail,j
      double precision::nnxx
!
!      --- input
      if (inter.le.0) then
!
!         --- reading everything from the file 'data'
       else
!
!         --- reading everything from the screen
         if (ne.eq.1) then
!            write(06,30)
   30       format(t5,' detailed questions? ( 1 -> yes, 0 -> no )')
!            read(05,*) idetail
            idetail=1
!            
            if (idetail.eq.1) then
               write(06,210)
  210          format(t5,'problem number :'/&
                     t20,' 1 :'&
                    /t28,' sound speed c=1 on (-a,a)'&
                    /t28,' transparent boundaries at a=2'&
                    /t28,' phi(t) = bump')
!
               write(06,211)
 211           format(t5,'problem number :'/&
                     t20,' 2 :'&
                    /t28,' sound speed c=2 on (-L,L), L=1'&
                    /t28,' transparent boundaries at a=2'&
                    /t28,' phi(t) = bump')
               write(06,*) "problem 3 periodic boundary condition"
            endif
!
            write(06,31)
   31       format(t5,'problem number :')
            read(05,*) nproblem 
            nexist=1
!
            if (nproblem.eq.1) then
             cc=1.d00
             rl=1.d00
             a =2.d00
            endif

            if (nproblem.eq.3) then
             cc=1.d00
             rl=1.d00
             a =2.d00
            endif

!
            if (nproblem.eq.2) then
             cc=2.d00
             rl=1.d00
             a =2.d00
            endif           
!
            write(6,*) "Post-Processing? 0:No, 1:Yes"
            read(05,*) npost
            write(06,310)
 310         format(t5,'Polynomial degree k:')
            read(05,*) k

            ks=2*k+1
            if (npost .eq. 1) then
             kt=ks
            else
            write(06,311)
 311        format(t5,'Order of the RK method:')
            read(05,*) kt
            endif
!            nbc=2
            write(06,*) "boudary treatment 0:No, 1:test1"
            read(05,*) nbc
!
!            --- computing nexist
!            --- ***************************
            if (nproblem.eq.1) nexist=1
            if (nproblem.eq.2) nexist=1
            if (nproblem.eq.3) nexist=1

!
            if (idetail.eq.1) then
               write(06,22)
   22          format(t5,'scheme number :'/&
                     t23,'  1: Upwinding')
            endif
!
            write(06,32)
   32       format(t5,'scheme number :')
             nsch=1
!            read(05,*) nsch                         
!
!            write(06,38)
!   38       format(t30,'1->save exact and approx solutions at the end ')
!            read(05,*) imp
!
            write(06,*) "Initial Time"
!            read(05,*) t0
             t0=1.d0
             t0=0.d0

            write(06,39)
   39       format(t5,'Final time T :')
            read(05,*) t
            do j=1,80
            tf(j)=(t*j)/80d0
            enddo
!
            write(06,40)
   40       format(t5,'cfl=|c|dt/dx :')
            read(05,*) cfl
         write(06,*) "level of refinement=?"
         read(05,*) nlevel
            allocate(errl2(nlevel+1))
            allocate(errl1(nlevel+1))
            allocate(errli(nlevel+1))
            allocate(ordl2(nlevel+1)) 
            allocate(ordl1(nlevel+1))
            allocate(ordli(nlevel+1))
            allocate(pordl2(nlevel+1))
            allocate(pordl1(nlevel+1))
            allocate(pordli(nlevel+1))
            allocate(constl2(nlevel+1))
            allocate(constl1(nlevel+1))
            allocate(constli(nlevel+1))
            allocate(serrli(nlevel+1))
            allocate(perrli(nlevel+1))
            allocate(perrl1(nlevel+1))
            allocate(perrl2(nlevel+1)) 
            allocate(ddx(nlevel+1))
            endif
!
         write(06,41)
   41    format(t5,'number of intervals in [-L,L] :')

!
         nnxx=2*2**(ne+1)
         dx=2.d00*rl/nnxx
         nx=int(a/rl)*nnxx
         write(6,*) 'nx=',nx
!
!         computing dt
         call cdt
          write(6,221) dx,dt
 221    format('dx=',e12.4,';','dt=', e12.4,';')  

!(cfl,k,kt)
         nt=int((t-t0)/dt)
!         write(6,*) nt

         return
!
         endif
      end subroutine input
