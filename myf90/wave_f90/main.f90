      program wave
!**************************
!*     --- COMMONS        *
!**************************
        use mod_wave
        implicit none
        integer::iread,inter,i
!*******************************
!*     --- INSTRUCTIONS        *
!*******************************
!
!      --- Self-introduction
!      --- *****************
      write(06,999)
  999 format(//t5,' This code solves numerically the'&
             /t5,' initial boundary value problem for the'&
             /t5,' Wave equation:'&
            //t5,' 1/c^2 du/dt - dv/dx = 0, in (0,T)x(-a,a)'&
            //t5,'       dv/dt - du/dx = 0, in (0,T)x(-a,a)'&
             /t5,' u(x,t=0)=phi(x), v(x,t=0)=-phi(x),  on (-a,a)'&
             /t5,' u+cv=0 at x=a (c=1 at x=a,-a),'&
             /t5,' u-cv=2 phi(x-t) at x=-a for t in (0,T)'///)
!      --- Entering Mode of Operation:
!      --- ***************************
         iread=05
!         write(06,01)
   01    format(t5,' enter the operation mode: ( 1 -> interactive,'/&
               t33,' 0 -> read data in file DATA)')
!         read(05,*) inter
            inter=1
         if (inter .eq.0) then
            iread=10
            open(unit=10,file='data')
         endif
!      --- PREPROCESSING.
!      --- **************
            open(unit=20,file='res')
         open(unit=12,file='sol.m')
         open(unit=16,file='err.m')
!          write(12,*) "figure(1)"
!          write(12,*) "set(axes,'fontsize', 20)"
!           write(12,*) "figure(2)"
!          write(12,*) "set(axes,'fontsize', 20)"
      ne=0
   99 ne=ne+1
      if (ne .gt.1) then
      write(6,*) ne
      deallocate(alpha)
      deallocate(x)
      deallocate(xt)
      deallocate(w)
      deallocate(wt)
      deallocate(pl)
      deallocate(plt)
!      deallocate(cbin)
!      deallocate(dl)
      deallocate(u)
      deallocate(v)   
      deallocate(fluxu)
      deallocate(fluxv)
      deallocate(f)
      deallocate(AA)
      deallocate(c11)
      deallocate(c12)  
      deallocate(c22)  
      deallocate(c) 
      endif 
!         --- Checking the current number of different meshes
         if (ne.gt.100) then
            write(06,98)
   98       format(t5,' Only 100 different meshes are allowed!')
            stop
         endif
!         --- data input.
         call input(inter)
!         --- initialization
         call claim
         call init
!      --- PROCESSING.
!      --- ***********
!         --- computation of the approximate solution uh.

         call calc
!         --- computation of the error [Ph(u)-uh].
!         --- Ph(u) is the interpolant of u used to
!         --- define the initial condition.
!         --- The error is computed only if the exact solution
!         --- has been programmed, i.e., if nexist = 1.

      

         if (nexist.eq.1) call errors
         if (npost.eq.1) then
         open(unit=14,file='sop.m')
         open(unit=15,file='uex.m')
         
         call superconv
!         call post
         call eno
         endif

!      --- FLOW-CONTROL
!      --- *************
!      --- ind=1 => we are going to refine the mesh,
!      --- all the remaining data is left unchanged;
!      --- ind=0 => no refinement => go to POST-PROCESSING
!      if (inter.gt.0) write(06,02)
!   02 format(t5,'another discretization ?')
!      read (iread,*) ind  
!      if (ind.ne.0) go to 99
           ddx(ne)=dx
!           write(20,*) dx,errli(ne),errl1(ne)
       if (ne .lt. nlevel) goto 99
!      --- POST-PROCESSING.
!      --- ****************
!         --- opening the approximate solution file "sol",
!         --- and saving uh in it.
!         if (imp.ne.0) then
!            open(unit=12,file='sol.m')
            call impvar
!            call errout
        write(12,51) 
        write(12,61) t
        write(12,31)
        write(12,41) nt
        write(12,52)
        write(12,62) t 
        write(12,32)    
        write(12,42) nt

        write(16,51)   
        write(16,61) t 
        write(16,31)   
        write(16,41) nt
        write(16,*) 'figure(3)'   
        write(16,63) t 
        write(16,33)   
        write(16,43) nt

   51  format('figure(1)')
   52  format('figure(2)')
   61  format("title('U at T=",f7.3,"');")
   62  format("title('V at T=",f7.3,"');")
   63  format("title('Error at T=",f7.3,"');")
   64  format("title('U_a at T=",f7.3,"');")

   31  format('axis([-2 2 -0.5 1.5]);')
   32  format('axis([-2 2 -1.5 0.5]);')
   33  format('axis([-2 2 -0.03 0.03]);')

   41  format("saveas(gcf, 'u",I4.4,".eps');")
   42  format("saveas(gcf, 'v",I4.4,".eps');")
   43  format("saveas(gcf, 'er",I4.4,".eps');")

!!         endif         

!         --- computing the rates of convergence and constants,
!         --- and opening the results file "res" and printing 
!         --- the results in it.
         if (nexist.eq.1) then
            call ordr
            call pordr
         endif
!           write(6,*) "      Li                       L1"  

           do i=1,ne
           write(6,21) ordli(i),ordl1(i),ordl2(i)

           enddo
           write(6,*) "post-processing"
!           write(6,*) "      Li                       L1"

           do i=1,ne
           write(6,21) pordli(i),pordl1(i),pordl2(i)
           enddo

! Printout the table of the result
           write(20,91) 
           write(20,112)
       do i=1,ne
           write(20,7) ddx(i),errli(i),ordli(i), &
            errl2(i),ordl2(i),errl1(i),ordl1(i) 
       enddo
          write(20,92) k,kt
           write(20,91) 
           write(20,112)
       do i=1,ne
           write(20,7) ddx(i),perrli(i),pordli(i), &
            perrl2(i),pordl2(i),perrl1(i),pordl1(i) 
       enddo
          write(20,93) k,kt

 7    format(f10.4,'&',e8.2,'&',f6.2,'&',e8.2,'&',f6.2,/ &
      '&',e8.2,'&',f6.2,'\\ \hline')
  91 format('\begin{table}[htb] \begin{center}', &
        '\begin{tabular}{|c|c|c|c|c|c|c|}\hline')
  92 format('\end{tabular}\caption{Rate of convergence for $P^', &
           i2.1,'$-elements, SSP-RK',I2.1,'.}\end{center}\end{table}')
  93 format('\end{tabular}\caption{Rate of convergence for $P^', &
           i2.1,'$-elements, SSP-RK',I2.1,'After Post-Processing', &
           '.}\end{center}\end{table}')
 112 format('$h$&$L^{\infty}$-error&order&$L2$-error&order&$L1$', &
       'error&order\\')

 21    format('Li=',f8.4,';','L1=', f8.4,';','L2=',f8.4,';')  

!   05 continue       
!
!            call impres
!!
!         endif     
!         write(6,*) nexist
!          write(6,*) "L1-Error",errl1
!          write(6,*) "Li-Error",errli
      deallocate(alpha)
      deallocate(x)
      deallocate(xt)
      deallocate(w) 
      deallocate(wt)
      deallocate(pl)
      deallocate(plt)
      deallocate(cbin)
      deallocate(dl)  
      deallocate(u)   
      deallocate(v)   
      deallocate(fluxu)
      deallocate(fluxv)
      deallocate(f)
      deallocate(AA)
      deallocate(c11)  
      deallocate(c12)  
      deallocate(c22)  
      deallocate(c)    
      deallocate(errl1)
      deallocate(errli)
            deallocate(ordl1)
            deallocate(ordli)
            deallocate(pordl1)
            deallocate(pordli)
            deallocate(constl1)
            deallocate(constli)
            deallocate(serrli)


!****************
!*    END       *
!****************
      end program wave
