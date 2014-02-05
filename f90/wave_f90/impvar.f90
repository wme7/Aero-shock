      subroutine impvar
!     *****************
!
      use mod_wave
!
!      printing out u 
!      --------------     
      write(12,1) k,a,nx
 1    format( &
      t7,'clear all',/ &
      t7,'k=', I5,';',/ &
      t7,'a=', f7.3,';',/ &
      t7,'nx=', I5,';',/ &
      t7,'dx=2*a/nx;',/ &
      t7,' pts=[-a:dx:a];',/)
!
      do i=1,nx
         write(12,*) 'u(',i,',:)=[',(u(i,j,0),j=1,k+1),'];'
      end do
      write(12,2)
 2    format(&
       t7,"n_int_pts=5*k+1;"/ &
       t7,"dx_int=2/n_int_pts;"/ &
       t7,"x=[-1:dx_int:1]"/ &             
       t7,"leg=zeros(k+1,n_int_pts+1);"/ &
       t7,"for m=0:k;"/ &
       t7,"    aux=legendre(m,x);"/ &
       t7,"    leg(m+1,:)=aux(1,:);"/ &
       t7," end;"/ &
       t7," figure(1)"/ &
       t7," clf"/ &
       t7,"set(axes,'fontsize', 20)"/ &
       t7," for i=1:nx;"/ &
       t7,"    dx_int_local=(pts(i+1)-pts(i))/n_int_pts;"/ &
       t7,"    x_local=[pts(i):dx_int_local:pts(i+1)];"/ &
       t7,"    f=u(i,:)*leg;"/ &
       t7,"    plot(x_local,f);"/ &
       t7,"    hold on"/ &
       t7," end;")
!
!      printing out v 
!      --------------     
      do i=1,nx
         write(12,*) 'v(',i,',:)=[',(v(i,j,0),j=1,k+1),'];'
      end do
      write(12,3)
 3    format(&
       t7,"n_int_pts=5*k+1;"/ &
       t7,"dx_int=2/n_int_pts;"/ &
       t7,"x=[-1:dx_int:1]"/ &             
       t7,"leg=zeros(k+1,n_int_pts+1);"/ &
       t7,"for m=0:k;"/ &
       t7,"    aux=legendre(m,x);"/ &
       t7,"    leg(m+1,:)=aux(1,:);"/ &
       t7," end;"/ &
       t7," figure(2)"/ &
       t7," clf"/ &
       t7,"set(axes,'fontsize', 20)"/ &
       t7," for i=1:nx;"/ &
       t7,"    dx_int_local=(pts(i+1)-pts(i))/n_int_pts;"/ &
       t7,"    x_local=[pts(i):dx_int_local:pts(i+1)];"/ &
       t7,"    f=v(i,:)*leg;"/ &
       t7,"    plot(x_local,f);"/ &
       t7,"    hold on"/ &
       t7," end;")
!
      return
      end subroutine impvar
