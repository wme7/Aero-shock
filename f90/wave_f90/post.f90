      subroutine post
!     ***************
!     post processing           
!     ***************
      use mod_wave
      implicit none
!      double precision,external::func
      double precision::fu,fv,ua,ue,aux,va,ve
      double precision::xic,aux1,xm,aux2,xl,aux3
      double precision,dimension(ks)::AAP
      integer::i,j,m,il,l,ip,n1
        n1=int((a-rl)/dx)

!      printing out u
!      --------------
      write(14,1) k,a,nx
      write(15,1) k,a,nx
      ip=0
 1    format( &
      t7,'clear all',/ &
      t7,'k=', I5,';',/ &
      t7,'a=', f7.3,';',/ &
      t7,'nx=', I5,';',/ &
      t7,'dx=2*a/nx;',/ &
      t7,' pts=[-a:dx:a];',/)
!      --- evaluating the fluxes
      do  i=1,nx-1
        call intflux(i,fu,fv,u,v,0)
        fluxu(i)=-fv
        fluxv(i)=-fu
      enddo

        call bdryflux(1,fu,fv,t,0)
        fluxu(0)=-fv
        fluxv(0)=-fu

        call bdryflux(2,fu,fv,t,0)
        fluxu(nx)=-fv
        fluxv(nx)=-fu

         perrl2(ne)=0.d0
         perrl1(ne)=0.d0
         perrli(ne)=0.d0
! 
        do i=1,ks
         aux2=0.d0
         aux3=0.d0
         il=-1
! Interpolate at nodes il+1,..., il+ks

          do j=1,ks
           aux=1.d0
           do l=1,j-1
           aux=aux*(j-l)*dx
           enddo
           do l=j+1,ks
           aux=aux*(j-l)*dx
           enddo
           AAP(j)=1.d0/aux
          enddo

         xic=-a+(2.d00*i-1.d00)*dx/2.d00

          do m=1,ks
          ua=0.d0
          va=0.d0
          xm=xic+x(m)*dx/2.d00
          call exac(ue,ve,xm)

          do j=1,ks
           aux=1.d0
           do l=1,j-1
          xl=-a+(il+l)*dx
           aux=aux*(xm-xl)
           enddo
           do l=j+1,ks
          xl=-a+(il+l)*dx
           aux=aux*(xm-xl)
           enddo

           ua=ua+AAP(j)*aux*fluxu(il+j)
           va=va+AAP(j)*aux*fluxv(il+j)
          enddo
         ip=ip+1
         write(14,*) 'x(',ip,')=[',xm,'];'
         write(14,*) 'ua(',ip,')=[',ua,'];'
         write(14,*) 'va(',ip,')=[',va,'];'
         write(15,*) 'x(',ip,')=[',xm,'];'
         write(15,*) 'ue(',ip,')=[',ue,'];'
         write(15,*) 've(',ip,')=[',ve,'];'


            aux1=dabs(ua-ue)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
            aux1=dabs(va-ve)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
          enddo
            perrl1(ne)=perrl1(ne)+aux2*dx
            perrl2(ne)=perrl2(ne)+aux3*dx
        enddo
! ***************************************
! One-Sided Interpolation
! ***************************************
! Before the discontinuity (x=-1)

        do i=ks+1,n1
         aux2=0.d0
         aux3=0.d0
         il=i-ks-1

          do j=1,ks
           aux=1.d0
           do l=1,j-1
           aux=aux*(j-l)*dx
           enddo
           do l=j+1,ks
           aux=aux*(j-l)*dx
           enddo
           AAP(j)=1.d0/aux
          enddo

         xic=-a+(2.d00*i-1.d00)*dx/2.d00

          do m=1,ks
          ua=0.d0
          ue=0.d0
          va=0.d0
          xm=xic+x(m)*dx/2.d00
           call exac(ue,ve,xm)
          do j=1,ks
           aux=1.d0
           do l=1,j-1
            xl=-a+(il+l)*dx
            aux=aux*(xm-xl)
           enddo
           do l=j+1,ks
            xl=-a+(il+l)*dx
            aux=aux*(xm-xl)
           enddo
           ua=ua+AAP(j)*aux*fluxu(il+j)
           va=va+AAP(j)*aux*fluxv(il+j)
          enddo
         ip=ip+1
         write(14,*) 'x(',ip,')=[',xm,'];'
         write(14,*) 'ua(',ip,')=[',ua,'];'
         write(14,*) 'va(',ip,')=[',va,'];'
         write(15,*) 'x(',ip,')=[',xm,'];'
         write(15,*) 'ue(',ip,')=[',ue,'];'
         write(15,*) 've(',ip,')=[',ve,'];'
            aux1=dabs(ua-ue)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
            aux1=dabs(va-ve)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1

          enddo
            perrl1(ne)=perrl1(ne)+aux2*dx
            perrl2(ne)=perrl2(ne)+aux3*dx
        enddo

        do i=n1+1,nx-ks
         aux2=0.d0
         aux3=0.d0
         il=i-1

          do j=1,ks
           aux=1.d0
           do l=1,j-1
           aux=aux*(j-l)*dx
           enddo
           do l=j+1,ks
           aux=aux*(j-l)*dx
           enddo
           AAP(j)=1.d0/aux
          enddo

         xic=-a+(2.d00*i-1.d00)*dx/2.d00

          do m=1,ks
          ua=0.d0
          ue=0.d0
          va=0.d0
          xm=xic+x(m)*dx/2.d00
           call exac(ue,ve,xm)
          do j=1,ks
           aux=1.d0
           do l=1,j-1
            xl=-a+(il+l)*dx
            aux=aux*(xm-xl)
           enddo
           do l=j+1,ks
            xl=-a+(il+l)*dx
            aux=aux*(xm-xl)
           enddo
           ua=ua+AAP(j)*aux*fluxu(il+j)
           va=va+AAP(j)*aux*fluxv(il+j)
          enddo
         ip=ip+1
         write(14,*) 'x(',ip,')=[',xm,'];'
         write(14,*) 'ua(',ip,')=[',ua,'];'
         write(14,*) 'va(',ip,')=[',va,'];'
         write(15,*) 'x(',ip,')=[',xm,'];'
         write(15,*) 'ue(',ip,')=[',ue,'];'
         write(15,*) 've(',ip,')=[',ve,'];'
            aux1=dabs(ua-ue)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
            aux1=dabs(va-ve)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1

          enddo
            perrl1(ne)=perrl1(ne)+aux2*dx
            perrl2(ne)=perrl2(ne)+aux3*dx
        enddo

        do i=nx-ks+1,nx
         aux2=0.d0
         aux3=0.d0
         il=nx-ks
          do j=1,ks
           aux=1.d0
           do l=1,j-1
           aux=aux*(j-l)*dx
           enddo
           do l=j+1,ks
           aux=aux*(j-l)*dx
           enddo
           AAP(j)=1.d0/aux
          enddo

         xic=-a+(2.d00*i-1.d00)*dx/2.d00

          do m=1,ks
          ua=0.d0
          ue=0.d0
          va=0.d0
          xm=xic+x(m)*dx/2.d00
          call exac(ue,ve,xm)
          do j=1,ks
           aux=1.d0
           do l=1,j-1
          xl=-a+(il+l)*dx
           aux=aux*(xm-xl)
           enddo
           do l=j+1,ks
          xl=-a+(il+l)*dx
           aux=aux*(xm-xl)
           enddo
           ua=ua+AAP(j)*aux*fluxu(il+j)
           va=va+AAP(j)*aux*fluxv(il+j)
          enddo
         ip=ip+1
         write(14,*) 'x(',ip,')=[',xm,'];'
         write(14,*) 'ua(',ip,')=[',ua,'];'
         write(14,*) 'va(',ip,')=[',va,'];'
         write(15,*) 'x(',ip,')=[',xm,'];'
         write(15,*) 'ue(',ip,')=[',ue,'];'
         write(15,*) 've(',ip,')=[',ve,'];'
            aux1=dabs(ua-ue)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
            aux1=dabs(va-ve)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
          enddo
            perrl1(ne)=perrl1(ne)+aux2*dx
            perrl2(ne)=perrl2(ne)+aux3*dx
        enddo
      write(14,3) 
 3    format(&
       t7,' figure(3)'/ &
       t7,' clf'/ &
       t7,'    plot(x,ua);'/ &
       t7,'    hold on'/ &
       t7,' figure(4)'/ &
       t7,' clf'/ &
       t7,'    plot(x,va);'/ &
       t7,'    hold on'/ &
       t7,' ') 


      write(15,4) 
 4    format(& 
       t7,' figure(5)'/ &
       t7,' clf'/ &
       t7,'    plot(x,ue);'/ &      
       t7,'    hold on'/ &
       t7,' figure(6)'/ &
       t7,' clf'/ &
       t7,'    plot(x,ve);'/ &
       t7,'    hold on'/ &
       t7,' ')     
          perrl2(ne)=dsqrt(perrl2(ne))/4.d0
          perrl1(ne)=perrl1(ne)/4.d0
        write(6,*) "post-processing" 
        write(6,5) perrli(ne),perrl1(ne),perrl2(ne)
 5    format('Li=',e12.4,';','L1=', e12.4,';','L2=',e12.4,';')  
      return
      end subroutine post
