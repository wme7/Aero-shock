      subroutine eno
!     ***************
!     post processing           
!     ***************
      use mod_wave
      implicit none
!      double precision,external::func
      double precision::fu,fv,ua,ue,aux,va,ve
      double precision::xic,aux1,xm,aux2,xl,aux3,xi
      double precision,dimension(ks)::AAP
       integer, dimension(ks)::IA
      double precision,dimension(ks,2)::BB
      integer::i,j,m,il,l,ip,ll,lo,in,ii
      integer:: n1,iel,ier
        n1=int((a-rl)/dx)  
!      --- evaluating the fluxes
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


        do i=1,ks
         aux2=0.d0
         aux3=0.d0
         il=-1
          do j=1,ks
           aux=1.d0
           do l=1,j-1
           aux=aux*(j-l)*dx
           enddo
           do l=j+1,ks
           aux=aux*(j-l)*dx
           enddo
!           write(6,*) aux
           AAP(j)=1.d0/aux
          enddo

         xic=-a+(2.d00*i-1.d00)*dx/2.d00

!        evaluating the function `func' at the quadrature points
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
!         write(14,*) 'va(',ip,')=[',va,'];'
         write(15,*) 'x(',ip,')=[',xm,'];'
         write(15,*) 'ue(',ip,')=[',ue,'];'
         write(15,*) 've(',ip,')=[',ve,'];'


            aux1=dabs(ua-ue)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
!            aux1=dabs(va-ve)
!            aux2=aux2+aux1*w(m)
!            aux3=aux3+aux1**2*w(m)
!            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
          enddo
            perrl1(ne)=perrl1(ne)+aux2*dx
            perrl2(ne)=perrl2(ne)+aux3*dx
        enddo
! *************************************
! ENO       
!**************************************
        do i=ks+1,nx-ks
        if (i .le. n1) then
        iel=1
        ier=n1
        else
        iel=n1+1
        ier=nx
        end if

         aux2=0.d0
         aux3=0.d0
         il=i-1
          in=i 
          IA(1)=il
          AAP(1)=(fluxu(il+1)-fluxu(il))/dx

          do j=2,ks
           do l=1,j+1
            BB(l,1)=(fluxu(il+l-1)-fluxu(il+l-2))/dx
           enddo

           lo=1
           ll=2
           do m=1,j-1
           do l=1,j+1-m
           BB(l,ll)=(BB(l+1,lo)-BB(l,lo))/dx/(dble(m+1))
           enddo
           l=lo
           lo=ll
           ll=l

           enddo
              IA(j)=in

           if (dabs(BB(2,lo)).ge.dabs(BB(1,lo))) then
            if (il-1 .lt. iel ) then
           AAP(j)=BB(2,lo)
           in=il+j
           else
           AAP(j)=BB(1,lo)
           il=il-1
           in=il
           end if
           else
            if (il+j .gt. ier) then
           AAP(j)=BB(1,lo)
           il=il-1
           in=il  
            else         
           AAP(j)=BB(2,lo)            
           in=il+j
            end if
           endif


          enddo

         xic=-a+(2.d00*i-1.d00)*dx/2.d00

          do m=1,ks
          ua=0.d0
          ue=0.d0
          va=0.d0
!			         ii=IA(3)

          xm=xic+x(m)*dx/2.d00
!           xm=-a+dble(ii)*dx
           call exac(ue,ve,xm)
          ii=IA(1)
          xi=-a+dble(ii)*dx
          aux=xm-xi
          ua=fluxu(ii)+AAP(1)*aux
!          write(6,*) fluxu(ii)
          do j=2,ks-1
          ii=IA(j)
          xi=-a+dble(ii)*dx
          aux=aux*(xm-xi)
          ua=ua+AAP(j)*aux
          enddo 
         ip=ip+1
         write(14,*) 'x(',ip,')=[',xm,'];'
         write(14,*) 'ua(',ip,')=[',ua,'];'
!         write(14,*) 'va(',ip,')=[',va,'];'
         write(15,*) 'x(',ip,')=[',xm,'];'
         write(15,*) 'ue(',ip,')=[',ue,'];'
         write(15,*) 've(',ip,')=[',ve,'];'
		
            aux1=dabs(ua-ue)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
!			 ii=IA(2)
!			 	aux=ua-fluxu(ii)
!			if (aux .gt. 1.d-8) then
		
!			 write(6,*) "i=",ii,"ua",ua,"ue",ue,"uhat",fluxu(ii-1),fluxu(ii),fluxu(ii+1)
!			 write(6,*) AAP
!			 write(6,*) IA
!			 endif
!            aux1=dabs(va-ve)
!            aux2=aux2+aux1*w(m)
!            aux3=aux3+aux1**2*w(m)
!            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1

          enddo
            perrl1(ne)=perrl1(ne)+aux2*dx
            perrl2(ne)=perrl2(ne)+aux3*dx
        enddo

!***********************************************************************88
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

!           write(6,*) aux
           AAP(j)=1.d0/aux
          enddo

         xic=-a+(2.d00*i-1.d00)*dx/2.d00

!        evaluating the function `func' at the quadrature points
          do m=1,ks
          ua=0.d0
          ue=0.d0
          va=0.d0
          xm=xic+x(m)*dx/2.d00
          call exac(ue,ve,xm)
!            ue=func(xm)
! xic+x(m)*dx/2.d00-t)
! +t**4
!             ue=ue+t**4
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
!         write(14,*) 'va(',ip,')=[',va,'];'
         write(15,*) 'x(',ip,')=[',xm,'];'
         write(15,*) 'ue(',ip,')=[',ue,'];'
         write(15,*) 've(',ip,')=[',ve,'];'
            aux1=dabs(ua-ue)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
!            aux1=dabs(va-ve)
!            aux2=aux2+aux1*w(m)
!            aux3=aux3+aux1**2*w(m)
!            if (aux1 .gt. perrli(ne)) perrli(ne)=aux1
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
       t7,' ')
!       t7,' figure(4)'/ &
!       t7,' clf'/ &
!       t7,'    plot(x,va);'/ &
!       t7,'    hold on'/ &
!       t7,' ') 


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
      end subroutine eno
