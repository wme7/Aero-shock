      subroutine superconv
!     ***************
!                
      use mod_wave
      implicit none
!      double precision,external::func
      double precision::fu,fv,ua,ue,ve,va
      double precision::xic,aux1
      integer::i,j,m
 
!      --- evaluating the fluxes
        j=0
      do 11 i=1,nx-1
        call intflux(i,fu,fv,u,v,j)
        fluxu(i)=-fv
        fluxv(i)=-fu
   11 continue

        call bdryflux(1,fu,fv,t,j)
        fluxu(0)=-fv
        fluxv(0)=-fu

        call bdryflux(2,fu,fv,t,j)
        fluxu(nx)=-fv
        fluxv(nx)=-fu



!      do 1 ne=1,2
!         serrl1(ne)=0.
         serrli(ne)=0.
!    1 continue

!      computing the Guassian points
!       call gauleg(-1.d00,1.d00,x,w,k+1)  

        do i=1,nx
!         aux2=0.d0
!         xic=-a+(2.d00*i-1.d00)*dx/2.d00
          xic=-a+i*dx
!        write(6,*) xic

!        evaluating the function `func' at the quadrature points
          ua=0
          ue=0

!            ue=func(xic+x(m)*dx/2.d00-t)+t**4
             call exac(ue,ve,xic)
!             ue=func(xic-t)
!             ue=ue+t**4
             ua=fluxu(i)
             va=fluxv(i)

            aux1=dabs(ua-ue)
            if (aux1 .gt. serrli(ne)) serrli(ne)=aux1
            aux1=dabs(va-ve)
            if (aux1 .gt. serrli(ne)) serrli(ne)=aux1
        enddo

        write(6,21) serrli(ne)
 21    format(' Li at node',e12.4,';')  
! ,"L1",serrl1(ne)
      return
      end subroutine superconv
