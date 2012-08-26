      subroutine errors
!     ***************
!                
       use mod_wave
       implicit none
!       double precision, external::func
       integer::i,m,j
       double precision::aux1,aux2,xic,aux3
       double precision::ua,va,ue,ve,xm
!
!      double precision pl(KK+1,KK+1)
!       double precision ffunc(KK+1),phi(NINT,KK+1)
!      do 1 ne=1,2
         errl2(ne)=0.
         errl1(ne)=0.
         errli(ne)=0.
!    1 continue
!
!      computing the Guassian points
!       call gauleg(-1.d00,1.d00,x,w,k+1)  
!
!      evaluating the Legendre polynomials at
!      the quadrature points
!        do i=1,k+1
!         pl(i,1)=1.d00
!         pl(i,2)=x(i)
!         do j=2,k+1
!          pl(i,j+1)=((2.d00*j-1.d00)*x(i)*pl(i,j)
!     *                -(j-1.d00)*pl(i,j-1))/j
!         enddo
!        enddo
!       write(6,*) x
!
!      computing the moments of the initial data in `func'
!      and storing them in `phi'

        do i=1,nx
         aux2=0.d0
         aux3=0.d0
         xic=-a+(2.d00*i-1.d00)*dx/2.d00

!        evaluating the function `func' at the quadrature points
          do m=1,ks
          ua=0.d0
          va=0.d0
          xm=xic+x(m)*dx/2.d00

          call exac(ue,ve,xm)

!            ue=func(xic+x(m)*dx/2.d00-t)
! +t**4
!         write(6,*) "ue",ue
!
            do j=0,k
            ua=ua+u(i,j+1,0)*pl(m,j+1)
            va=va+v(i,j+1,0)*pl(m,j+1)
            end do
            aux1=dabs(ua-ue)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. errli(ne)) errli(ne)=aux1
            aux1=dabs(va-ve)
            aux2=aux2+aux1*w(m)
            aux3=aux3+aux1**2*w(m)
            if (aux1 .gt. errli(ne)) errli(ne)=aux1
          enddo
            errl1(ne)=errl1(ne)+aux2*dx
            errl2(ne)=errl2(ne)+aux3*dx
        enddo
            errl1(ne)=errl1(ne)/4.d0
            errl2(ne)=dsqrt(errl2(ne))/4.d0
!
        write(6,1) errli(ne),errl1(ne),errl2(ne)
 1    format('Li=',e12.4,';','L1=', e12.4,';','L2=',e12.4,';')
!        i=1/dx/2
!        xm=-a+i*dx
!          call exac(ue,ve,xm)
!        write(6,*) "i=",i,"ue=",ue
      return
      end subroutine errors
