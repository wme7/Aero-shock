      subroutine euler(tn,j)
c     **********************
c
c       compute u(.,.,j+1) from u(.,.,j)
c
      include 'com' 
c
c      --- evaluating the fluxes
      do 11 i=1,nx-1
        call intflux(i,fu,u,j)
        fluxu(i)=fu
   11 continue
c
        call bdryflux(1,fu,tn,j)
        fluxu(0)=fu
c
        call bdryflux(2,fu,tn,j)
        fluxu(nx)=fu
c     
c      --- computing uo at t=tn1
c      --- the interior values
      do 12 i=1,nx
         do 13 l=0,k
c
c         computing u
         aux=0.d00
         do 14 m=0,l-1
            aux=aux+u(i,m+1,j)*(1-(-1)**(m+l))
 14      continue
         u(i,l+1,j+1)=u(i,l+1,j)+dt/dx*dble(2*l+1)*
     *(aux - fluxu(i)+fluxu(i-1)*(-1)**l)
c    
 13      continue
 12   continue
c             
c
      end
