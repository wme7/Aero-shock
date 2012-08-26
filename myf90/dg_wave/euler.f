      subroutine euler(tn,j)
c     **********************
c
c       compute u(.,.,j+1) from u(.,.,j)
c       compute v(.,.,j+1) from v(.,.,j)
c
      include 'com' 
c
c      --- evaluating the fluxes
      do 11 i=1,nx-1
        call intflux(i,fu,fv,u,v,j)
        fluxu(i)=fu
        fluxv(i)=fv
   11 continue
c
        call bdryflux(1,fu,fv,u,v,tn,j)
        fluxu(0)=fu
        fluxv(0)=fv
c
        call bdryflux(2,fu,fv,u,v,tn,j)
        fluxu(nx)=fu
        fluxv(nx)=fv
c     
c      --- computing uo at t=tn1
c      --- the interior values
      do 12 i=1,nx
         do 13 l=0,k
c
c         computing u
         aux=0.d00
         do 14 m=0,l-1
            aux=aux+v(i,m+1,j)*(1-(-1)**(m+l))
 14      continue
         u(i,l+1,j+1)=u(i,l+1,j)+c(i)*c(i)*dt/dx*dble(2*l+1)*
     *(-aux - fluxu(i)+fluxu(i-1)*(-1)**l)
c    
c         computing v
         aux=0.d00
         do 15 m=0,l-1
            aux=aux+u(i,m+1,j)*(1-(-1)**(m+l))
 15      continue
         v(i,l+1,j+1)=v(i,l+1,j)+dt/dx*dble(2*l+1)*
     *(-aux - fluxv(i)+fluxv(i-1)*(-1)**l)
 13      continue
 12   continue
c             
c
      end
