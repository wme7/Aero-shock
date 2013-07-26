! Use discontinuous Galerkin, polynomial truncation for the nonlinear
! terms, to compute 1D Boltzmann-BGK Equation (conservation laws) - July 18, 2013
! The main modification is to add an outside loop with do is=1, nv (isigma from 1 to nv) enddo
program main
  use data_module
  integer kcount, lll, kk0, i
  real burgex_val
  character*32 yc,yc1

  call setup

  open(109,file='weno_pt'//'_p'//char(48+mp)//'_flux'//char(kflux+48)//'.plt')
  indexmin=1000
  indexmax=0
  kcount=0
  call initdata
  call init
  ! call coeff
  ! pause
  ! call smooth_coef
  ! begin time iteration

  ! determine dt
  call setdt

  do
     if((t.ge.tprint-1.e-10).or.(kcount.ge.kcmax)) exit

     if(kcount/1*1.eq.kcount) then
        ! write(6,*) kcount,' t= ',t
        ! pause
     end if

     do is=1,nv
        ! Runge-Kutta for all discrete velocity v_{\sigma}
        call rk(is)
        do lll=1,index
           write(109,*) is,x(indx1(lll)),t
        enddo
     enddo

     call macrop
     print *,kcount,t,dt

     t=t+dt
     kcount=kcount+1
     call setdt
  end do
  call macrop

  ! print and compute errors
  do is=1,nv
     call limit(is)
     do lll=1,index
        write(109,*) is,x(indx1(lll)),t
     enddo
  enddo

  call tec_output

  write(*,*) n, indexmin,indexmax
  write(*,*) 'points:',index,(x(indx1(i)),i=1,index)

102 format(i6,1x,3('&',1x, es12.2e2,1x,'&' 1x,f8.2 ,1x)&
       ,'&',1x, i6,1x,'&' 1x,i6 ,1x&
       ,'\\',1x,'\hline')
103 format(i6,1x,3('&',1x,es12.2E2,1x,'&',1x),'&',1x, i6,1x&
       ,'&' 1x,i6 ,1x,'\\',1x,'\hline')
124 format(i6,1x,3( es12.5e2,1x))
125 format(i6,1x,3( f12.5,1x))
  stop
end program main


!!-------------------------------------------
!!
!!-------------------------------------------
subroutine tec_output
  use data_module
  real*8 rho,uuu,mx,me,tt,pp

  !write sol
  open(1,file='sol_cell_avg_'//'_p'//char(48+mo)//'_flux'//char(kflux+48)//'.dat')
  WRITE(1,*) 'VARIABLES = "x","rho","p","t","u" '

  do i=0,n
     rho=0.
     mx=0.
     me=0.
     do is=1,nv
        rho = rho + cc(is) * u(0,i,is)
        mx = mx + cc(is) * u(0,i,is) * vv(is)
        me  = me + cc(is) * u(0,i,is) * (0.5 * vv(is) * vv(is))
     enddo

     uuu = mx/rho
     tt = 4.0d0*me/rho - 2.0d0*uuu**2
     pp = 0.5d0*rho*tt

     write(1,123) x(i),rho,pp,tt,uuu
  enddo
  close(1)

123 format(4(1x,f16.6))
end subroutine tec_output


!!-------------------------------------------
!!
!!-------------------------------------------
subroutine setup
  use data_module
  ! set up the problem and the initial condition

  pi=4.0*atan(1.0)
  sr=sqrt(5.0)/10.0

  open(unit= 24, file="proj.in", form="FORMATTED", action="READ")
  read(24,nml=proj_list)

  mp=mo-1
  if(xmmm.lt.0.10) then
     mm=0
  elseif(xmmm.lt.2.0) then
     mm=1
  elseif(xmmm.lt.20.0) then
     mm=2
  elseif(xmmm.lt.60.) then
     mm=3
  elseif(xmmm.lt.145.0) then
     mm=4
  else
     mm=5
  endif

  kcmax=1000000

  return
end subroutine setup

!!-------------------------------------------
!! Init
!!-------------------------------------------
subroutine init
  use data_module

  dimension a(0:10), temp(0:n)
  ! set up the initial condition
  ! u0(x)=.5+sin(pi*x)
  gh = (/-10.1591092462,-9.52090367701,-8.9923980014, &
                    -8.52056928412,-8.08518865425,-7.6758399375,-7.2862765944, &
                    -6.91238153219,-6.55125916706,-6.20077355799,-5.85929019639, &
                    -5.52552108614,-5.19842653458,-4.87715007747,-4.56097375794, &
                    -4.24928643596,-3.94156073393,-3.63733587617,-3.33620465355, &
                    -3.03780333823,-2.74180374807,-2.44790690231,-2.15583787123, &
                    -1.86534153123,-1.57617901198,-1.28812467487,-1.00096349956, &
                    -0.714488781673,-0.428500064221,-0.142801238703,0.142801238703, &
                    0.428500064221,0.714488781673,1.00096349956,1.28812467487, &
                    1.57617901198,1.86534153123,2.15583787123,2.44790690231, &
                    2.74180374807,3.03780333823,3.33620465355,3.63733587617, &
                    3.94156073393,4.24928643596,4.56097375794,4.87715007747, &
                    5.19842653458,5.52552108614,5.85929019639,6.20077355799, &
                    6.55125916706,6.91238153219,7.2862765944,7.6758399375, &
                    8.08518865425,8.52056928412,8.9923980014,9.52090367701, &
                    10.1591092462/)
            w = (/1.10958724797e-45,2.43974758815e-40, &
                    3.77162672712e-36,1.33255961176e-32,1.71557314767e-29, &
                    1.02940599717e-26,3.34575695575e-24,6.5125672575e-22, &
                    8.15364047302e-20,6.92324790958e-18,4.15244410969e-16, &
                    1.81662457626e-14,5.94843051606e-13,1.48895734906e-11, &
                    2.89935901281e-10,4.45682277523e-9,5.47555461928e-8, &
                    5.4335161342e-7,4.39428693627e-6,0.0000291874190416, &
                    0.000160277334682,0.000731773556966,0.00279132482895, &
                    0.00893217836031,0.0240612727661,0.0547189709322,0.105298763698, &
                    0.171776156919,0.237868904959,0.279853117523,0.279853117523, &
                    0.237868904959,0.171776156919,0.105298763698,0.0547189709322, &
                    0.0240612727661,0.00893217836031,0.00279132482895, &
                    0.000731773556966,0.000160277334682,0.0000291874190416, &
                    4.39428693627e-6,5.4335161342e-7,5.47555461928e-8, &
                    4.45682277523e-9,2.89935901281e-10,1.48895734906e-11, &
                    5.94843051606e-13,1.81662457626e-14,4.15244410969e-16, &
                    6.92324790958e-18,8.15364047302e-20,6.5125672575e-22, &
                    3.34575695575e-24,1.02940599717e-26,1.71557314767e-29, &
                    1.33255961176e-32,3.77162672712e-36,2.43974758815e-40, &
                    1.10958724797e-45/)

  if (ighq .eq. 1) go to 10
  do is=1,nv
     cc(is) =w(is)
     vv(is) = gh(is)
  end do
  go to 30
10 continue
  v1  = -10
  v2  = 10
  dv  = (v2-v1)/(nv-1.)
  do is =1,nv
     vv(is) = v1 + dv * (is-1.)
  end do
  do is = 2, (nv-1)
     cc(is) = 64./45. * dv
     if(mod(is,4) .eq.1) cc(is) = 28./45. * dv
     if (mod(is,4) .eq.3) cc(is) = 24./45. * dv
  end do
  cc(1)    = 14./45. * dv
  cc(nv)   = cc(1)
30 continue
  pi=4.0*atan(1.0)
  xleft=  .0
  xright= 1.0
  xlen=xright-xleft

  dxuni=xlen/n
  dxmin=1.e10
  do i=0,n
     temp(i)=xleft+dxuni*i
  enddo

  do i=1,n
     x(i)=(temp(i)+temp(i-1))/2.0
     dx(i)=temp(i)-temp(i-1)
     dxmin=min(dxmin,dx(i))
  enddo
  do i=0,md
     x(-i)=x(n-i)-xlen
     dx(-i)=dx(n-i)
     x(n+i)=x(i)+xlen
     dx(n+i)=dx(i)
  enddo
  !Initial condition for shock tube-Density: rho; Mean Velocity: uxm; Pressure: pre; Temperature:tem;
  ! u(kk,i,is) will be the unknowns, the degree of freedom

  rl=0.125
  ul=0.
  pl=0.1
  rr=1.0
  ur=0.
  pr=1.

  do i = 1, n
     if (x(i) .le. 0.5)then
        do kk=0,mp
           rki(kk,i) = rl
           uki(kk,i) = ul
           pki(kk,i) = pl
           tki(kk,i) = 2.0d0*pl/rl
        end do
     else
        do kk=0,mp
           rki(kk,i) = rr
           uki(kk,i) = ur
           pki(kk,i) = pr
           tki(kk,i) = 2.0d0*pr/rr
        end do
     end if
     do kk=0,mp
        do is = 1, nv
           pres       = (vv(is)-uki(kk,i))**2/tki(kk,i)
           u(kk,i,is)   = rki(kk,i)*exp(-pres)/sqrt(pi*tki(kk,i))
        end do
     end do
  end do

  do is=1,nv
     do i=0,n+1
        do kk=0,mp
           a(kk)=(u(kk,i,is)*fle(kk,gau(1,1)) + u(kk,i,is)*fle(kk,gau(2,1)))*gau(1,2)&
                +(u(kk,i,is)*fle(kk,gau(3,1)) + u(kk,i,is)*fle(kk,gau(4,1)))*gau(3,2)&
                +(u(kk,i,is)*fle(kk,gau(5,1)) + u(kk,i,is)*fle(kk,gau(6,1)))*gau(6,2)

           ! a(kk)=(u0(x(i)+gau(1,1)*dx(i),is)*fle(kk,gau(1,1))&
           !             +u0(x(i)+gau(2,1)*dx(i),is)*fle(kk,gau(2,1)))*gau(1,2)&
           !             +(u0(x(i)+gau(3,1)*dx(i),is)*fle(kk,gau(3,1))&
           !             +u0(x(i)+gau(4,1)*dx(i),is)*fle(kk,gau(4,1)))*gau(3,2)&
           !             +(u0(x(i)+gau(5,1)*dx(i),is)*fle(kk,gau(5,1))&
           !             +u0(x(i)+gau(6,1)*dx(i),is)*fle(kk,gau(6,1)))*gau(6,2)
        enddo
        ! take care of the mass matrix
        do kk=0,mp
           if(kk.eq.0) then
              u(kk,i,is)=ai(kk,kk)*a(kk)
           else
              u(kk,i,is)=ai(kk,kk)*a(kk)
           endif
        enddo
     enddo
  enddo
  kcount=0
  t=0.

  return
end subroutine init

!!-------------------------------------------
!!
!!-------------------------------------------
!function u0(x00)
!  pi=4.0*atan(1.0)
!  u0=sin(pi*x00)+0.5
!  return
!end function u0

!!-------------------------------------------
!!
!!-------------------------------------------
subroutine res(is)
  use data_module

  dimension h(0:md,0:n+1,nv),flx(0:n+1),a(0:md),b(0:md)&
       ,un(-1:n+1),up(-1:n+1)

  ! compute the residue
  ! Boltzmann-BGK equation

  ! compute the maximum f'(u)
  sr=sqrt(5.0)/10.0
  sr7=sqrt(21.0)/14.0
  if(kflux.eq.2) then
     amax=0.0
     do i=1,n
        amax=max(amax,abs(vv(is)))
     enddo
  end if

  ! amax=1.0
  ! compute the contribution from the cell boundary
  do  i=0,n+1
     un(i)=eval(u(0,i,is),mp,0.5)
     up(i-1)=eval(u(0,i,is),mp,-0.5)
  enddo

  do i=0,n
     if(kflux.eq.1) then
        ! Roe flux
        !        if(un(i)+up(i).gt.0.) then
        !           flx(i)=vv(is)*un(i)
        !           flx(i)=flux(un(i))
        flx(i)=0.5*(vv(is)*up(i) + vv(is)*un(i) - amax*(up(i)-un(i)))
        !
        !           flx(i)=vv(is)*up(i)
        !           flx(i)=flux(up(i))
        ! flx(i)=un(i)
     else if(kflux.eq.2) then
        ! global LF flux
        flx(i)=0.5*(vv(is)*up(i) + vv(is)*un(i) - amax*(up(i)-un(i)))
        !        flx(i)=0.5*(flux(up(i))+flux(un(i))-amax*(up(i)-un(i)))
     else
        ! local LF flux
        amax=max(abs(un(i)),abs(up(i)))
        flx(i)=0.5*(vv(is)*up(i) + vv(is)*un(i) - amax*(up(i)-un(i)))
        !        flx(i)=0.5*(flux(up(i))+flux(un(i))-amax*(up(i)-un(i)))
     end if
  enddo

  do i=1,n
     do k=0,mp
        h(k,i,is)=flx(i-1)*bl(k)-flx(i)*br(k)
     enddo
  enddo

  ! compute the contribution from the volume integral by Gauss-Lobatto quadrature
  do i=1,n
     if (mp.eq.1) then
        h(1,i,is)=h(1,i,is)+(vv(is)*up(i-1) + vv(is)*eval(u(0,i,is),mp,0.)*4.0 + vv(is)*un(i))/6.0
        !        h(1,i,is)=h(1,i,is)+(flux(up(i-1))+flux(eval(u(0,i,is),mp,0))*4.0+flux(un(i)))/6.0
     else if(mp.eq.2)  then
        do kk=1,mp
           h(kk,i,is)=h(kk,i,is)&
                +( vv(is)*up(i-1)*fled(kk,-0.5) + vv(is)*un(i)*fled(kk,0.5)&
                +( vv(is)*eval(u(0,i,is),mp,-sr)*fled(kk,-sr)&
                + vv(is)*eval(u(0,i,is),mp,sr)*fled(kk,sr) )*5.0 )/12.0
           !          h(kk,i,is)=h(kk,i,is)&
           !               +(flux(up(i-1))*fled(kk,-0.5)+flux(un(i))*fled(kk,0.5)&
           !               +(flux(eval(u(0,i,is),mp,-sr))*fled(kk,-sr)&
           !               +flux(eval(u(0,i,is),mp,sr))*fled(kk,sr))*5.0)/12.0
        enddo
     else if(mp.eq.3) then
        do kk=1,mp
           h(kk,i,is)=h(kk,i,is)&
                +( (vv(is)*up(i-1)*fled(kk,-0.5) + vv(is)*un(i)*fled(kk,0.5) )*9.0&
                +( vv(is)*eval(u(0,i,is),mp,-sr7)*fled(kk,-sr7)&
                + vv(is)*eval(u(0,i,is),mp,sr7)*fled(kk,sr7) )*49.0&
                + vv(is)*eval(u(0,i,is),mp,0.)*fled(kk,0.)*64.0)/180.0
           !           h(kk,i,is)=h(kk,i,is)&
           !                +((flux(up(i-1))*fled(kk,-0.5)+flux(un(i))*fled(kk,0.5))*9.0&
           !                +(flux(eval(u(0,i,is),mp,-sr7))*fled(kk,-sr7)&
           !                +flux(eval(u(0,i,is),mp,sr7))*fled(kk,sr7))*49.0&
           !                +flux(eval(u(0,i,is),mp,0))*fled(kk,0.0)*64.0)/180.0
        enddo
     else
        do kk=1,mp
           h(kk,i,is)=h(kk,i,is)&
                +(vv(is)*up(i-1)*fled(kk,-.5) + vv(is)*un(i)*fled(kk,.5))*gauss(1,2)&
                +( vv(is)*eval(u(0,i,is),mp,gauss(3,1))*fled(kk,gauss(3,1))&
                + vv(is)*eval(u(0,i,is),mp,gauss(4,1))*fled(kk,gauss(4,1)) )*gauss(3,2)&
                +( vv(is)*eval(u(0,i,is),mp,gauss(5,1))*fled(kk,gauss(5,1))&
                + vv(is)*eval(u(0,i,is),mp,gauss(6,1))*fled(kk,gauss(6,1)) )*gauss(5,2)
           !           h(kk,i,is)=h(kk,i,is)&
           !                +(flux(up(i-1))*fled(kk,-.5)+flux(un(i))*fled(kk,.5))*gauss(1,2)&
           !               +(flux(eval(u(0,i,is),mp,gauss(3,1)))*fled(kk,gauss(3,1))&
           !               +flux(eval(u(0,i,is),mp,gauss(4,1)))*fled(kk,gauss(4,1)))*gauss(3,2)&
           !               +(flux(eval(u(0,i,is),mp,gauss(5,1)))*fled(kk,gauss(5,1))&
           !               +flux(eval(u(0,i,is),mp,gauss(6,1)))*fled(kk,gauss(6,1)))*gauss(5,2)

        enddo
     endif
  enddo

  ! take care of the mass matrix
  do i=1,n
     do kk=0,mp
        hg(kk,i,is)=ai(kk,kk)*h(kk,i,is)/dx(i)
        !     if(kk.eq.4) hg(kk,i,is)=0.0
     enddo
  enddo
  return
end subroutine res

!!-------------------------------------------
!! Find macroscopic properties and equilibrium distribution
!!-------------------------------------------
subroutine macrop
  use data_module
  do i = 1, n
     do kk=0,mp
        skr = 0.
        sku = 0.
        ske = 0.
        do is = 1, nv
           skr = skr + cc(is) * u(kk,i,is)
           sku = sku + cc(is) * u(kk,i,is) * vv(is)
           ske = ske + cc(is) * u(kk,i,is) * (0.5 * vv(is) * vv(is))
        enddo
        rki(kk,i) = skr
        uki(kk,i) = sku/skr
        eki(kk,i) = ske
        tki(kk,i) = 4.0d0*eki(kk,i)/rki(kk,i) - 2.0d0*uki(kk,i)*uki(kk,i)
        pki(kk,i) = 0.5d0*rki(kk,i)*tki(kk,i)

     enddo
  enddo
  do is = 1, nv
     do i = 1, n
        do kk=0,mp
           pres  = (vv(is)-uki(kk,i))**2/tki(kk,i)
           ueq(kk,i,is) = rki(kk,i)*exp(-pres)/sqrt(pi*tki(kk,i))
        enddo
     enddo
  enddo

end subroutine macrop
!!-------------------------------------------
!!Calculate macroscopic properties and equilibrium distribution
!!-------------------------------------------

!!-------------------------------------------
!!
!!-------------------------------------------
subroutine setdt
  use data_module

  ! set up dt
  amax=0.
  do iv=1,nv
     amax=max(amax,abs(vv(iv)))
  enddo

  if(mp.eq.3) then
     if(mp.le.2.0+1.0e-5) then
        rr=1.0
     else
        rr=(mp+1.0)/mt
     endif
  else
     if(mp.le.3.0+1.0e-5) then
        rr=1.0
     else
        rr=(mp+1.0)/mt
     endif
  endif
  dt=cflc*dxmin**rr/amax
  if(t+dt.gt.tprint) dt=tprint-t
  return
end subroutine setdt

!!-------------------------------------------
!!
!!-------------------------------------------
subroutine bc(is)
  use data_module

  ! set up the boundary condition keep constant state extrapolation
  ! periodic
  do k=0,mp
     do j=0,md
        u(k,-j,is)=u(k,1,is)
        u(k,n+j+1,is)=u(k,n,is)
     enddo
  enddo

  return
end subroutine bc

!!-------------------------------------------
!!
!!-------------------------------------------
!!function flux(x)
!  ! the function of flux
!  flux=x*vv(is)
!  return
!end function flux

!!-------------------------------------------
!!
!!-------------------------------------------
function burgex(te,xe)
  use data_module
  ! prepare two matrices for the Burgers' equation exact solution
  dxx=1.0/199.0
  do  i=1,200
     xx(i)=2.*dxx*float(i-1)-1.0
     uu(i)=.5+sin(pi*xx(i))
  enddo

  xt0=(2*xe-1.0)-.5*te
  if(xt0.gt.1.) xt0=xt0-2.*int((xt0+1.)/2.)
  if(xt0.lt.-1.) xt0=xt0+2.*int((1.-xt0)/2.)
  ay=abs(xt0)
  i0=1
  do i=1,200
     xt=xx(i)+uu(i)*te
     if(ay.lt.xt) goto 7
     i0=i
  enddo
7 us=uu(i0)
  un=us
  do i=1,400
     us=un
     x0=ay-us*te
     un=us-(us-sin(pi*x0))/(1.+pi*cos(pi*x0)*te)
     if(abs(un-us).lt.1.e-15) goto 15
  enddo
  ! write(6,*) 'did not converge at (xe,te) = ',xe,te
  ! write(6,*) ' final newton residue = ', un-us
15 y=sign(1.0,xt0)*un
  if(abs(ay-1.0).lt.1.e-15) y=0.0
  if(abs(te).lt.1.e-15) y=sin(pi*xe)
  burgex=y*0.5+.250

  return
end function burgex

!!-------------------------------------------
!!
!!-------------------------------------------
function eval(a,m,x0)
  use data_module

  ! evaluate the value of the polynomial of degree m with coefficient a,
  ! at the location x=x0
  ! a(0)+a(1)*x0+...

  dimension a(0:md)

  eval=0.0
  do i=m,0,-1
     eval=eval+fle(i,x0)*a(i)
  enddo

  return
end function eval


!!-------------------------------------------
!!
!!-------------------------------------------
function fle(k,x)
  ! the function of Legendre polynomial
  if(k.eq.0) then
     fle=1.0
  elseif (k.eq.1) then
     fle=x
  elseif(k.eq.2) then
     fle=x*x-1./12.0
  elseif (k.eq.3) then
     fle=x*x*x-0.15*x
  elseif(k.eq.4) then
     fle=(x**2-3./14.0)*x*x+3.0/560.0
  endif
  return
end function fle

!!-------------------------------------------
!!
!!-------------------------------------------
function fled(k,x)
  ! the function of derivative of Legendre polynomial
  if (k.eq.1) then
     fled=1
  elseif(k.eq.2) then
     fled=2.0*x
  elseif (k.eq.3) then
     fled=3.0*x*x-0.15
  elseif(k.eq.4) then
     fled=(4.0*x*x-3.0/7.0)*x
  endif
  return
end function fled

!!-------------------------------------------
!!
!!-------------------------------------------
subroutine initdata
  ! use dfport
  ! use dflib
  use data_module

  dimension a(0:md),b(0:md),c(0:md),aic(0:20,0:20,0:10),temp(0:n)

  ! set up the necessary data before setting the initial condition
  do k=0,mp
     bl(k)=fle(k,-0.5)
     br(k)=fle(k,0.5)
  enddo

  ! read in the inverse of the mass matrix:
  aic(0,0,0)=1.
  aic(0,0,1)=1.
  aic(0,1,1)=0.0
  aic(1,1,1)=12.0
  aic(0,0,2)=1.0
  aic(0,1,2)=0.0
  aic(0,2,2)=0.0
  aic(1,1,2)=12.0
  aic(1,2,2)=0.0
  aic(2,2,2)=180.0
  aic(0,0,3)=1.0
  aic(0,1,3)=0.0
  aic(0,2,3)=0.0
  aic(0,3,3)=0.0
  aic(1,1,3)=12.0
  aic(1,2,3)=0.0
  aic(1,3,3)=0.0
  aic(2,2,3)=180.0
  aic(2,3,3)=0.0
  aic(3,3,3)=2800.0
  aic(0,0,4)=1
  aic(0,1,4)=0.0
  aic(0,2,4)=0
  aic(0,3,4)=0.0
  aic(0,4,4)=0
  aic(1,1,4)=12.0
  aic(1,2,4)=0.0
  aic(1,3,4)=0
  aic(1,4,4)=0.0
  aic(2,2,4)=180
  aic(2,3,3)=0.0
  aic(2,4,4)=0
  aic(3,3,4)=2800.0
  aic(3,4,4)=0.0
  aic(4,4,4)=44100.0

  do i0=1,4
     do i01=0,i0
        do i02=i01+1,i0
           aic(i02,i01,i0)=aic(i01,i02,i0)
        enddo
     enddo
  enddo

  mpp=max(2,mp/2*2)
  do j=0,mp
     do i=0,mp
        ai(i,j)=aic(i,j,mp)
     enddo
  enddo

  aic(0,0,0)=1.
  aic(0,0,1)=1.
  aic(0,1,1)=0.0
  aic(1,1,1)=12.0
  aic(0,0,2)=2.250
  aic(0,1,2)=0.0
  aic(0,2,2)=-15.0
  aic(1,1,2)=12.0
  aic(1,2,2)=0.0
  aic(2,2,2)=180.0
  aic(0,0,3)=2.25
  aic(0,1,3)=0.0
  aic(0,2,3)=-15.0
  aic(0,3,3)=0.0
  aic(1,1,3)=75.0
  aic(1,2,3)=0.0
  aic(1,3,3)=-420.0
  aic(2,2,3)=180.0
  aic(2,3,3)=0.0
  aic(3,3,3)=2800.0
  aic(0,0,4)=225.0/64.0
  aic(0,1,4)=0.0
  aic(0,2,4)=-525.0/8.0
  aic(0,3,4)=0.0
  aic(0,4,4)=945.0/4.0
  aic(1,1,4)=75.0
  aic(1,2,4)=0.0
  aic(1,3,4)=-420.0
  aic(1,4,4)=0.0
  aic(2,2,4)=2205.0
  aic(2,3,4)=0.0
  aic(2,4,4)=-9450.0
  aic(3,3,4)=2800.0
  aic(3,4,4)=0.0
  aic(4,4,4)=44100.0
  do i0=1,4
     do i01=0,i0
        do i02=i01+1,i0
           aic(i02,i01,i0)=aic(i01,i02,i0)
        enddo
     enddo
  enddo
  mpp=max(2,mp/2*2)

  do j=0,mp
     do i=0,mp
        aii(i,j)=aic(i,j,mp)
     enddo
  enddo

  ! the points of 10th order Gauss-Lobatto quadrature
  gau(1,1)=-0.5
  gau(2,1)=0.5
  gau(3,1)=-sqrt(147.-42.*sqrt(7.))/42.0
  gau(4,1)= sqrt(147.-42.*sqrt(7.))/42.0
  gau(5,1)=-sqrt(147.+42.*sqrt(7.))/42.0
  gau(6,1)= sqrt(147.+42.*sqrt(7.))/42.0

  ! coefficients of 10th order gau-Lobatto quadrature
  gau(1,2)=1.0/30.0
  gau(2,2)=gau(1,2)
  gau(5,2)=(-7.+5.*sqrt(7.))*sqrt(7.)*(7.+sqrt(7.))/840.0
  gau(6,2)=gau(5,2)
  gau(3,2)=(7.+5.*sqrt(7.))*sqrt(7.)/(7.+sqrt(7.))/20.0
  gau(4,2)=gau(3,2)
  if(mp.eq.1) then
     ! the points of 4h order Gauss-Lobatto quadrature
     gauss(1,1)=-0.5
     gauss(2,1)=0.5
     gauss(3,1)=0.0
     ! coefficients of 4h order Gauss-Lobatto quadrature
     gauss(1,2)=1.0/6.0
     gauss(2,2)=gauss(1,2)
     gauss(3,2)=2.0/3.0

  endif
  if(mp.eq.2) then
     ! the points of 6th order Gauss-Lobatto quadrature
     gauss(1,1)=-0.5
     gauss(2,1)=0.5
     gauss(3,1)=-sqrt(5.)/10.0
     gauss(4,1)= sqrt(5.0)/10.0
     ! coefficients of 6th order Gauss-Lobatto quadrature
     gauss(1,2)=1.0/12.0
     gauss(2,2)=gauss(1,2)
     gauss(3,2)=5.0/12.0
     gauss(4,2)=gauss(3,2)
     ! gauss(5,2)=64.0/180.0
  endif
  if(mp.eq.3) then
     ! the points of 8th order Gauss-Lobatto quadrature
     gauss(1,1)=-0.5
     gauss(2,1)=0.5
     gauss(3,1)=-sqrt(21.)/14.0
     gauss(4,1)= sqrt(21.0)/14.0
     gauss(5,1)=0.0
     ! coefficients of 8th order Gauss-Lobatto quadrature
     gauss(1,2)=1.0/20.0
     gauss(2,2)=gauss(1,2)
     gauss(3,2)=49.0/180.0
     gauss(4,2)=gauss(3,2)
     gauss(5,2)=64.0/180.0
  endif
  if(mp.eq.4) then
     ! the points of 10th order Gauss-Lobatto quadrature
     gauss(1,1)=-0.5
     gauss(2,1)=0.5
     gauss(5,1)=-sqrt(147.-42.*sqrt(7.))/42.0
     gauss(6,1)= sqrt(147.-42.*sqrt(7.))/42.0
     gauss(3,1)=-sqrt(147.+42.*sqrt(7.))/42.0
     gauss(4,1)= sqrt(147.+42.*sqrt(7.))/42.0

     ! coefficients of 10th order Gauss-Lobatto quadrature
     gauss(1,2)=1.0/30.0
     gauss(2,2)=gauss(1,2)
     gauss(3,2)=(-7.+5.*sqrt(7.))*sqrt(7.)*(7.+sqrt(7.))/840.0
     gauss(4,2)=gauss(3,2)
     gauss(5,2)=(7.+5.*sqrt(7.))*sqrt(7.)/(7.+sqrt(7.))/20.0
     gauss(6,2)=gauss(5,2)
  endif

11 format(6e13.5)
1001 format(5(f12.8,2x))

  return
end subroutine initdata

!!-------------------------------------------
!!
!!-------------------------------------------
subroutine limit(is)
  use data_module

  dimension b(0:md),du(-1:n+1,nv),am(10),temp(10)
  !     xmmm=-2.0/3.0
  !     limit the solution
  index=0
  if(xm2.lt.-1.e-4.or.mp.eq.0) then
     !     do nothing
     return
  else
     !     compute the first order difference of the mean for limiting

     do i=-1,n
        du(i,is)=u(0,i+1,is)-u(0,i,is)
     enddo
     do i=0,n
        xmm=xmmm*dx(i)**2

        !     compute the cell boundary values:
        uleft=eval(u(0,i,is),mp,0.5)
        uright=eval(u(0,i,is),mp,-0.5)

        !     limit the values uleft and uright
        am(1)=uleft-u(0,i,is)
        am(2)=du(i-1,is)
        am(3)=du(i,is)
        sg=sign(1.,am(1))
        amc=am(1)
        if(abs(amc).le.xmm)  goto 31
        do j=2,3
           if(amc*am(j).lt.0.) then
              !     amc=0.
              goto 330
           else
              if(abs(am(j)).lt.abs(amc)) then
                 !     amc=am(j)
                 goto 330
              end if
           end if
        enddo

31      continue
        ! if(abs(amc-am(1)).gt.1e-6) goto 100
        ! if(amc.ne.am(1)) goto 100
        am(1)=u(0,i,is)-uright
        am(2)=du(i-1,is)
        am(3)=du(i,is)
        amc=am(1)
        if(abs(amc).le.xmm) cycle
        do j=2,3
           if(amc*am(j).lt.0.) then
              !     amc=0.
              goto 330
           else
              if(abs(am(j)).le.abs(amc)) then
                 !     amc=am(j)
                 goto 330
              end if
           end if
        enddo
        cycle
        !     100  if(abs(amc-am(1)).gt.1e-6)  then
330     index=index+1
        call wenorecon(i,is)
        indexmin=min(indexmin,index)
        indexmax=max(indexmax,index)
        indx1(index)=i
     end do
  end if

  if(index.gt.n/20) then
     !     write(*,*) index
  endif
  return
end subroutine limit

!!-------------------------------------------
!!
!!-------------------------------------------
subroutine rk(is)
  use data_module
  dimension temp1(0:mp,n,3)
  ! third order Runge-Kutta
  ! set up boundary condition for u
  if(mt.eq.3) then
     call bc(is)

     ! limit the solution
     call limit(is)
     call bc(is)

     do i=1,n
        !  uold(i)=u(0,i,is)
        do k=0,mp
           v(k,i,is)=u(k,i,is)
        enddo
     enddo

     ! compute the residue
     call res(is)
     ! call  wenorecon
     do i=1,n
        do k=0,mp
           u(k,i,is)=v(k,i,is)+dt*hg(k,i,is)
           temp1(k,i,1)=u(k,i,is)
        enddo
     enddo

     ! call bc
     ! call limit
     call bc(is)
     call res(is)
     ! call  wenorecon

     do i=1,n
        do k=0,mp
           u(k,i,is)=0.75*v(k,i,is)+0.25*(temp1(k,i,1)+dt*hg(k,i,is))
           temp1(k,i,2)=u(k,i,is)
        enddo
     enddo

     ! call bc
     ! call limit
     call bc(is)
     call res(is)
     ! call  wenorecon

     do i=1,n
        do k=0,mp
           u(k,i,is)=(v(k,i,is)+2.0*(temp1(k,i,2)+dt*hg(k,i,is)))/3.0
        enddo
     enddo
  else
     ! set up boundary condition for u
     call bc(is)

     ! limit the solution
     call limit(is)
     call bc(is)

     do i=1,n
        ! uold(i)=u(0,i)
        do k=0,mp
           v(k,i,is)=u(k,i,is)
        enddo
     enddo

     ! compute the residue
     call res(is)
     ! call  wenorecon
     do i=1,n
        do k=0,mp
           u(k,i,is)=v(k,i,is)+0.5*dt*hg(k,i,is)
           temp1(k,i,1)=u(k,i,is)
        enddo
     enddo

     call bc(is)
     call limit(is)
     call bc(is)
     call res(is)
     ! call  wenorecon

     do i=1,n
        do k=0,mp
           u(k,i,is)=v(k,i,is)+0.5*dt*hg(k,i,is)
           temp1(k,i,2)=u(k,i,is)
        enddo
     enddo

     call bc(is)
     call limit(is)
     call bc(is)
     call res(is)
     ! call  wenorecon

     do i=1,n
        do k=0,mp
           u(k,i,is)=v(k,i,is)+dt*hg(k,i,is)
           temp1(k,i,3)=u(k,i,is)
        enddo
     enddo

     call bc(is)
     call limit(is)
     call bc(is)
     call res(is)
     ! call  wenorecon
     do i=1,n
        do k=0,mp
           u(k,i,is)=(-v(k,i,is)+temp1(k,i,1)+2*temp1(k,i,2)&
                +temp1(k,i,3)+0.5*dt*hg(k,i,is))/3.
        enddo
     enddo
  endif
  return
end subroutine rk

!!-------------------------------------------
!! base function for interpolant of WENO reconstruction
!!-------------------------------------------
function qll(x0,k1,k2,l,i)
  use data_module
  ! dimension a(-10:10)
  x00=x(i)+x0*dx(i)
  qll=0.0
  do kk=l,k2
     t1=1.0
     do j=k1-1,k2
        if(j.ne.kk) t1=t1*(x(kk+i)+0.5*(dx(kk+i)-dx(j+i))-x(j+i))
     enddo
     t2=0.0
     do j=k1-1,k2
        if(j.ne.kk) then
           t3=1.0
           do j1=k1-1,k2
              if(j1.ne.j.and.j1.ne.kk) t3=(x00-x(j1+i)-0.5*dx(j1+i))*t3
           enddo
           t2=t2+t3
        endif
     enddo
     qll=qll+t2/t1
  enddo
  qll=qll*dx(i+l)
  return
end function qll

!!-------------------------------------------
!! polynomial function of reconstruction with cells: kk+k1---kk+k2
!!-------------------------------------------
function pll(x0,k1,k2,i,is)
  use data_module
  ! dimension a(-10:10)
  pll=0.0
  do l=k1,k2
     pll=pll+qll(x0,k1,k2,l,i)*u(0,l+i,is)
  enddo
  return
end function pll

!!-------------------------------------------
!! set the coefficent of function pll, and compute the smooth indicators (is indicates isigma)
!!-------------------------------------------
subroutine smooth(i,is,s)
  use data_module
  dimension s(0:10,nv),aa(0:5,2)
  if(mp.eq.1) then
     aa(1,1)=(pll(gau(1,1),-1,0,i,is)*gau(1,1)&
          +pll(gau(2,1),-1,0,i,is)*gau(2,1))*gau(1,2)&
          +(pll(gau(3,1),-1,0,i,is)*gau(3,1)&
          +pll(gau(4,1),-1,0,i,is)*gau(4,1))*gau(3,2)&
          +(pll(gau(5,1),-1,0,i,is)*gau(5,1)&
          +pll(gau(6,1),-1,0,i,is)*gau(6,1))*gau(5,2)
     s(1,is)=aa(1,1)**2
     aa(1,1)=(pll(gau(1,1),0,1,i,is)*gau(1,1)&
          +pll(gau(2,1),0,1,i,is)*gau(2,1))*gau(1,2)&
          +(pll(gau(3,1),0,1,i,is)*gau(3,1)&
          +pll(gau(4,1),0,1,i,is)*gau(4,1))*gau(3,2)&
          +(pll(gau(5,1),0,1,i,is)*gau(5,1)&
          +pll(gau(6,1),0,1,i,is)*gau(6,1))*gau(5,2)
     s(2,is)=aa(1,1)**2
  elseif(mp.eq.2) then
     do l=0,2
        aa(0,1)=(pll(gau(1,1),l-2,l,i,is)&
             +pll(gau(2,1),l-2,l,i,is))*gau(1,2)&
             +(pll(gau(3,1),l-2,l,i,is)&
             +pll(gau(4,1),l-2,l,i,is))*gau(3,2)&
             +(pll(gau(5,1),l-2,l,i,is)&
             +pll(gau(6,1),l-2,l,i,is))*gau(5,2)
        aa(1,1)=(pll(gau(1,1),l-2,l,i,is)*gau(1,1)&
             +pll(gau(2,1),l-2,l,i,is)*gau(2,1))*gau(1,2)&
             +(pll(gau(3,1),l-2,l,i,is)*gau(3,1)&
             +pll(gau(4,1),l-2,l,i,is)*gau(4,1))*gau(3,2)&
             +(pll(gau(5,1),l-2,l,i,is)*gau(5,1)&
             +pll(gau(6,1),l-2,l,i,is)*gau(6,1))*gau(5,2)
        aa(2,1)=(pll(gau(1,1),l-2,l,i,is)*gau(1,1)**2&
             +pll(gau(2,1),l-2,l,i,is)*gau(2,1)**2)*gau(1,2)&
             +(pll(gau(3,1),l-2,l,i,is)*gau(3,1)**2&
             +pll(gau(4,1),l-2,l,i,is)*gau(4,1)**2)*gau(3,2)&
             +(pll(gau(5,1),l-2,l,i,is)*gau(5,1)**2&
             +pll(gau(6,1),l-2,l,i,is)*gau(6,1)**2)*gau(5,2)
        do ll=1,2
           aa(ll,2)=aa(0,1)*aii(ll,0)+aa(1,1)*aii(ll,1)+aa(2,1)*aii(ll,2)
        enddo
        s(l+1,is)=aa(1,2)**2+13.0/3.0*aa(2,2)**2
     enddo
  elseif(mp.eq.3) then
     do l=0,3
        aa(0,1)=(pll(gau(1,1),l-3,l,i,is)&
             +pll(gau(2,1),l-3,l,i,is))*gau(1,2)&
             +(pll(gau(3,1),l-3,l,i,is)&
             +pll(gau(4,1),l-3,l,i,is))*gau(3,2)&
             +(pll(gau(5,1),l-3,l,i,is)&
             +pll(gau(6,1),l-3,l,i,is))*gau(5,2)
        aa(1,1)=(pll(gau(1,1),l-3,l,i,is)*gau(1,1)&
             +pll(gau(2,1),l-3,l,i,is)*gau(2,1))*gau(1,2)&
             +(pll(gau(3,1),l-3,l,i,is)*gau(3,1)&
             +pll(gau(4,1),l-3,l,i,is)*gau(4,1))*gau(3,2)&
             +(pll(gau(5,1),l-3,l,i,is)*gau(5,1)&
             +pll(gau(6,1),l-3,l,i,is)*gau(6,1))*gau(5,2)
        aa(2,1)=(pll(gau(1,1),l-3,l,i,is)*gau(1,1)**2&
             +pll(gau(2,1),l-3,l,i,is)*gau(2,1)**2)*gau(1,2)&
             +(pll(gau(3,1),l-3,l,i,is)*gau(3,1)**2&
             +pll(gau(4,1),l-3,l,i,is)*gau(4,1)**2)*gau(3,2)&
             +(pll(gau(5,1),l-3,l,i,is)*gau(5,1)**2&
             +pll(gau(6,1),l-3,l,i,is)*gau(6,1)**2)*gau(5,2)
        aa(3,1)=(pll(gau(1,1),l-3,l,i,is)*gau(1,1)**3&
             +pll(gau(2,1),l-3,l,i,is)*gau(2,1)**3)*gau(1,2)&
             +(pll(gau(3,1),l-3,l,i,is)*gau(3,1)**3&
             +pll(gau(4,1),l-3,l,i,is)*gau(4,1)**3)*gau(3,2)&
             +(pll(gau(5,1),l-3,l,i,is)*gau(5,1)**3&
             +pll(gau(6,1),l-3,l,i,is)*gau(6,1)**3)*gau(5,2)
        do ll=1,3
           aa(ll,2)=aa(0,1)*aii(ll,0)+aa(1,1)*aii(ll,1)&
                +aa(2,1)*aii(ll,2)+aa(3,1)*aii(ll,3)
        enddo
        s(l+1,is)=3129./80.*aa(3,2)**2+0.5*aa(1,2)*aa(3,2)&
             +13./3.*aa(2,2)**2+aa(1,2)**2
     enddo
  else
     do l=0,4
        aa(0,1)=(pll(gau(1,1),l-4,l,i,is)&
             +pll(gau(2,1),l-4,l,i,is))*gau(1,2)&
             +(pll(gau(3,1),l-4,l,i,is)&
             +pll(gau(4,1),l-4,l,i,is))*gau(3,2)&
             +(pll(gau(5,1),l-4,l,i,is)&
             +pll(gau(6,1),l-4,l,i,is))*gau(5,2)
        aa(1,1)=(pll(gau(1,1),l-4,l,i,is)*gau(1,1)&
             +pll(gau(2,1),l-4,l,i,is)*gau(2,1))*gau(1,2)&
             +(pll(gau(3,1),l-4,l,i,is)*gau(3,1)&
             +pll(gau(4,1),l-4,l,i,is)*gau(4,1))*gau(3,2)&
             +(pll(gau(5,1),l-4,l,i,is)*gau(5,1)&
             +pll(gau(6,1),l-4,l,i,is)*gau(6,1))*gau(5,2)
        aa(2,1)=(pll(gau(1,1),l-4,l,i,is)*gau(1,1)**2&
             +pll(gau(2,1),l-4,l,i,is)*gau(2,1)**2)*gau(1,2)&
             +(pll(gau(3,1),l-4,l,i,is)*gau(3,1)**2&
             +pll(gau(4,1),l-4,l,i,is)*gau(4,1)**2)*gau(3,2)&
             +(pll(gau(5,1),l-4,l,i,is)*gau(5,1)**2&
             +pll(gau(6,1),l-4,l,i,is)*gau(6,1)**2)*gau(5,2)
        aa(3,1)=(pll(gau(1,1),l-4,l,i,is)*gau(1,1)**3&
             +pll(gau(2,1),l-4,l,i,is)*gau(2,1)**3)*gau(1,2)&
             +(pll(gau(3,1),l-4,l,i,is)*gau(3,1)**3&
             +pll(gau(4,1),l-4,l,i,is)*gau(4,1)**3)*gau(3,2)&
             +(pll(gau(5,1),l-4,l,i,is)*gau(5,1)**3&
             +pll(gau(6,1),l-4,l,i,is)*gau(6,1)**3)*gau(5,2)
        aa(4,1)=(pll(gau(1,1),l-4,l,i,is)*gau(1,1)**4&
             +pll(gau(2,1),l-4,l,i,is)*gau(2,1)**4)*gau(1,2)&
             +(pll(gau(3,1),l-4,l,i,is)*gau(3,1)**4&
             +pll(gau(4,1),l-4,l,i,is)*gau(4,1)**4)*gau(3,2)&
             +(pll(gau(5,1),l-4,l,i,is)*gau(5,1)**4&
             +pll(gau(6,1),l-4,l,i,is)*gau(6,1)**4)*gau(5,2)
        do ll=1,4
           aa(ll,2)=aa(0,1)*aii(ll,0)+aa(1,1)*aii(ll,1)&
                +aa(2,1)*aii(ll,2)+aa(3,1)*aii(ll,3)+aa(4,1)*aii(ll,4)
        enddo
        s(l+1,is)=87617./140.*aa(4,2)**2+4.25*aa(2,2)*aa(4,2)&
             +3129./80.*aa(3,2)**2+0.5*aa(1,2)*aa(3,2)&
             +13./3.*aa(2,2)**2+aa(1,2)**2
     enddo
  endif
  return
end subroutine smooth


!!-------------------------------------------
!! reconstruction polynomial on trouble cell by WENO
!!-------------------------------------------
subroutine wenorecon(i,is)
  use data_module
  dimension s(0:10,nv),ww(6),pp1(6),a(-10:10),temp(6)
  weps=1.0e-6
  sr=sqrt(5.0)
  call smooth(i,is,s)
  if(mp.eq.1) then
     do l=1,2
        ko=mp
        x1=gauss(l,1)
        do i0=0,ko
           do j=0,ko
              k1=j-ko
              k2=j
              ll=i0+j-ko
              coef9(i0+1,j+1,l)=qll(x1,k1,k2,ll,i)
           enddo
        enddo
        do kk=ko,0,-1
           ttc=0.0
           do j=ko,kk+1,-1
              ttc=ttc+rco9(j+1,l,1)*coef9(ko+kk-j+1,j+1,l)
           enddo
           rco9(kk+1,l,1)=(qll(x1,-ko,ko,kk,i)-ttc)/coef9(ko+1,kk+1,l)
        enddo
     enddo

     do l=1,2
        temp(l)=0.0
        sum=0.0
        do j=1,ko+1
           sum=rco9(j,l,1)/(weps+s(j,is))**2+sum
        enddo
        do j=1,ko+1
           ww(j)=rco9(j,l,1)/(weps+s(j,is))**2/sum
        enddo

        do ll=1,ko+1
           pp1(ll)=0.0
           do j=1,ko+1
              pp1(ll)=pp1(ll)+coef9(j,ll,l)*u(0,i-ko-2+j+ll,is)
           enddo
        enddo
        do ll=1,ko+1
           temp(l)=temp(l)+ww(ll)*pp1(ll)
        enddo

     enddo
     u(1,i,is)=(temp(2)-temp(1))

  elseif(mp.eq.2) then
     do l=1,4
        ko=mp
        x1=gauss(l,1)
        do i0=0,ko
           do j=0,ko
              k1=j-ko
              k2=j
              ll=i0+j-ko
              coef9(i0+1,j+1,l)=qll(x1,k1,k2,ll,i)
           enddo
        enddo
        do kk=ko,0,-1
           ttc=0.0
           do j=ko,kk+1,-1
              ttc=ttc+rco9(j+1,l,1)*coef9(ko+kk-j+1,j+1,l)
           enddo
           rco9(kk+1,l,1)=(qll(x1,-ko,ko,kk,i)-ttc)/coef9(ko+1,kk+1,l)
        enddo
        !     write(*,*) x1
        !     write(*,*) coef9,rco9
        !     pause
     enddo

     do l=1,4
        temp(l)=0.0
        sum=0.0
        do j=1,ko+1
           sum=rco9(j,l,1)/(weps+s(j,is))**2+sum
        enddo
        do j=1,ko+1
           ww(j)=rco9(j,l,1)/(weps+s(j,is))**2/sum
        enddo

        do ll=1,ko+1
           pp1(ll)=0.0
           do j=1,ko+1
              pp1(ll)=pp1(ll)+coef9(j,ll,l)*u(0,i-ko-2+j+ll,is)
           enddo
        enddo
        do ll=1,ko+1
           temp(l)=temp(l)+ww(ll)*pp1(ll)
        enddo

     enddo

     u(1,i,is)=(temp(2)-temp(1)+sr*(temp(4)-temp(3)))/2.0
     u(2,i,is)=(temp(2)+temp(1)-(temp(4)+temp(3)))*2.5

  elseif(mp.eq.3) then
     do l=1,4
        ko=mp
        x1=gauss(l,1)
        do i0=0,ko
           do j=0,ko
              k1=j-ko
              k2=j
              ll=i0+j-ko
              coef9(i0+1,j+1,l)=qll(x1,k1,k2,ll,i)
           enddo
        enddo
        do kk=ko,0,-1
           ttc=0.0
           do j=ko,kk+1,-1
              ttc=ttc+rco9(j+1,l,1)*coef9(ko+kk-j+1,j+1,l)
           enddo
           rco9(kk+1,l,1)=(qll(x1,-ko,ko,kk,i)-ttc)/coef9(ko+1,kk+1,l)
        enddo
        !     write(*,*) x1
        !     write(*,*) coef9,rco9
        !     pause
     enddo

     do l=1,4
        temp(l)=0.0
        sum=0.0
        do j=1,ko+1
           sum=rco9(j,l,1)/(weps+s(j,is))**2+sum
        enddo
        do j=1,ko+1
           ww(j)=rco9(j,l,1)/(weps+s(j,is))**2/sum
        enddo

        do ll=1,ko+1
           pp1(ll)=0.0
           do j=1,ko+1
              pp1(ll)=pp1(ll)+coef9(j,ll,l)*u(0,i-ko-2+j+ll,is)
           enddo
        enddo
        do ll=1,ko+1
           temp(l)=temp(l)+ww(ll)*pp1(ll)
        enddo

     enddo

     temp(5)=(180.0*u(0,i,is)-9.0*(temp(1)+temp(2))&
          -(temp(3)+temp(4))*49.0)/64.0
     !     reconstruction of u(1,i,is), u(2,i,is),u(3,i,is)
     u(1,i,is)=((fle(1,gauss(1,1))*temp(1)&
          +fle(1,gauss(2,1))*temp(2))*gauss(1,2)&
          +(fle(1,gauss(3,1))*temp(3)&
          +fle(1,gauss(4,1))*temp(4))*gauss(3,2)&
          +(fle(1,gauss(5,1))*temp(5))*gauss(5,2))*12.0
     u(2,i,is)=((fle(2,gauss(1,1))*temp(1)&
          +fle(2,gauss(2,1))*temp(2))*gauss(1,2)&
          +(fle(2,gauss(3,1))*temp(3)&
          +fle(2,gauss(4,1))*temp(4))*gauss(3,2)&
          +(fle(2,gauss(5,1))*temp(5))*gauss(5,2))*180.0
     u(3,i,is)=((fle(3,gauss(1,1))*temp(1)&
          +fle(3,gauss(2,1))*temp(2))*gauss(1,2)&
          +(fle(3,gauss(3,1))*temp(3)&
          +fle(3,gauss(4,1))*temp(4))*gauss(3,2)&
          +(fle(3,gauss(5,1))*temp(5))*gauss(5,2))*2800.0

  else

     do l=1,6
        ko=mp
        x1=gauss(l,1)
        do i0=0,ko
           do j=0,ko
              k1=j-ko
              k2=j
              ll=i0+j-ko
              coef9(i0+1,j+1,l)=qll(x1,k1,k2,ll,i)
           enddo
        enddo
        do kk=ko,0,-1
           ttc=0.0
           do j=ko,kk+1,-1
              ttc=ttc+rco9(j+1,l,1)*coef9(ko+kk-j+1,j+1,l)
           enddo
           rco9(kk+1,l,1)=(qll(x1,-ko,ko,kk,i)-ttc)/coef9(ko+1,kk+1,l)
        enddo
        !     write(*,*) x1
        !     write(*,*) coef9,rco9
        !     pause
     enddo

     do l=1,4
        temp(l)=0.0
        sum=0.0
        do j=1,ko+1
           sum=rco9(j,l,1)/(weps+s(j,is))**2+sum
        enddo
        do j=1,ko+1
           ww(j)=rco9(j,l,1)/(weps+s(j,is))**2/sum
        enddo

        do ll=1,ko+1
           pp1(ll)=0.0
           do j=1,ko+1
              pp1(ll)=pp1(ll)+coef9(j,ll,l)*u(0,i-ko-2+j+ll,is)
           enddo
        enddo
        do ll=1,ko+1
           temp(l)=temp(l)+ww(ll)*pp1(ll)
        enddo
     enddo


     do l=5,6
        sigma(1,l)=0.0
        sigma(2,l)=0.0
        do j=1,ko+1
           rco9(j,l,2)=(3.0*abs(rco9(j,l,1))+rco9(j,l,1))/2.0
           rco9(j,l,1)=(3.0*abs(rco9(j,l,1))-rco9(j,l,1))/2.0
           sigma(1,l)=sigma(1,l)+rco9(j,l,1)
           sigma(2,l)=sigma(2,l)+rco9(j,l,2)
        enddo
        do j=1,ko+1
           rco9(j,l,1)=rco9(j,l,1)/sigma(1,l)
           rco9(j,l,2)=rco9(j,l,2)/sigma(2,l)
        enddo
     enddo
     do l=5,6
        temp(l)=0.0
        sum=0.0
        do j=1,ko+1
           sum=rco9(j,l,1)/(weps+s(j,is))**2+sum
        enddo
        do j=1,ko+1
           ww(j)=rco9(j,l,1)/(weps+s(j,is))**2/sum
        enddo

        do ll=1,ko+1
           pp1(ll)=0.0
           do j=1,ko+1
              pp1(ll)=pp1(ll)+coef9(j,ll,l)*u(0,i-ko-2+j+ll,is)
           enddo
        enddo
        do ll=1,ko+1
           temp(l)=temp(l)+ww(ll)*pp1(ll)
        enddo
        sum=0.0
        do j=1,ko+1
           sum=rco9(j,l,2)/(weps+s(j,is))**2+sum
        enddo
        do j=1,ko+1
           ww(j)=rco9(j,l,2)/(weps+s(j,is))**2/sum
        enddo
        temp(l)=-sigma(1,l)*temp(l)
        do ll=1,ko+1
           temp(l)=temp(l)+ww(ll)*pp1(ll)*sigma(2,l)
        enddo
     enddo
     !     reconstruction of u(1,i,is), u(2,i,is),u(3,i,is),u(4,i,is)
     u(1,i,is)=((fle(1,gauss(1,1))*temp(1)&
          +fle(1,gauss(2,1))*temp(2))*gauss(1,2)&
          +(fle(1,gauss(3,1))*temp(3)&
          +fle(1,gauss(4,1))*temp(4))*gauss(3,2)&
          +(fle(1,gauss(5,1))*temp(5)&
          +fle(1,gauss(6,1))*temp(6))*gauss(5,2))*12.0
     u(2,i,is)=((fle(2,gauss(1,1))*temp(1)&
          +fle(2,gauss(2,1))*temp(2))*gauss(1,2)&
          +(fle(2,gauss(3,1))*temp(3)&
          +fle(2,gauss(4,1))*temp(4))*gauss(3,2)&
          +(fle(2,gauss(5,1))*temp(5)&
          +fle(2,gauss(6,1))*temp(6))*gauss(5,2))*180.0
     u(3,i,is)=((fle(3,gauss(1,1))*temp(1)&
          +fle(3,gauss(2,1))*temp(2))*gauss(1,2)&
          +(fle(3,gauss(3,1))*temp(3)&
          +fle(3,gauss(4,1))*temp(4))*gauss(3,2)&
          +(fle(3,gauss(5,1))*temp(5)&
          +fle(3,gauss(6,1))*temp(6))*gauss(5,2))*2800.0
     u(4,i,is)=((fle(4,gauss(1,1))*temp(1)&
          +fle(4,gauss(2,1))*temp(2))*gauss(1,2)&
          +(fle(4,gauss(3,1))*temp(3)&
          +fle(4,gauss(4,1))*temp(4))*gauss(3,2)&
          +(fle(4,gauss(5,1))*temp(5)&
          +fle(4,gauss(6,1))*temp(6))*gauss(5,2))*44100.0

  endif
  return
end subroutine wenorecon
