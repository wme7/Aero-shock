!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                            Solving 1D B-BGK  Equation using CPR
!    (Original code for 1-D linear wave by M. L. Yu, ISU&KU; modified & extended by JY Yang)
!Brief introduction:
!
!Equation: linear advection du/dt+a1*du/dx=0; df^{(0)}/dt + v_{\sigam} df^{(0)}/dx =0
!Boundary: fixed (0 at present)
!
!parameters:
!Ttime  : total time steps
!delt   : time step
!kstage : RK stage, 2 for RK2, 3 for RK3
!np     : order of accuracy
!nc     : No. of cells
!nv     :No. of discrete velocity points
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!                            
             module mod_solu
			 implicit none
             
			 !**********grid information *************!
			 type tgrid
                 integer::nc                    !nc:no. of cell; nf:no. of face(node for 1D)
	   integer::np                    !order of accuracy
	   integer::nf                    !nf: no. of faces
                 integer::nv                    !nv: no. of discrete velocity points
                 integer ighq                   !ighq=0:Gauss-Hermite quadrature; ighq=1:Newton-Cotes
             end type tgrid
             real,pointer::feq(:,:,:)          !feq(np,nc,nv): distribution function
              real,pointer::gh(:),ww(:)
             real,pointer::vv(:),cc(:)
 
             real,pointer::x(:)             !x(nf): cell interfaces
			 real,pointer::xc(:,:)          !xc(np,nc): coordinates for solution/flux points
			 real::Ja 
			 real::delx
                                      real::a1
			 real,pointer::Coe_grid(:)                   !Coe_grid(np) for flux reconstruction
			 real,pointer::Coe_RadauL(:),Coe_RadauR(:)
			 real,pointer::Coe_DLag(:,:)               !Coe_DLag(np,np) for flux derivative interpolation on CV(SPs)
			 real,pointer::Coe_Lag(:,:)                !Coe_Lag(np,2) for flux interpolation on interfaces(FPs)
             !**********solutions related******************************!
		     real,pointer::QC(:,:),Res(:,:)                      !QC,Res(np,nc)
			 real,pointer::Flux(:,:)                             !Flux(np,nc)  
                                real,pointer::rki(:,:),pki(:,:),uki(:,:),tki(:,:)                      !rho, u, p, t(np,nc)
          
             !**********control parameters*****************************!
			 real::PI
			 real,parameter::gama=1.4d0
             !		 real,parameter::a1=1.0d0
                                     
             !**********RK related*******************************!
	!		 integer,parameter::Ttime=300
	!		 real,parameter::delt=0.01
			 real,parameter::kstage=3
			 
			 end module mod_solu

!  Main Program df^{(0)}/dt + v_{\sigam} df^{(0)}/dx =0

program CPRoneDBGK
			 use mod_solu
			 implicit none

			 type(tgrid)::grid
			 integer::i,j,itime,ic,iv
			 real::ResSum
			 integer::np,nc,nf,nv       
             
             PI=atan(1.0)*4.0d0
             !********* specify simulation information *****************! 
             np=4
			 nc=100
			 nf=nc+1
             grid%nc = nc  
		     grid%np = np
			 grid%nf = nf
			 allocate(x(nf),xc(np,nc))
			 allocate(Coe_grid(np),Coe_RadauL(np),Coe_RadauR(np),Coe_DLag(np,np),Coe_Lag(np,2))
			 allocate(QC(np,nc),Res(np,nc))
			 allocate(Flux(np,nc))

			 
			 !*************** geometry representation******************!
   xleft=0.0
   xright=1.0
   xlen = xright - xleft
   delx = xlen/nc
   tprint=0.1
			! delx=20.0/nc

 			 if(np==2) then                     !Guass-Lergendre points plus two ends
 		       Coe_grid(1)=-sqrt(1.0/3.0)
 			   Coe_grid(2)=sqrt(1.0/3.0)
 		     else if(np==3) then
 		       Coe_grid(1)=-sqrt(15.0)/5.0
 			   Coe_grid(2)=0.0
               Coe_grid(3)=sqrt(15.0)/5.0
 		     else if(np==4) then
 		       Coe_grid(1)=-sqrt(525.0+70.0*sqrt(30.))/35.0
 			   Coe_grid(2)=-sqrt(525.0-70.0*sqrt(30.))/35.0
               Coe_grid(3)=sqrt(525.0-70.0*sqrt(30.))/35.0
               Coe_grid(4)=sqrt(525.0+70.0*sqrt(30.))/35.0
 		     else if(np==5) then
 			   Coe_grid(1)=-1.0
 			   Coe_grid(2)=-sqrt(15.0)/5.0
 			   Coe_grid(3)=0.0
               Coe_grid(4)=sqrt(15.0)/5.0
               Coe_grid(5)=1.0
 			 end if

			 if(np==2) then
                           Coe_RadauR(1) = -1.3660254037844386
                           Coe_RadauR(2) = 0.3660254037844386
			 else if(np==3) then
                           Coe_RadauR(1) = -2.6618950038622251
                           Coe_RadauR(2) = 0.75
                           Coe_RadauR(3) = -0.3381049961377749
			 else if(np==4) then
                           Coe_RadauR(1) = -4.3891529665310851
                           Coe_RadauR(2) = 1.2476247709889344
                           Coe_RadauR(3) = -0.6145280959667939
                           Coe_RadauR(4) = 0.3274848629375160
                         else if(np==5) then
                           Coe_RadauR(1) = -6.5480456079367451
                           Coe_RadauR(2) = 1.8660779160350236
                           Coe_RadauR(3) = -0.9375
                           Coe_RadauR(4) = 0.5598111202653962
                           Coe_RadauR(5) = -0.3222878728081190
			 end if 

             do j=1,np
                Coe_RadauL(j) = -Coe_RadauR(np-j+1)
             end do            
             do i=1,nc+1
			 !   x(i)=-10.0+delx*(i-1) ! Cell Interfaces or boundaries
                               x(i)=xleft + delx*(i-1)
			 end do

			 do ic=1,nc            !Quadrature points
			    do i=1,np
				   xc(i,ic)=(Coe_grid(i)+1.0)*(x(ic+1)-x(ic))/2.0+x(ic)
				end do
			 end do
			 
             Ja = delx/2.0			 
 
             do ic=1,nc
			    do j=1,np               !every solution point
				   do i=1,np            !every solution point based DLagrange
				      call DLagrange(coe_DLag(i,j),Coe_grid(j),Coe_grid(1),i,np)
				   end do
				end do
			 end do
                       
			 do ic=1,nc
			    do j=1,np               !every CVs (solution points)
				   call Lagrange(coe_Lag(j,1),-1.0,Coe_grid(1),j,np)
				   call Lagrange(coe_Lag(j,2),1.0,Coe_grid(1),j,np)
				end do
			 end do
			 
				             
			 !***********************initialization*********************!
                  call initial
			 !************ sinusoidal waves or exp wave *******************!
!			 do ic=1,nc
!			    do j=1,np
!			       !QC(j,ic)=sin(xc(j,ic))
!			       !   QC(j,ic)=exp(-1.*xc(j,ic)**2)
!                              	    end do
!			 end do
! Determine delt
                    call setdelt
!
! Do time integration with delt time step
!
 do 
     if((t.ge.tprint-1.e-10).or.(kcount.ge.kcmax)) exit

     if(kcount/1*1.eq.kcount) then
        ! write(6,*) kcount,' t= ',t  
        ! pause
     end if


			  !******************************!
              open(unit=10,file='Residual.dat')
  do itime=1,Ttime
           do iv=1,nv
                 a1= vv(iv)
                   do ic=1,nc
	      do j=1,np
                       QC(j,ic)= feq(j,ic,iv) 
                   end do
                   end do

	          if(kstage==2) then
			       call RK2(grid,ResSum)
				 else if(kstage==3) then
				   call RK3(grid,ResSum)
				 end if
				 write(10,*) itime,ResSum
                                 write(*,*) itime,ResSum 
                      do ic=1,nc
	         do j=1,np
                          feq(j,ic,iv) = QC(j,ic) 
                     end do
                     end do
                enddo
       call macrov
!!!
 print *,kcount,t,dt
     
     t=t+delt
     kcount=kcount+1
     call setdelt
  end do

			  close(10)
301           format(e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x,e12.5,1x)              

			  open(unit=10,file='solution.dat')
			  do i=1,nc
			     do j=1,np
       			        write(10,301) xc(j,i),rki(j,i),pki(j,i),uki(j,i) ,tki(j,i)
			     end do
			  end do
			  close(10)

			 stop
			 end
end program CPRoneDBGK

      !**************************** subroutines ****************************!
	  	   
           subroutine DLagrange(output,input,Coe,ind,np)
		   implicit none

		    integer::np,i,j
			integer::ind
			real::input,Lag,output
			real::Coe(np)
            
			output=0.0
            do j=1,np
			   Lag=1.0
			   do i=1,np
			      if(i/=ind .and. i/=j) then 
			         Lag=(input-Coe(i))/(Coe(ind)-Coe(i))*Lag
			      end if
			   end do
			   if(j/=ind) output=output+Lag/(Coe(ind)-Coe(j))
			end do
            
			end subroutine DLagrange

	  !***********************************************************************!

	      subroutine Lagrange(output,input,Coe,ind,np)
		  implicit none

		    integer::np,i
			integer::ind
			real::input,output
			real::Coe(np)

			output = 1.0d0
			do i=1,np
			   if(i/=ind) then
			     output = output*(input-Coe(i))/(Coe(ind)-Coe(i))
			   end if
            end do

			end subroutine Lagrange
     !*********************************************************************************************************************!
subroutine initial
  use data_module

  !dimension a(0:10), temp(0:n) 
  ! set up the initial condition
  ! u0(x)=.5+sin(pi*x)
  gh(1:5) =(/-5.38748089001,-4.60368244955,-3.94476404012,-3.34785456738,-2.78880605843/)
  gh(6:10) =(/-2.25497400209,-1.73853771212,-1.2340762154,-0.737473728545,-  0.245340708301/)
  gh(11:15) =(/0.245340708301,0.737473728545,1.2340762154,1.73853771212,2.25497400209/)
  gh(16:20) =(/2.78880605843,3.34785456738,3.94476404012,4.60368244955,5.38748089001/)
  ww(1:5)  =(/0.898591961453,0.704332961176,0.62227869619,0.575262442852,0.544851742366/)
  ww(6:10) = (/0.524080350949,0.509679027117,0.499920871336,0.493843385272,0.490921500667/)
  ww(11:15) =(/0.490921500667,0.493843385272,0.499920871336,0.509679027117,0.524080350949/)
  ww(16:20) =(/0.544851742366,0.575262442852,0.62227869619,0.704332961176,0.898591961453/)
  ighq=2
  if (ighq .eq. 1) go to 10
  do iv=1,nv
     cc(iv) =ww(iv)
     vv(iv) = gh(iv)
  end do
  go to 30
10 continue
  v1  = -10
  v2  = 10
  dv  = (v2-v1)/(nv-1.)
  do iv =1,nv
     vv(iv) = v1 + dv * (iv-1.)
  end do
  do iv = 2, (nv-1)
     cc(iv) = 64./45. * dv
     if(mod(iv,4) .eq.1) cc(iv) = 28./45. * dv
     if (mod(iv,4) .eq.3) cc(iv) = 24./45. * dv
  end do
  cc(1)    = 14./45. * dv
  cc(nv)   = cc(1)
30 continue
  pi=4.0*atan(1.0)
  !Initial condition for shock tube-Density: rho; Mean Velocity: uxm; Pressure: pre; Temperature:tem;
  ! u(kk,i,is) will be the unknowns, the degree of freedom

  rl=0.125
  ul=0.
  pl=0.1
  rr=1.0
  ur=0.
  pr=1.

  do ic = 1, nc
     if (x(i) .le. 0.5)then
        do j=1,np
           rki(j,ic) = rl
           uki(j,ic) = ul
           pki(j,ic) = pl
           tki(j,ic) = 2.0d0*pl/rl
        end do
     else
        do j=1,np
           rki(j,ic) = rr
           uki(j,ic) = ur  
           pki(j,ic) = pr
           tki(j,ic) = 2.0d0*pr/rr
        end do
     end if
     do j=1,np       
        do iv = 1, nv
           pres       = (vv(iv)-uki(j,ic))**2/tki(j,ic)
           feq(j,ic,iv)   = rki(j,ic)*exp(-pres)/sqrt(pi*tki(j,ic)) 
        end do
     end do
  end do

!  do iv=1,nv
!     do i=0,n+1
!        do kk=0,mp
!           a(kk)=(u(kk,i,iv)*fle(kk,gau(1,1)) + u(kk,i,iv)*fle(kk,gau(2,1)))*gau(1,2)&
!                +(u(kk,i,iv)*fle(kk,gau(3,1)) + u(kk,i,iv)*fle(kk,gau(4,1)))*gau(3,2)&
!                +(u(kk,i,iv)*fle(kk,gau(5,1)) + u(kk,i,iv)*fle(kk,gau(6,1)))*gau(6,2)
 !       enddo
        ! take care of the mass matrix
 !       do kk=0,mp
!           if(kk.eq.0) then
 !             u(kk,i,iv)=ai(kk,kk)*a(kk)
!           else
!              u(kk,i,iv)=ai(kk,kk)*a(kk)
!           endif
!        enddo
!     enddo
!  enddo
  kcount=0
  t=0.

  return
end subroutine initial
!!-------------------------------------------
!! Find macroscopic properties and equilibrium distribution
!!-------------------------------------------
subroutine macrov
  use data_module
  
  !evaluate the density distribution func at position j
 ! do ic=1,nc
  !     do iv=1,nv
  !        do j=1,np
   !          f(j,ic,iv)=eval(u(0,ic,iv),mp,xi_out(j))
  !        enddo
 !      enddo
 !   enddo

  !compute macro var
  do ic = 1, nc
     do j=1,np
        skr = 0.
        sku = 0.
        ske = 0.
        do iv = 1, nv
           skr = skr + cc(iv) * feq(j,ic,iv)
           sku = sku + cc(iv) * feq(j,ic,iv) * vv(iv)
           ske = ske + cc(iv) * feq(j,ic,iv) * (0.5 * vv(iv) * vv(iv))
        enddo
        rki(j,ic) = skr
        uki(j,ic) = sku/skr  
        eki(j,ic) = ske  
        tki(j,ic) = 4.0d0*eki(j,ic)/rki(j,ic) - 2.0d0*uki(j,ic)*uki(j,ic)
        pki(j,ic) = 0.5d0*rki(j,ic)*tki(j,ic) 

     enddo
  enddo

!  do iv = 1, nv
!     do ic = 1, nc
!        do j=1,np
!           pres  = (vv(iv)-uki(j,ic))**2/tki(j,ic)
!           feq(j,ic,iv) = rki(j,ic)*exp(-pres)/sqrt(pi*tki(j,ic)) 
!       enddo
!     enddo
!  enddo

end subroutine macrov
!!-------------------------------------------
!!
!!-------------------------------------------
subroutine setdelt
  use data_module

  ! set up delt
  cflc=0.5
  dxmin=delx
  amax=0.
  do iv=1,nv
     amax=max(amax,abs(vv(iv)))
  enddo
    tmt =3.0
  if(np.eq.3) then
     if(np.le.2.0+1.0e-5) then
        rr=1.0
     else
        rr=(np+1.0)/tmt
     endif
  else
     if(np.le.3.0+1.0e-5) then
        rr=1.0
     else
        rr=(np+1.0)/tmt
     endif
  endif
  !dt=cflc*dxmin**rr/amax
  delt=cflc*dxmin**rr/amax
  if(t+delt.gt.tprint) delt=tprint-t
  return
end subroutine setdelt

!!-------------------------------------------
!!-------------------------------------------
subroutine tec_output(iter)
  use data_module
  integer iter
  real*8 rho,uuu,mx,me,tt,pp,pos_x

  character(len=30)::it_str
  write(it_str,'(I4)') iter

  !!---------------------------------------
  ! write average sol at the cell center pt
  open(1,file='sol_cell_avg_'//'_p'//char(48+mo)//'_flux'//char(kflux+48)//'_iter'//it_str//'.dat')
  WRITE(1,*) 'VARIABLES = "x","rho","p","t","u" '

  do i=1,n
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

  !!---------------------------------------
  ! write sol at the equally spaced cell pt
  open(2,file='sol_cell_'//'_p'//char(48+mo)//'_flux'//char(kflux+48)//'_iter'//it_str//'.dat')
  WRITE(2,*) 'VARIABLES = "x","rho","p","t","u" '

  do i=1,nc
     do j=1,np
        pos_x=x(i)+xi_out(j)*dx(i)
        write(2,123) pos_x,rki(j,i),pki(j,i),tki(j,i),uki(j,i)
     enddo 
  enddo
  close(2)
  
123 format(f21.16 f21.16 f21.16 f21.16 f21.16)
end subroutine tec_output

 subroutine FluxRec(grid)
			 use mod_solu
			 implicit none
             
			 type(tgrid)::grid
			 integer::np,nc,nf
			 integer::ic,i,j
             real::QCL,QCR,Cs
			 real,pointer::QcSVB(:,:),FluxSVB(:),Flux_com(:)
			 real,pointer::DFlux(:,:)                  !DFlux at SPs
			 real,pointer::Flux_cor_l(:),Flux_cor_r(:),Flux_bond(:,:)
             
			 np=grid%np
			 nc=grid%nc
			 nf=grid%nf

             
			 allocate(QcSVB(2,nf),FluxSVB(nf),Flux_com(nf))
			 allocate(DFlux(np,nc))
			 allocate(Flux_cor_l(nc),Flux_cor_r(nc))
			 allocate(Flux_bond(2,nc))
			 QcSVB=0.0
			 !************* reconstrcut all fluxes from interpolation******************!
			 do ic=1,nc
			    do i=1,np            !every interface(FPs)
				!   Flux(i,ic) = a1 * QC(i,ic) / Ja  (or set al = vsgm(iv))
                                                      Flux(i,ic) = a1 * QC(i,ic) / Ja

				   if(i==1) then                   !store the SV interface valuces
					 QcSVB(1,ic) = 0.d0
					 do j=1,np
					    QcSVB(1,ic) = QcSVB(1,ic) + Coe_Lag(j,1) * QC(j,ic)
					 end do

				   else if (i==np) then
                     QcSVB(2,ic+1) = 0.d0
					 do j=1,np
                        QcSVB(2,ic+1) = QcSVB(2,ic+1) + Coe_Lag(j,2) * QC(j,ic)
					 end do
				   end if
				end do
				Flux_bond(1,ic) = 0.
				Flux_bond(2,ic) = 0.
				do j=1,np
                   Flux_bond(1,ic) = Flux_bond(1,ic) + Coe_Lag(j,1) * Flux(j,ic)
                   Flux_bond(2,ic) = Flux_bond(2,ic) + Coe_Lag(j,2) * Flux(j,ic)
				end do
			  end do
			  
			  
			  !************ use Riemman solver to reconstruct SV boundary fluxes*********!
              do i=1,nf
			     if(i==1) then    !left boundary, periodic
				   QcSVB(2,i)=0.d0
                   do j=1,np
					  QcSVB(2,i) = QcSVB(2,i) + Coe_Lag(j,2) * QC(j,nc)
				   end do
				 else if(i==nf) then   !right boundary
                   QcSVB(1,i)=0.d0
				   do j=1,np
				      QcSVB(1,i) = QcSVB(1,i) + Coe_Lag(j,1) * QC(j,1)
				   end do
				 end if

				   QCL=QcSVB(2,i)

				   QCR=QcSVB(1,i)

				   
                   Cs=abs(a1)

                   FluxSVB(i) = 0.5d0*(QCL*a1+QCR*a1-Cs*(QCR-QCL))
			 end do
			 
			 Flux_com(1:nf) = FluxSVB(1:nf) / Ja
			 
			 do ic=1,nc
			    Flux_cor_l(ic)=Flux_com(ic)-Flux_bond(1,ic)
				Flux_cor_r(ic)=Flux_com(ic+1)-Flux_bond(2,ic)
			 end do
             
			 do ic=1,nc
			    do i=1,np
			       DFlux(i,ic)=0.0
				   do j=1,np
				      DFlux(i,ic)=DFlux(i,ic)+Coe_DLag(j,i)*Flux(j,ic)
				   end do
				end do
			 end do
             
             do ic=1,nc                        !update CVs using correction via flux reconstructions
			    do i=1,np
				   Res(i,ic)=DFlux(i,ic) + Flux_cor_l(ic)*Coe_RadauR(i) + Flux_cor_r(ic)*Coe_RadauL(i)
				end do
				
			 end do

			 Res = -Res
             
			 deallocate(QcSVB,FluxSVB,DFlux,Flux_com,Flux_cor_l,Flux_cor_r,Flux_bond)

			 end subroutine FluxRec

		 !***********************RK2***********************************!
		 subroutine RK2(grid,ResSum)
		 use mod_solu
		 implicit none
		 
		 type(tgrid)::grid
		 integer::ic,i
		 integer::np,nc
		 integer::count
		 real,pointer::Qco(:,:)
		 real::ResSum
		 
		 np=grid%np
		 nc=grid%nc
		 allocate(Qco(np,nc))
		 
		 Qco=QC
		 
		 call FluxRec(grid)
		 do ic=1,nc
		    do i=1,np
		       QC(i,ic)=Qco(i,ic)+delt*Res(i,ic)
			end do
		 end do
		 
		 call FluxRec(grid)
		 do ic=1,nc
		    do i=1,np
		       QC(i,ic)=0.5d0*(Qco(i,ic)+QC(i,ic)+delt*Res(i,ic))
			end do
		 end do
		 
		 ResSum=0.0
		 count=0
         do ic=1,nc
		    do i=1,np
			   count=count+1
		       ResSum=ResSum+abs(QC(i,ic)-Qco(i,ic))
			end do
		 end do
         ResSum=ResSum/count
		 !write(*,*) ResSum

		 deallocate(Qco)
		 end subroutine RK2

		 !***********************RK3***********************************!
		 subroutine RK3(grid,ResSum)
		 use mod_solu
		 implicit none
		 
		 type(tgrid)::grid
		 integer::ic,i
		 integer::np,nc
		 integer::count
		 real,pointer::Qco(:,:)
		 real::ResSum
		 
		 np=grid%np
		 nc=grid%nc
		 allocate(Qco(np,nc))
		 
		 Qco=QC
		 
		 call FluxRec(grid)
		 do ic=1,nc
		    do i=1,np
		       QC(i,ic)=Qco(i,ic)+delt*Res(i,ic)
			end do
		 end do
		 
		 call FluxRec(grid)
		 do ic=1,nc
		    do i=1,np
		       QC(i,ic)=0.75d0*Qco(i,ic)+0.25d0*QC(i,ic)+0.25d0*delt*Res(i,ic)
			end do
		 end do

		 call FluxRec(grid)
		 do ic=1,nc
		    do i=1,np
		       QC(i,ic)=1.0d0/3.0d0*Qco(i,ic)+2.0d0/3.0d0*QC(i,ic)+2.0d0/3.0d0*delt*Res(i,ic)
			end do
		 end do
		 
		 ResSum=0.0
		 count=0
         do ic=1,nc
		    do i=1,np
			   count=count+1
		       ResSum=ResSum+abs(QC(i,ic)-Qco(i,ic))
			end do
		 end do
         ResSum=ResSum/count
		 !write(*,*) ResSum

		 deallocate(Qco)
		 end subroutine RK3
		 




		 
		                
		        
		    
