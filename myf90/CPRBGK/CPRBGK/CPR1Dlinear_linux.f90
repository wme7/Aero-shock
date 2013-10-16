!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                            Solving 1D Wave Equation using CPR
!                                  (M. L. Yu, ISU&KU)
!Brief introduction:
!
!Equation: linear advection du/dt+a1*du/dx=0
!Boundary: fixed (0 at present)
!
!parameters:
!Ttime  : total time steps
!delt   : time step
!kstage : RK stage, 2 for RK2, 3 for RK3
!np     : order of accuracy
!nc     : No. of cells
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!                            
             module mod_solu
			 implicit none
             
			 !**********grid information *************!
			 type tgrid
                integer::nc                    !nc:no. of cell; nf:no. of face(node for 1D)
				integer::np                    !order of accuracy
				integer::nf                    !nf: no. of faces
             end type tgrid
             
             real,pointer::x(:)             !x(nf): cell interfaces
			 real,pointer::xc(:,:)          !xc(np,nc): coordinates for solution/flux points
			 real::Ja 
			 real::delx
			 real,pointer::Coe_grid(:)                   !Coe_grid(np) for flux reconstruction
			 real,pointer::Coe_RadauL(:),Coe_RadauR(:)
			 real,pointer::Coe_DLag(:,:)               !Coe_DLag(np,np) for flux derivative interpolation on CV(SPs)
			 real,pointer::Coe_Lag(:,:)                !Coe_Lag(np,2) for flux interpolation on interfaces(FPs)
             !**********solutions related******************************!
		     real,pointer::QC(:,:),Res(:,:)                      !QC,Res(np,nc)
			 real,pointer::Flux(:,:)                             !Flux(np,nc)            
             !**********control parameters*****************************!
			 real::PI
			 real,parameter::gama=1.4d0
			 real,parameter::a1=1.0d0
             !**********RK related*******************************!
			 integer,parameter::Ttime=300
			 real,parameter::delt=0.01
			 real,parameter::kstage=3
			 
			 end module mod_solu



			 program CPRoneDLinear
			 use mod_solu
			 implicit none

			 type(tgrid)::grid
			 integer::i,j,itime,ic
			 real::ResSum
			 integer::np,nc,nf       
             
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
			 delx=20.0/nc

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
			    x(i)=-10.0+delx*(i-1)
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
			 !************ sinusoidal waves or exp wave *******************!
			 do ic=1,nc
			    do j=1,np
			       !QC(j,ic)=sin(xc(j,ic))
				   QC(j,ic)=exp(-1.*xc(j,ic)**2)
				end do
			 end do

			  !******************************!
              open(unit=10,file='Residual.dat')
			  do itime=1,Ttime	     
			     if(kstage==2) then
			       call RK2(grid,ResSum)
				 else if(kstage==3) then
				   call RK3(grid,ResSum)
				 end if
				 write(10,*) itime,ResSum
                                 write(*,*) itime,ResSum 
              end do
			  close(10)
301           format(e12.5,1x,e12.5,1x)              

			  open(unit=10,file='solution.dat')
			  do i=1,nc
			     do j=1,np
			        write(10,301) xc(j,i),QC(j,i)
				 end do
			  end do
			  close(10)

			 stop
			 end
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
		 




		 
		                
		        
		    
