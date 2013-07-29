module dg_solvers
  implicit none

contains

  subroutine compute_rhs_dg_ntp(rhs,q,u,v,ksi_x,ksi_y,&
       eta_x,eta_y,jac,psi,psi_ksi,psi_eta,nelem,npts,nqs)
    implicit none
    integer,intent(in)::nqs,npts,nelem
    real(kind=4),dimension(nelem,npts),intent(in)::q,u,v
    real(kind=4),dimension(nelem,nqs),intent(in)::ksi_x,ksi_y,eta_x,eta_y,jac

    real(kind=4),dimension(npts,nqs),intent(in)::psi,psi_ksi,psi_eta
    real(kind=4),dimension(nelem,npts),intent(inout)::rhs
    character(len=36):: fmt_string,fmt_string_nq
    ! local variables and arrays
    integer :: ie,k,i
    
    real(kind=4) :: wq,e_x,e_y,n_x,n_y,h_k,u_k,v_k,q_k,h_e,h_n,dhdx_k,dhdy_k

    ! initialize rhs

    rhs = 0.
    do ie = 1,nelem

       ! loop integration points
       do k = 1,nqs
          wq = jac(ie,k)
          e_x = ksi_x(ie,k)
          e_y = ksi_y(ie,k)
          n_x = eta_x(ie,k)
          n_y = eta_y(ie,k)

          !interpolate at integration points
          u_k = 0.; v_k = 0.; q_k = 0.;
          do i = 1,npts
             h_k = psi(i,k)
             u_k = u_k + h_k*u(ie,i)
             v_k = v_k + h_k*v(ie,i)
             q_k = q_k + h_k*q(ie,i)
          end do !i

          ! interpolate at integration points..more
          do i = 1,npts
             h_e = psi_ksi(i,k)
             h_n = psi_eta(i,k)
             dhdx_k = h_e*e_x+h_n*n_x
             dhdy_k = h_e*e_y+h_n*n_y
             rhs(ie,i) = rhs(ie,i)+wq*q_k*(dhdx_k*u_k+dhdy_k*v_k)
          end do!i
       end do !k
    end do !ie

!!$write(fmt_string,*) '(1X,', npts, 'ES12.3)'
!!$write(fmt_string_nq,*) '(1X,',nqs,'ES12.3)'
!!$
!!$
!!$write (*,*) 'output from compute RHS:'
!!$write(*,fmt_string) ((rhs(i,ie), ie = 1,npts),i = 1,nelem)





  end subroutine compute_rhs_dg_ntp

  subroutine compute_flux_dg_ntp(rhs,q,u,v,psideh,&
       nx,ny,jac_side,psi,nside,ngl,nq,imapl,imapr,nelem,npts)
    implicit none
    integer,intent(in)::nside,ngl,nq,nelem,npts
    real(kind=4),dimension(nelem,npts),intent(in)::q,u,v
    real(kind=4),dimension(nside,nq),intent(in)::nx,ny,jac_side
    real(kind=4),dimension(ngl,nq),intent(in)::psi
    integer,dimension(4,2,ngl),intent(in)::imapl,imapr
    integer,dimension(nside,4),intent(in)::psideh
    real(kind=4),dimension(nelem,npts),intent(inout)::rhs

    ! local arrays and variables
    real(kind=4),dimension(ngl)::ql,qr,ul,vl
    integer::is,iel,ilocl,l,il,jl,kl,ier,ilocr,ir,jr,kr,i
    real(kind=4) :: wq,nxl,nyl,nxr,nyr,qlq_k,qrq_k,u_k,v_k,ul_k,vl_k
    real(kind=4) :: ur_k,vr_k,unl,unr,claml,clamr,clam,fxl,fyl,fxr,fyr
    real(kind=4) :: flux_ql,flux_qr,flux_q,diss_q,h_i
    character(len=36):: fmt_string,fmt_string_nq
    integer :: ie, nqs
    nqs = nq*nq
    ! initialize local arrays
    ql = 0.; qr = 0.; ul = 0.; vl = 0.;

    do is = 1,nside
       !store left side variables
       iel = psideh(is,3)
       if(iel .ne. -6) then
          ilocl = psideh(is,1)
          do l = 1,ngl
             !get pointers
             il = imapl(ilocl,1,l)
             jl = imapl(ilocl,2,l)
             kl = (jl-1)*ngl+il
             !left element
             ql(l)=q(iel,kl)
             ul(l)=u(iel,kl)
             vl(l)=v(iel,kl)
          end do!l

          ! store right side variables
          ier=psideh(is,4)
          if(ier .ne.0) then
             ilocr = psideh(is,2)
             do l = 1,ngl
                ! get pointers
                ir = imapr(ilocr,1,l)
                jr = imapr(ilocr,2,l)
                kr = (jr-1)*ngl + ir
                ! right element
                qr(l)=q(ier,kr)
             end do !l
          endif !if ier...

          ! do gauss-lobatto integration
          do l = 1,nq
             wq = jac_side(is,l)
             !store normal vectors
             nxl = nx(is,l)
             nyl = ny(is,l)
             nxr = -nxl
             nyr = -nyl 

             ! interpolate onto Quadrature points
             qlq_k=0.; qrq_k=0.; u_k = 0.; v_k = 0.;
             do i = 1,ngl
                qlq_k=qlq_k+psi(i,l)*ql(i)
                qrq_k=qrq_k+psi(i,l)*qr(i)
                u_k=u_k+psi(i,l)*ul(i)
                v_k=v_k+psi(i,l)*vl(i)
             end do!i
             ul_k=u_k; vl_k = v_k;
             ur_k = u_k; vr_k = v_k;

             ! computer rusanov flux constant
             unl=nxl*ul_k+nyl*vl_k
             unr=nxl*ur_k+nyl*vr_k
             claml=abs(unl)
             clamr=abs(unr)
             clam = max(claml,clamr)

             ! flux variables
             fxl=qlq_k*ul_k
             fyl=qlq_k*vl_k
             fxr=qrq_k*ur_k
             fyr=qrq_k*vr_k

             !normal flux component
             flux_ql=(nxl*fxl+nyl*fyl)
             flux_qr=(nxr*fxr+nyr*fyr)
             flux_q=flux_ql-flux_qr

             !dissipation term
             diss_q = clam*(qrq_k-qlq_k)

             !construct rusanov flux
             flux_q=0.5*(flux_q-diss_q)

             !loop through side interpolation points
             do i = 1,ngl
                h_i=psi(i,l)

                !left side
                il=imapl(ilocl,1,i)
                jl=imapl(ilocl,2,i)
                kl=(jl-1)*ngl+il
                rhs(iel,kl)=rhs(iel,kl)-wq*h_i*flux_q

                !right side
                if(ier .gt. 0) then
                   ir = imapr(ilocr,1,i)
                   jr = imapr(ilocr,2,i)
                   kr=(jr-1)*ngl+ir
                   rhs(ier,kr)=rhs(ier,kr)+wq*h_i*flux_q
                endif !(ier...
             end do !i
          end do !l
       end if !(iel...

    end do !is

!!$write(fmt_string,*) '(1X,', npts, 'ES12.3)'
!!$write(fmt_string_nq,*) '(1X,',nqs,'ES12.3)'
!!$
!!$
!!$write (*,*) 'output from compute flux:'
!!$write(*,fmt_string) ((rhs(i,ie), ie = 1,npts),i = 1,nelem)

  end subroutine compute_flux_dg_ntp


  subroutine dg_advection_2d_non_tensor_product(ntime,dt,kstages,q0,u0,v0, &
       psi,ksi2d_x,ksi2d_y,eta2d_x,eta2d_y,jac2d,psi2d,psi2d_ksi,psi2d_eta, &
       intma2d,psideh, jac_side, imapl,imapr, nx,ny, nelem,npts,nqs,ngl,nq, &
       Mmatrix_inv,nside)
    implicit none

    ! input arguments
    integer, intent(in)::nside,nq,ngl,nqs,npts,nelem,ntime,kstages
    real(kind=4),intent(in)::dt
    integer,dimension(4,3,ngl),intent(in) :: imapl,imapr
    integer,dimension(nside,4),intent(in) :: psideh
    integer,dimension(nelem,npts),intent(in)::intma2d
    real(kind=4),dimension(ngl,nq),intent(in) :: psi
    real(kind=4),dimension(npts,nqs),intent(in)::psi2d,psi2d_ksi,psi2d_eta
    real(kind=4),dimension(nelem,nqs),intent(in)::ksi2d_x,ksi2d_y,eta2d_x,eta2d_y
    real(kind=4),dimension(nelem,nqs),intent(in)::jac2d
    real(kind=4),dimension(nside,nq),intent(in)::jac_side,nx,ny
    real(kind=4),dimension(nelem,npts,npts),intent(in)::Mmatrix_inv
    ! in/out arguments
    real(kind=4),dimension(nelem,npts),intent(inout)::q0,u0,v0

    ! local variables and arrays
    integer :: itime, ik,e,ie,i
    real(kind=4),dimension(npts,npts) :: Mtemp
    real(kind=4),dimension(npts) :: rtemp
    real(kind=4),dimension(nelem,npts)::rhs
    real(kind=4),dimension(nelem,npts)::qp,q1
    real(kind=4) :: a0,a1,beta,dtt
character(len=36):: fmt_string,fmt_string_nq
integer, parameter :: time_step_rep_freq = 100
    ! initialize the state vector
    rhs = 0.; q1 = q0; qp = q0;
    do itime = 1,ntime

       if(mod(itime,time_step_rep_freq) .eq. 0) then
          write (*,100) itime
100       format(' ','Commencing time step', I6)
       endif

       do ik = 1,kstages

          select case (kstages)
          case (2)
             select case (ik)
             case (1)
                a0 = 1.; a1 = 0.; beta = 1.;
             case (2)
                a0 = 0.5; a1 = 0.5; beta = 0.5;
             end select !ik
          case (3)
             select case (ik)
             case (1)
                a0 = 1.; a1 = 0.; beta = 1.;
             case (2)
                a0 = 3./4.; a1 = 1./4.; beta = 1./4.
             case (3)
                a0 = 1./3.; a1 = 2./3.; beta = 2./3.
             end select !ik
          end select !kstages
          dtt = dt*beta;

          ! compute RHS
          call compute_rhs_dg_ntp(rhs,qp,u0,v0,ksi2d_x,ksi2d_y, &
               eta2d_x,eta2d_y,jac2d,psi2d,psi2d_ksi,psi2d_eta,&
               nelem,npts,nqs)
          ! compute/apply flux
          call compute_flux_dg_ntp(rhs,qp,u0,v0,psideh,nx,ny,&
               jac_side,psi,nside,ngl,nq,imapl,imapr,nelem,npts)
          ! apply inverse mass matrix
          do e = 1,nelem
             Mtemp = Mmatrix_inv(e,:,:)
             rtemp = rhs(e,:)
             rtemp = matmul(Mtemp,rtemp)
             rhs(e,:) = rtemp
          end do !e

!!$write(fmt_string,*) '(1X,', npts, 'ES12.3)'
!!$write(fmt_string_nq,*) '(1X,',nqs,'ES12.3)'
!!$
!!$
!!$write (*,*) 'output after inverting mass matrix:'
!!$write(*,fmt_string) ((rhs(i,ie), ie = 1,npts),i = 1,nelem)


!!$write(*,*) 'dtt = ',dtt
          ! update inter-stage
          qp = a0*q0 + a1*q1+dtt*rhs

!!$write (*,*) 'qp after interstage:'
!!$write(*,fmt_string) ((qp(i,ie), ie = 1,npts),i = 1,nelem)


          ! update
          q1 = qp

!!$write (*,*) 'q1 after interstage:'
!!$write(*,fmt_string) ((q1(i,ie), ie = 1,npts),i = 1,nelem)

       end do !ik

       ! update q0
       q0 = qp

    end do!itime

  end subroutine dg_advection_2d_non_tensor_product

end module dg_solvers
