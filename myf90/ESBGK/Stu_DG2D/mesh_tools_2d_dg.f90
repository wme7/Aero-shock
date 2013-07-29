module meshtools_2d_dg

  implicit none

contains

  subroutine create_grid_2d_dg(coord,intma,bsido, &
       iperiodic,xgl, npoin, nelem, nboun, &
       nelx, nely, ngl)
    implicit none
    integer, intent(in) :: npoin, nelem, nboun
    integer, intent(in) :: nelx, nely, ngl
    real(kind=4),dimension(ngl),intent(in) :: xgl
    real(kind=4),dimension(npoin,2), intent(inout) :: coord
    integer, dimension(nelem,ngl,ngl),intent(inout)::intma
    integer, dimension(nboun,4), intent(inout)::bsido
    integer, dimension(npoin), intent(inout)::iperiodic

    integer, dimension(npoin,npoin)::node

    real(kind=4) :: xmin, xmax, ymin, ymax, dx, dy
    real(kind=4) :: x,y,x0,y0
    integer :: nop, nx, ny
    integer :: ip,jj,j,k,l,l1,l2,ii,i,j1,ie,ib,i1,i2,ip1,ip2

!write(*,*) 'entering the subroutine'
    node = 0

    xmin = -1.; xmax = 1.; ymin = -1.; ymax = 1.;
    dx = (xmax-xmin)/real(nelx)
    dy = (ymax-ymin)/real(nely)
    nop = ngl-1;
    nx = nelx*nop+1
    ny = nely*nop+1

!write(*,*) 'generate coord'
    ! generate coord
    ip = 0; jj = 0;
    do k = 1,nely
       y0=ymin+real(k-1)*dy
       if(k .eq. 1) then
          l1 = 1
       else
          l1 = 2
       endif !if(k .eq...

       do l = l1,ngl
          y = (xgl(l)+1.)*dy/2. + y0
          jj=jj+1
          ii=0
          do i=1,nelx
             x0 = xmin + real(i-1)*dx
             if (i .eq. 1) then
                j1 = 1
             else 
                j1 = 2
             endif !if(i .eq. ...

             do j = j1,ngl
                ii = ii+1
                ip = ip+1
                x = (xgl(j)+1.)*dx/2. + x0
                coord(ip,1) = x
                coord(ip,2) = y
                node(ii,jj) = ip
             end do !j
          end do !i
       end do !l
    end do !k

!write(*,*) 'generating intma'

    ! generate intma
    ie = 0
    do k = 1,nely
       do i = 1,nelx
          ie = ie+1
          do l = 1,ngl
             jj = (ngl-1)*(k-1)+l
             do j = 1,ngl
                ii = (ngl-1)*(i-1)+j
                ip = node(ii,jj)
                intma(ie,j,l) = ip
             end do !j
          end do !l
       end do !i
    end do !k

!write(*,*) 'generating bsido'
    ! generate bsido
    ib = 0
    do i = 1,nelx
       ie = i
       ib = ib+1
       i1 = (i-1)*(ngl-1)+1
       i2 = (i-1)*(ngl-1) + ngl
       ip1 = node(i1,1)
       ip2 = node(i2,1)
       bsido(ib,1)=ip1
       bsido(ib,2)=ip2
       bsido(ib,3)=ie
       bsido(ib,4)=6
    end do

    ! right boundary
    do i = 1,nely
       ie = nelx*i
       ib = ib+1
       i1=(i-1)*(ngl-1)+1
       i2=(i-1)*(ngl-1)+ngl
       ip1=node(nx,i1)
       ip2=node(nx,i2)
       bsido(ib,1)=ip1
       bsido(ib,2)=ip2
       bsido(ib,3)=ie
       bsido(ib,4)=6
    end do

    ! top boundary
    do i = nelx,1,-1
       ie=nelem-(nelx-i)
       ib=ib+1
       i1=(i-1)*(ngl-1)+ngl
       i2=(i-1)*(ngl-1)+1
       ip1=node(i1,ny)
       ip2=node(i2,ny)
       bsido(ib,1)=ip1
       bsido(ib,2)=ip2
       bsido(ib,3)=ie
       bsido(ib,4)=6
    end do

    ! left boundary
    do i = nely,1,-1
       ie=nelx*(i-1)+1
       ib=ib+1
       i1=(i-1)*(ngl-1)+ngl
       i2=(i-1)*(ngl-1)+1
       ip1=node(1,i1)
       ip2=node(1,i2)
       bsido(ib,1)=ip1
       bsido(ib,2)=ip2
       bsido(ib,3)=ie
       bsido(ib,4)=6
    end do
!write(*,*) 'periodicity'

    ! periodicity
    do i = 1,npoin
       iperiodic(i) = i
    end do

!write(*,*) 'x-periodicity'
    ! x-periodicity
    do i = 1,ny
       i1 = node(1,i)
       i2 = node(nx,i)
       iperiodic(i2) = i1
    end do
!write(*,*) 'y-periodicity'
    ! y-periodicity
    do i = 1,nx
       i1 = node(i,1)
       i2 = node(i,ny)
       iperiodic(i2)=i1
    end do

  end subroutine create_grid_2d_dg

  subroutine map_deriv_2d(f_ksi,f_eta,psi,dpsi,f,ngl,nq)
    implicit none
    integer, intent(in) :: ngl,nq
    real(kind=4), dimension(ngl,ngl),intent(in) :: f
    real(kind=4), dimension(ngl,nq),intent(in)::psi,dpsi
    real(kind=4), dimension(nq,nq),intent(inout) :: f_ksi,f_eta

    integer :: l,k,i,j
    real(kind=4) :: sum_ksi, sum_eta

    f_ksi = 0.; f_eta = 0.

    do l = 1,nq
       do k = 1,nq
          sum_ksi = 0.
          sum_eta = 0.

          do j = 1,ngl
             do i = 1,ngl
                sum_ksi = sum_ksi + dpsi(i,k)*psi(j,l)*f(i,j)
                sum_eta = sum_eta + psi(i,k)*dpsi(j,l)*f(i,j)
             end do ! i
          end do ! j
          f_ksi(k,l) = sum_ksi
          f_eta(k,l) = sum_eta

       end do ! k
    end do ! l



  end subroutine map_deriv_2d


  subroutine metrics_2d(ksi_x,ksi_y,eta_x,eta_y,jac, &
       coord,intma,psi,dpsi, &
       wnq,npoin,nelem,ngl,nq)

    implicit none
    integer, intent(in) :: ngl, nq, nelem, npoin
    real(kind=4),dimension(nq), intent(in) :: wnq
    integer, dimension(nelem,ngl,ngl), intent(in) :: intma
    real(kind=4),dimension(ngl,nq),intent(in) :: psi, dpsi
    real(kind=4),dimension(npoin,2),intent(in) :: coord
    real(kind=4),dimension(nelem,nq,nq),intent(inout):: ksi_x,ksi_y,eta_x,eta_y,jac

    ! local helper arrays and index variables
    real(kind=4),dimension(ngl,ngl) :: x, y
    real(kind=4),dimension(nq,nq) :: x_ksi, x_eta, y_ksi,y_eta
    integer :: ie, j, i, ip
    real(kind=4) :: xjac

    ! initialize the variables
    x_ksi = 0.; x_eta = 0.;
    y_ksi = 0.; y_eta = 0.;
    ksi_x = 0.; ksi_y = 0.;
    eta_x = 0.; eta_y = 0.;
    jac = 0.;

    x = 0.
    y = 0.



    ! loop thru the elements
    do ie = 1,nelem

       ! store element variables
       do j = 1,ngl
          do i = 1,ngl
             ip = intma(ie,i,j)
             x(i,j) = coord(ip,1)
             y(i,j) = coord(ip,2)
          end do !i
       end do !j


       ! construct mapping derivatives
       call map_deriv_2d(x_ksi,x_eta,psi,dpsi,x,ngl,nq)
       call map_deriv_2d(y_ksi,y_eta,psi,dpsi,y,ngl,nq)


       ! construct inverse mapping
       do j = 1,nq
          do i = 1,nq
             xjac = x_ksi(i,j)*y_eta(i,j) - x_eta(i,j)*y_ksi(i,j)
             ksi_x(ie,i,j)=1./xjac*y_eta(i,j)
             ksi_y(ie,i,j)=-1./xjac*x_eta(i,j)
             eta_x(ie,i,j)=-1./xjac*y_ksi(i,j)
             eta_y(ie,i,j)=1./xjac*x_ksi(i,j)
             jac(ie,i,j)=wnq(i)*wnq(j)*abs(xjac)
          end do !i
       end do !j
    end do !ie

  end subroutine metrics_2d

  subroutine create_side_2d(iside,jeside,intma,bsido, &
       npoin,nelem,nboun,nside,ngl)
    implicit none
    ! arguments
    integer, intent(in) :: npoin,nelem,nboun,nside,ngl
    integer,dimension(nboun,4),intent(in) :: bsido
    integer,dimension(nelem,ngl,ngl),intent(in)::intma
    integer,dimension(nside,4),intent(inout) ::iside
    integer,dimension(nelem,4),intent(inout) :: jeside

    ! local arrays
    integer, dimension(npoin) :: lwher,lhowm
    integer,dimension(5*npoin) :: icone
    integer,dimension(4)::inode,jnode
    ! a whole butt-ton of local integer variables and counters.
    integer ::in,ie,ip,iloca,iloc1,iele,iwher,ip1,in1,ipt,j
    integer :: iold,in2,ip2,is,jloca,il,ir,is1,is2,i1,i2,iel,jnod 
    integer :: ib,ibe,ibc,ilb,irb,ier

    ! initialize some local arrays
    lwher = 0
    lhowm = 0
    icone = 0
    inode = 0
    jnode = 0

    ! fix lnode
    inode(1) = 1
    inode(2) = ngl
    inode(3) = ngl
    inode(4) = 1
    jnode(1) = 1
    jnode(2) = 1
    jnode(3) = ngl
    jnode(4) = ngl



    ! count how many elements own each node
    do in =1,4
       do ie = 1,nelem
          ip = intma(ie,inode(in),jnode(in))
          lhowm(ip)=lhowm(ip)+1
       end do
    end do


    ! track elements owning each node
    lwher(1) = 0
    do ip =2,npoin
       lwher(ip)=lwher(ip-1)+lhowm(ip-1)
    end do



    !another tracker array
    lhowm = 0

    do in = 1,4
       do ie = 1,nelem

          ip = intma(ie,inode(in),jnode(in))

          lhowm(ip)=lhowm(ip)+1

          jloca=lwher(ip)+lhowm(ip)

          icone(jloca)=ie
       end do
    end do


    ! loop over the nodes
    iloca = 0
    do ip = 1,npoin
       iloc1 = iloca
       iele=lhowm(ip)
       if(iele .ne. 0) then
          iwher=lwher(ip)


          ! loop over those elements surrounding node ip
          ip1 = ip
          do iel = 1,iele
             ie = icone(iwher+iel)
             ! find out position of ip in intma
             do in = 1,4
                in1 = in
                ipt = intma(ie,inode(in),jnode(in))
                if(ipt .eq. ip) then
                   exit
                endif
             end do !in
             !check edge element ie which claims ip
             j = 0
             do jnod = 1,3,2
                iold = 0
                j=j+1
                in2 =in+jnod
                if(in2 .gt. 4) then
                   in2 = in2-4
                endif
                ip2 = intma(ie,inode(in2),jnode(in2))
                if(ip2 .ge. ip1) then
                   ! check whether side is old or new
                   if(iloca .ne. iloc1) then
                      do is = (iloc1+1),iloca
                         ! iside(is,2)
                         jloca=is
                         if(iside(is,2).eq.ip2) then
                            iold=1
                            exit
                         endif
                      enddo !is
                   endif !iloca

                   if(iold .eq. 0) then
                      iloca=iloca+1
                      iside(iloca,1)=ip1
                      iside(iloca,2)=ip2
                      iside(iloca,2+j)=ie
                   elseif (iold .eq. 1) then
                      iside(jloca,2+j)=ie
                   endif!iold
                endif!ip2
             enddo!jnod
          enddo !iel

          !perform some shifting to order the nodes of a side in CCW
          do is = (iloc1+1),iloca
             if(iside(is,3).eq. 0) then
                iside(is,3)=iside(is,4)
                iside(is,4) = 0
                iside(is,1) = iside(is,2)
                iside(is,2) = ip1
             endif!iside
          enddo!is
       endif !if iele
    enddo !ip

    if(iloca .ne. nside) then
       write (*,*) 'Error in SIDE. iloca nside = ', iloca, nside
    endif

    !reset the boundary markers
    do is = 1,nside
       if(iside(is,4) .eq. 0) then
          il=iside(is,1)
          ir = iside(is,2)
          ie = iside(is,3)
          do ib = 1,nboun
             ibe = bsido(ib,3)
             ibc=bsido(ib,4)
             if(ibe .eq. ie) then
                ilb = bsido(ib,1)
                irb = bsido(ib,2)
                if((ilb .eq. il) .and. (irb .eq. ir)) then
                   iside(is,4) = -ibc
                   exit
                endif !ilb
             endif !ibe
          enddo !ib
       endif !iside
    enddo ! is


    ! form element/side connectivity array
    do is = 1,nside
       iel = iside(is,3)
       ier = iside(is,4)
       is1=iside(is,1)
       is2=iside(is,2)
       !left side
       do in = 1,4
          i1 = intma(iel,inode(in),jnode(in))
          in1 = in+1
          if(in1 .gt. 4) then
             in1 = 1
          endif
          i2 = intma(iel,inode(in1),jnode(in1)) 
          if ((is1 .eq. i1) .and. (is2 .eq. i2)) then
             jeside(iel,in) = is
          endif !is1
       enddo !in

       ! right side
       if(ier .gt. 0) then
          do in = 1,4
             i1 = intma(ier,inode(in),jnode(in))
             in1 = in+1
             if(in1 .gt. 4) then
                in1 = 1
             endif
             i2 = intma(ier,inode(in1),jnode(in1))
             if((is1 .eq. i2) .and. (is2 .eq. i1)) then
                jeside(ier,in) = is
             endif !is1
          enddo !in
       endif !ier
    enddo !is

  end subroutine create_side_2d

  subroutine create_side_2d_dg(psideh,imapl,imapr,iside, &
       intma,nside,nelem,ngl)
    implicit none
    integer, intent(in) :: ngl,nelem,nside
    integer, dimension(nelem,ngl,ngl),intent(in)::intma
    integer, dimension(nside,4),intent(in) :: iside
    integer, dimension(4,2,ngl),intent(inout) :: imapr, imapl
    integer, dimension(nside,4), intent(inout) :: psideh

    ! local arrays
    integer, dimension(4):: inode, jnode

    ! local index variables and other temps
    integer :: l,i,ip1,ip2,iel,ier,j,j1,j2,jp1,jp2

    psideh = 0
    imapl = 0
    imapr = 0

    inode = 0
    jnode = 0

    ! construct boundary pointer
    inode(1) = 1
    inode(2) = ngl
    inode(3) = ngl
    inode(4) = 1
    jnode(1) = 1
    jnode(2) = 1
    jnode(3) = ngl
    jnode(4) = ngl

    ! construct IMAP arrays
    do l = 1,ngl
       imapl(1,1,l) = l
       imapl(1,2,l) = 1
       imapr(1,1,l) = ngl+1-l
       imapr(1,2,l) = 1

       imapl(2,1,l) = ngl
       imapl(2,2,l) = l
       imapr(2,1,l)=ngl
       imapr(2,2,l)=ngl+1-l

       imapl(3,1,l) = ngl+1-l
       imapl(3,2,l)=ngl
       imapr(3,1,l)=l
       imapr(3,2,l)=ngl

       imapl(4,1,l)=1
       imapl(4,2,l)=ngl+1-l
       imapr(4,1,l)=1
       imapr(4,2,l)=l
    end do

    ! loop through the sides
    do i = 1,nside
       ip1=iside(i,1)
       ip2=iside(i,2)
       iel=iside(i,3)
       ier=iside(i,4)

       ! check for position on left element
       do j=1,4
          j1=j
          j2=j+1
          if(j2 .gt. 4) then
             j2 = 1
          end if
          jp1=intma(iel,inode(j1),jnode(j1))
          jp2=intma(iel,inode(j2),jnode(j2))

          if((ip1 .eq. jp1) .and. (ip2 .eq. jp2)) then
             psideh(i,1)=j
             exit
          end if
       end do

       !check for position on Right Element
       if(ier .gt. 0) then
          do j = 1,4
             j1 = j
             j2 = j+1
             if(j2 .gt. 4) then
                j2 = 1
             endif
             jp1 = intma(ier,inode(j1),jnode(j1))
             jp2 = intma(ier,inode(j2),jnode(j2))

             if((ip1 .eq. jp2) .and. (ip2 .eq. jp1)) then
                psideh(i,2)=j

                exit
             endif
          end do
       endif

       ! store elements into PSIDEH
       psideh(i,3) = iel
       psideh(i,4) = ier

    end do! i


  end subroutine create_side_2d_dg

  subroutine compute_normals_2d(nx,ny,jac_side,psideh,&
       intma,coord,nside,ngl,nq,wnq,psi,dpsi,nelem,npoin)
    implicit none

    integer, intent(in) :: nq,ngl,nside,nelem,npoin
    real(kind=4),dimension(ngl,nq),intent(in) :: psi,dpsi
    real(kind=4),dimension(nq),intent(in):: wnq
    real(kind=4),dimension(npoin,2),intent(in)::coord
    integer,dimension(nelem,ngl,ngl),intent(in)::intma
    integer,dimension(nside,4),intent(in)::psideh
    real(kind=4),dimension(nside,nq)::nx,ny,jac_side

    ! local arrays
    real(kind=4),dimension(ngl,ngl) :: x,y
    real(kind=4),dimension(nq,nq) :: x_ksi,x_eta,y_ksi,y_eta

    ! local aux variables
    integer :: is, ilocl,ilocr,iel,ier,i,j,ip,l
    real(kind=4) :: wq, rnx

    ! initialize variables
    nx = 0.; ny = 0.; jac_side = 0.;
    x = 0.; y = 0.; x_ksi = 0.; x_eta = 0.;
    y_ksi = 0.; y_eta = 0.


    ! loop thru the sides
    do is = 1,nside

       ! store left and right elements
       ilocl=psideh(is,1)
       ilocr=psideh(is,2)
       iel=psideh(is,3)
       ier=psideh(is,4)

       ! store element variables
       do j = 1,ngl
          do i = 1,ngl
             ip = intma(iel,i,j)
             x(i,j) = coord(ip,1)
             y(i,j) = coord(ip,2)
          end do
       end do

       ! construct mapping derivatives
       call map_deriv_2d(x_ksi,x_eta,psi,dpsi,x,ngl,nq)
       call map_deriv_2d(y_ksi,y_eta,psi,dpsi,y,ngl,nq)

       ! compute normals
       do l = 1,nq
          wq = wnq(l)
          select case (ilocl)
          case (1)
             i = l
             j = 1
             nx(is,l)=y_ksi(i,j)
             ny(is,l)=-x_ksi(i,j)
             jac_side(is,l)=wq*sqrt(x_ksi(i,j)**2 + y_ksi(i,j)**2)

          case (2)
             i = nq
             j = l
             nx(is,l)=y_eta(i,j)
             ny(is,l)=-x_eta(i,j)
             jac_side(is,l)=wq*sqrt(x_eta(i,j)**2+y_eta(i,j)**2)

          case (3)
             i=nq+1-l
             j=nq
             nx(is,l)=-y_ksi(i,j)
             ny(is,l)=x_ksi(i,j)
             jac_side(is,l)=wq*sqrt(x_ksi(i,j)**2+y_ksi(i,j)**2)

          case (4)
             i=1
             j=nq+1-l
             nx(is,l)=-y_eta(i,j)
             ny(is,l)=x_eta(i,j)
             jac_side(is,l)=wq*sqrt(x_eta(i,j)**2+y_eta(i,j)**2)


          end select
       end do !l

       ! normalize normals
       do l = 1,nq
          rnx=sqrt(nx(is,l)**2 + ny(is,l)**2)
          nx(is,l)=nx(is,l)/rnx
          ny(is,l)=ny(is,l)/rnx
       end do
    end do ! is

  end subroutine compute_normals_2d

  subroutine create_periodicity_2d(psideh,iside,coord,nside,nboun,npoin)
    implicit none
    integer, intent(in) :: nside,nboun,npoin
    real(kind=4),dimension(npoin,2),intent(in)::coord
    integer,dimension(nside,4),intent(inout)::iside
    integer,dimension(nside,4),intent(inout)::psideh

    ! local arrays
    real(kind=4),parameter :: tol = 1e-6

    integer, dimension(nboun/4) :: ileft,iright,itop,ibot
    ! note the hack...I believe this will only really work on domains where nelx = nely

    integer :: nleft,nright,ntop,nbot

    real(kind=4) :: xmax,xmin,ymax,ymin,xm,ym

    integer :: j,is,ier,i1,i2,i,isr,isl
real(kind=4) :: yl1,yr2,xl1,xr2

    ! initialize
    nleft = 0; nright = 0; ntop = 0; nbot = 0;
    ileft = 0; iright = 0; itop = 0; ibot = 0;

    ! find extrema of domain
    xmax = maxval(coord(:,1)); xmin = minval(coord(:,1));
    ymax = maxval(coord(:,2)); ymin = minval(coord(:,2));

    !loop thru sides and extract Left, Right, Bot and top
    do is = 1,nside
       ier = psideh(is,4)
       if(ier .eq. -6) then
          i1 = iside(is,1); i2 = iside(is,2)
          xm = 0.5*(coord(i1,1)+coord(i2,1))
          ym = 0.5*(coord(i1,2)+coord(i2,2))

          ! check grid point
          if(abs(xm-xmin) .lt. tol) then
             nleft = nleft + 1
             ileft(nleft) = is
          elseif (abs(xm-xmax) .lt. tol) then
             nright = nright + 1
             iright(nright)=is
          elseif(abs(ym-ymin) .lt. tol) then
             nbot=nbot+1
             ibot(nbot)=is
          elseif (abs(ym-ymax) .lt. tol) then
             ntop = ntop + 1
             itop(ntop) = is
          else
             write (*,*) 'No match in PERIODIC_BCS for is ier = '
             write (*,*) is, ier
          endif
       endif !ier
    enddo

!!$write(*,*) 'nbot = ', nbot
!!$write(*,*) 'ntop = ', ntop
!!$write(*,*) 'nleft = ', nleft
!!$write(*,*) 'nright = ', nright
!!$
!!$write (*,*) 'itop = ', itop
!!$write (*,*) 'ibot = ', ibot
!!$write (*,*) 'ileft = ', ileft
!!$write (*,*) 'iright = ', iright

    ! loop through periodic BCs
    ! first - do left and right
    
    do i = 1,nleft
       isl=ileft(i)
       i1=iside(isl,1)
       yl1=coord(i1,2)
       !search for corresponding right edge
       do j=1,nright
          isr=iright(j)
          i2=iside(isr,2)
          yr2=coord(i2,2)
          if(abs(yl1-yr2) .lt. tol) then
!!$write(*,*) 'found edge. isl=',isl,'isr=',isr
             psideh(isl,2)=psideh(isr,1)
             psideh(isl,4)=psideh(isr,3)
             psideh(isr,3)=-6
             iside(isl,4)=iside(isr,3)
             iside(isr,3)=-6
             exit
          endif !abs(yl2...
       end do ! j
    end do !i

    ! second, do top and bottom
    do i = 1,nbot
       isl=ibot(i)
       i1=iside(isl,1)
       xl1=coord(i1,1)

       ! search for corresponding top edge
       do j=1,ntop
          isr = itop(j)
          i2=iside(isr,2)
          xr2=coord(i2,1)
          if(abs(xl1-xr2) .lt. tol) then
!!$write(*,*) 'found edge. isl=',isl,'isr=',isr
             psideh(isl,2)=psideh(isr,1)
             psideh(isl,4)=psideh(isr,3)
             psideh(isr,3)=-6
             iside(isl,4)=iside(isr,3)
             iside(isr,3)=-6
             exit
          endif
       enddo
    enddo


  end subroutine create_periodicity_2d

  subroutine reshape_1d_to_2d(psi2d,psi2d_ksi,psi2d_eta,&
       ksi2d_x,ksi2d_y,eta2d_x,eta2d_y, &
       jac2d,intma2d,npts,nqs,psi,dpsi,&
       ksi_x,ksi_y,eta_x,eta_y,jac,intma,&
       nelem,ngl,nq)
    integer,intent(in) :: nq,ngl,nelem,npts,nqs
    real(kind=4),dimension(nelem,nq,nq),intent(in)::ksi_x,ksi_y,eta_x,eta_y,jac
    real(kind=4),dimension(ngl,nq),intent(in)::psi,dpsi
    integer,dimension(nelem,ngl,ngl),intent(in)::intma
    real(kind=4),dimension(npts,nqs),intent(inout)::psi2d,psi2d_ksi,psi2d_eta
    real(kind=4),dimension(nelem,nqs),intent(inout)::ksi2d_x,ksi2d_y
    real(kind=4),dimension(nelem,nqs),intent(inout)::eta2d_x,eta2d_y
    real(kind=4),dimension(nelem,nqs),intent(inout)::jac2d
    integer,dimension(nelem,npts),intent(inout)::intma2d

    ! local variables
    integer :: e,k,j,i,n,m,l

    ! initialize arrays
    psi2d=0.; psi2d_ksi =0.; psi2d_eta =0.

    do k=1,ngl
       do j = 1,ngl
          i=(k-1)*ngl+j
          do n = 1,nq
             do m = 1,nq
                l=(n-1)*nq + m
                psi2d(i,l)=psi(j,m)*psi(k,n)
                psi2d_ksi(i,l)=dpsi(j,m)*psi(k,n)
                psi2d_eta(i,l)=psi(j,m)*dpsi(k,n)
             end do !m
          end do !n
       end do!j
    end do !k

    !initialize the next set of arrays
    ksi2d_x = 0.; ksi2d_y = 0.;
    eta2d_x = 0.; eta2d_y = 0.;
    jac2d = 0.
    do e = 1,nelem
       do k = 1,nq
          do j = 1,nq
             i=(k-1)*nq+j
             ksi2d_x(e,i)=ksi_x(e,j,k)
             ksi2d_y(e,i)=ksi_y(e,j,k)
             eta2d_x(e,i)=eta_x(e,j,k)
             eta2d_y(e,i)=eta_y(e,j,k)
             jac2d(e,i)=jac(e,j,k)
          end do !j
       end do !k
    end do !e

    ! intma
    intma2d = 0
    do e=1,nelem
       do k = 1,ngl
          do j = 1,ngl
             i=(k-1)*ngl+j
             intma2d(e,i)=intma(e,j,k)
          end do !j
       end do !k
    end do !e

  end subroutine reshape_1d_to_2d

  subroutine create_Mmatrix2d_dg_NonTensorProduct(Mmatrix, &
       jac2d,psi2d,nelem,npts,nqs)
    integer, intent(in) :: nqs,npts,nelem
    real(kind=4),dimension(npts,nqs),intent(in) :: psi2d
    real(kind=4),dimension(nelem,nqs),intent(in)::jac2d
    real(kind=4),dimension(nelem,npts,npts),intent(inout) :: Mmatrix

    ! local variables
    integer :: ie,k,i,j
    real(kind=4) :: wq,h_i,h_j

    ! initialize Mmatrix
    Mmatrix = 0.

    do ie = 1,nelem
       do k = 1,nqs
          wq = jac2d(ie,k)
          do i = 1,npts
             h_i = psi2d(i,k)
             do j = 1,npts
                h_j = psi2d(j,k)
                Mmatrix(ie,i,j)=Mmatrix(ie,i,j)+wq*h_i*h_j
             end do ! j
          end do ! i 
       end do ! k
    end do ! ie


  end subroutine create_Mmatrix2d_dg_NonTensorProduct

end module meshtools_2d_dg
