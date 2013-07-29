module problem_setup_2d_dg
  implicit none

contains


  subroutine exact_solution(qe,ue,ve,coord,npoin,&
       time,icase)
    implicit none
    integer, intent(in) :: npoin, icase
    real(kind=4),intent(in) :: time
    real(kind=4),dimension(npoin,2),intent(in)::coord
    real(kind=4),dimension(npoin),intent(inout)::qe,ue,ve

    ! local constants
    real(kind=4)::w,visc,sigma0,sigma,rc,a,den,timec,x,y,r,xx,yy
    real(kind=4)::xmin,xmax,ymin,ymax,xl,yl,xm,ym,xc,yc
    integer :: ip

    ! initialize
    qe = 0.; ue = 0.; ve = 0.;

    !set local constants
    w = 1.0; visc = 0.0;
    xmin = minval(coord(:,1));
    xmax = maxval(coord(:,1));
    ymin = minval(coord(:,2));
    ymax = maxval(coord(:,2));
    xl = xmax - xmin;
    yl = ymax - ymin;
    xm = 0.5*(xmax+xmin)
    ym = 0.5*(ymax+ymin)
    xc = xmin+0.25*xl
    yc = ymin+0.5*yl;
    sigma0 = 0.125
    rc =0.25
    sigma=sqrt(sigma0**2+2*visc*time)
    a = 1./(1.+2.*visc*time/(sigma0**2))
    den = 2.*(sigma0**2+2*visc*time)
    timec=time-floor(time)

    ! generate grid points
    do ip=1,npoin
       x = coord(ip,1)
       y = coord(ip,2)
       r = sqrt((x-xc)**2 + (y-yc)**2)
       xx = x-xc*cos(time)-yc*sin(time)
       yy = y+xc*sin(time)-yc*cos(time)
       qe(ip) = a*exp(-(xx**2 + yy**2)/den)
       select case (icase)

       case (1) ! gaussian in CCW direction
          ue(ip)=w*(y-ym)
          ve(ip)=-w*(x-xm)

       case (2) ! gaussian along x
          ue(ip)=w*xl
          ve(ip)=0.
          qe(ip)=a*exp(-(x**2 + y**2)/den)

       case (3) ! gaussian along y
          ue(ip) = 0.
          ve(ip)=w*yl
          qe(ip) = a*exp(-(x**2 +y**2)/den)

       case (4) ! gaussian along x-y diagonal
          ue(ip)=w*xl
          ve(ip)=w*yl
          qe(ip)=a*exp(-(x**2+y**2)/den)

       case (5) !square in ccw diretion
          qe(ip) = 0
          if(r .le. rc) then
             qe(ip) = 1.
          endif
          ue(ip) = w*(y-ym)
          ve(ip) = -w*(x-xm)

       case (6) ! square along x
          qe(ip) = 0.
          if (r .le. rc) then
             qe(ip) = 1.
          endif
          ue(ip) = w*xl
          ve(ip) = 0.

       end select

    end do



  end subroutine exact_solution


end module problem_setup_2d_dg
