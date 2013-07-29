! DG&WENO for Classical Botlzmann Equation
! By Manuel Diaz, NTU, 2013.06.16
! manuel.ade'at'gmail.com
!
program main
        ! Load modules
        use IDmodule            ! Create ID names for output files
        use mathmodule          ! Linear Algebra functions
        use dispmodule          ! Matlab 'disp' function for displaying arrays
        use tecplotmodule       ! Write to tecplot functions
        use quadraturemodule    ! Gauss Hermite and Newton cotes abscissas and weights
        use initmodule          ! Load 1d initial condition to algorithm

        ! No implicit definitions allowed
        implicit none

        ! Name of simulation
        character(len=7) :: name = "SBBGK1d"

        ! Define manually inputs
        !integer :: theta = 0
        !integer :: quad = 1
        !integer :: nx = 100
        !integer :: P_deg = 3
        !integer :: RK_stages = 3
        !integer :: IC_case  = 1
        !integer :: fmodel = 1
        !integer :: f_case = 1
        !integer :: method = 1
        !real    :: tau = 1.0/10000

        ! Inputs to be loaded Inputs (comment if parameters where to be setup manualy)
        integer :: theta,nx,quad,P_deg,RK_stages,IC_case,fmodel,f_case,method,kflux
        real    :: tau,tEnd,MM

        ! Constants and loop variables
        real, parameter :: PI = 4*atan(1.0);    ! 3.141529...
        integer :: i,j,l,m

        ! IDs for files
        character(len=100) :: IDn, IDf

        ! Quadrature Parameters / Velocity space points
            ! Gauss Hermite Quadrature (GH)
            integer, parameter :: ngh = 10      ! GH points (actual number of points)
            ! Newton Cotes Quadrature (NC)
            integer, parameter :: nnc = 10      ! NC points (desired number of points)
            integer, parameter :: nodes = 4     ! Number of nodes for NC formula (degree)
            integer, parameter :: nnc_real = 1 + (nodes-1)*ceiling(real(nnc)/(nodes-1))
                                                ! NC points (the actual number of points)
        ! Quadrature weight and abscissas
        integer :: nv                           ! number of velocity points in our problem
        real :: k                               ! quadrature constant
        real, allocatable :: w(:),v(:)          ! Size will depend on the choose quadrature
        real, allocatable :: ww(:,:),vv(:,:)    ! array of velocity points

        ! Array & IC variables
        real, parameter ::x_left=0.0,x_right=1.0! domain boundaries, D:[x_left,x_right]
        real, allocatable :: x(:)               ! the computational x-domain.
        real, allocatable :: xc(:),dx(:)        ! the center nodes of every interval I:[x_{i},x_{i+1}]
        real :: dxmin                           ! smallest interval size.
        real :: rl,ul,pl,rr,ur,pr               ! left and right IC values.
        real, allocatable :: r(:,:),u(:,:)      ! classical macroscopic quantities
        real, allocatable :: p(:,:),t(:,:)      ! r:density, u:x-velocity, p:pressure, t: temperature
        real, allocatable :: rki(:,:),uki(:,:)  ! classical macro-properties as degrees of freedom: (P_deg,nx)
        real, allocatable :: pki(:,:),tki(:,:)  ! r:density, u:x-velocity, p:pressure, t: temperature
        real, allocatable :: f(:,:,:),feq(:,:)  ! PDF arrays

        ! Name list (comment if parameters where to be setup manualy)
        namelist /parameters_list/ name,theta,nx,quad,P_deg,RK_stages,IC_case,fmodel,f_case,method,kflux,tau,tEnd,MM
        ! this list need not to be in order

        !*********************
        ! Main program start:
        !*********************

        ! Read file with parameters
        open(10, file="configuration.in", form="FORMATTED", action="READ")
        read(10, nml=parameters_list)

        ! Checking parametes read
        write(*, nml=parameters_list)

        ! Create IDs for simulation and result file
        call ID_name(name,theta,nx,P_deg,RK_stages,tau,IC_case,fmodel,f_case,method,IDn,IDf)

        ! Checking IDs (uncomment to diplay variables)
        !print *, 'IDn: ',IDn
        !print *, 'IDf: ',IDf

        Print *, 'Build Nv points using Selected a quadrature method'; print *, ' ';
        select case (quad)
            case (1)
                ! allocate w and nv
                allocate( w(ngh) ); allocate( v(ngh) );
                ! compute w and nv
                call gauss_hermite(ngh,v,w); k = 1.0;
                ! add factor e^(x^2) to weights
                w = w*exp(v*v)
                ! Print message
                print *, 'Number of GH points to be used: ',ngh; print *, ' ';
                nv = ngh
            case (2)
                ! allocate w and nv
                allocate( w(nnc_real) ); allocate( v(nnc_real) );
                ! compute w and nv
                call newton_cotes(-2.0,2.0,nnc_real,nodes,v,w,k)
                ! Print message
                print *, 'Number of NC points to be used',nnc_real; print *, ' ';
                nv = nnc_real
            case default
                print *, 'quadrature not available'
        end select

        ! Using DOM
        allocate(vv(nx,nv)); do i = 1,nx; vv(i,:) = v(:); end do; !call disp('vv =',vv,DIGMAX=3,LBOUND=LBOUND(vv))
        allocate(ww(nx,nv)); do i = 1,nx; ww(i,:) = v(:); end do; !call disp('ww =',ww,DIGMAX=3,LBOUND=LBOUND(ww))

        ! Chek point (uncomment to diplay variables)
        !call disp('nv = ',nv); print *, ' ';
        !call disp('w = ',w); print *, ' ';
        !call disp('k = ',k); print *, ' ';

        print *, 'Build Grid and Load initial condition (IC)'; print *, ' ';
        ! Build Grid Domain, define Initial Conditions, and load them into the space domain.
        allocate(x(nx));    x = linspace(x_left,x_right,nx)
        allocate(dx(nx-1)); dx = x(2:nx)-x(1:nx-1); dxmin = minval(dx)
        allocate(xc(nx-1)); xc = x(1:nx-1) + (x(2:nx)-x(1:nx-1))/2
        call Initial_Condition(IC_case,rl,ul,pl,rr,ur,pr)
        allocate(r(nx,nv)); r = rr; r(1:ceiling(real(nx)/2),:) = rl
        allocate(u(nx,nv)); u = ur; u(1:ceiling(real(nx)/2),:) = ul
        allocate(p(nx,nv)); p = pr; p(1:ceiling(real(nx)/2),:) = pl
        allocate(t(nx,nv)); t = 2.0*p/r;

        ! Chek point (uncomment to diplay variables)
        !call disp(' (rho,u,p) at the left  = ',reshape((/rl,ul,pl/),(/1,3/))  )
        !call disp(' (rho,u,p) at the right = ',reshape((/rr,ur,pr/),(/1,3/))  )
        !call disp('u = ',u,DIGMAX=4,LBOUND=LBOUND(u))
        !call disp('vv = ',vv,DIGMAX=3,LBOUND=LBOUND(vv))
        print *, ' ';

        ! Use Macroscopic Variables to proyect into probability distribution space by equilibrium PDF
        allocate(feq(nx,nv)); feq = (r/sqrt(PI*t))*exp(-(vv-u)**2/t); call disp('feq = ',feq,DIGMAX=3,LBOUND=LBOUND(feq))

        ! Use the equilibrium of IC variables as the IC for 'f',
        allocate(f(P_deg,nx,nv)); !f = feq;


end program
