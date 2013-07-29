program stuDG
  use lgl_interp
  use meshtools_2d_dg
  use problem_setup_2d_dg
  use lin_alg
  use dg_solvers


  implicit none
  ! Basic input variables
  integer, parameter :: nel = 100
  integer, parameter :: nop = 2
  logical, parameter :: exact_integration = .false.
  integer, parameter :: ntime = 1000
  real(kind=8),parameter :: time_final = 1.0e-0
  integer, parameter :: icase = 2
  integer, parameter :: kstages = 3
  
  integer, parameter:: solver =  1
  ! solver = 1 --> use FORTRAN 90 solver
  ! solver = 2 --> use C/C++ solver
  ! sovler = 3 --> use C/C++/CUDA solver

  ! variables derived from basic input

  real(kind=4) :: dt = time_final/dble(ntime)
  real(kind=4) :: c
  real(kind=4) :: pi
  integer :: noq
  integer :: nelx, nely
  integer :: ngl, nq
  integer :: nelem, npoin, nboun, nside
  integer :: array_stat
  real(kind=4), allocatable, dimension(:) :: xgl, wgl
  real(kind=4), allocatable, dimension(:) :: xnq, wnq
  real(kind=4), allocatable, dimension(:,:) :: psi, dpsi
  real(kind=4),allocatable, dimension(:,:) :: coord
  integer, allocatable, dimension(:,:,:):: intma
  integer, allocatable, dimension(:,:) :: bsido
  integer, allocatable, dimension(:) :: iperiodic 
  real(kind=4),allocatable,dimension(:,:,:) :: ksi_x, ksi_y
  real(kind=4),allocatable,dimension(:,:,:) :: eta_x,eta_y,jac
  integer, allocatable, dimension(:,:) :: iside, jeside,psideh
  integer, allocatable, dimension(:,:,:) :: imapl,imapr
  real(kind=4), allocatable, dimension(:,:) :: nx,ny,jac_side

  ! variables for conversion from tensor-product form
  ! to non-tensor-product form

  real(kind=4),allocatable,dimension(:,:) :: psi2d,psi2d_ksi,psi2d_eta
  real(kind=4),allocatable,dimension(:,:) :: ksi2d_x,ksi2d_y,eta2d_x,eta2d_y
  real(kind=4),allocatable,dimension(:,:) :: jac2d
  integer,allocatable,dimension(:,:) :: intma2d
  integer :: npts, nqs

  real(kind=4),allocatable,dimension(:,:,:) :: Mmatrix,Mmatrix_inv
  real(kind=4),allocatable,dimension(:,:):: Mtemp,Mtemp_inv
  real(kind=4),allocatable,dimension(:) :: rtemp
  ! solution and problem parameter variables
  real(kind=4),allocatable,dimension(:)::qa,ua,va
  real(kind=4) :: time

  ! re-shaped arrays for solver
  real(kind=4),allocatable,dimension(:,:) :: q0,qe,ue,ve


  ! aux variables for stuDG
  integer :: i, j, e, ip,ie
  real(kind=4) :: top, bot, l2_norm;

real(kind=4) :: time_start, time_end, time_per_timestep

character(len=36) :: fmt_string

  ! Begin execution section


  write (*,*) 'Beginning simulation'
  write (*,*) 'Number of elements = ', nel*nel
  write (*,*) 'Order of interpolation = ', nop
  write (*,*) 'Exact integration = ', exact_integration
  write (*,*) 'Time final = ', time_final
  write (*,*) 'Time step size = ', dt
  write (*,*) 'Using ',kstages,' stage RK integration'

  pi = 4.0*atan(1.0)

  select case (icase)
  case (1)
     c = 2.
  case (2)
     c = 1.
  case (3)
     c = 1.
  case (4)
     c = 1.
  case (5)
     c = 2.
  case (6)
     c = 1.

  end select

  if (exact_integration) then
     noq = nop + 1
  else
     noq = nop
  endif

  !ntime = ceiling(time_final/dt)
  nelx = nel
  nely = nel
  nelem = nelx*nely
  ngl = nop + 1
  nq = noq + 1

  npoin = (nop*nelx+1)*(nop*nely+1)
  nboun = 2*nelx + 2*nely
  nside = 2*nelem + nelx + nely

  ! get the interpolation points and weights
  allocate(xgl(ngl),wgl(ngl),stat=array_stat)
  if (array_stat .ne. 0) then
     write (*,*) 'Failure to allocate xgl and wgl arrays'
  endif

write (*,*) 'Computing LGL points and generating Lagrange basis'

  !Compute LGL Points and weights (for interpolation)
  call legendre_gauss_lobatto(xgl,wgl,ngl)


  allocate(xnq(nq),wnq(nq),psi(ngl,nq),dpsi(ngl,nq),stat=array_stat)
  if(array_stat .ne. 0) then
     write (*,*) 'Failure to allocate xnq, wnq, psi, dpsi'
  endif

  !Compute Legendre Cardinal functions and derivatives
  !(for mapping computational domain and integrating)
  call lagrange_basis(psi,dpsi,xnq,wnq,ngl,nq,xgl)

  allocate(coord(npoin,2),intma(nelem,ngl,ngl), &
       bsido(nboun,4), iperiodic(npoin),stat=array_stat)
  if(array_stat .ne. 0) then
     write (*,*) 'Failure to allocate coord,intma,bsido,iperiodic'
  endif

 

write (*,*) 'Generating the mesh'
  ! create the grid  
  call create_grid_2d_dg(coord,intma,bsido, &
       iperiodic, xgl, npoin,nelem,nboun, &
       nelx, nely, ngl)

 allocate(ksi_x(nelem,nq,nq),ksi_y(nelem,nq,nq), &
       eta_x(nelem,nq,nq),eta_y(nelem,nq,nq), &
       jac(nelem,nq,nq),stat=array_stat)


write (*,*) 'Computing Metric terms'
  ! compute metric terms
  call metrics_2d(ksi_x,ksi_y,eta_x,eta_y,jac, &
       coord, intma, psi, dpsi, &
       wnq, npoin,nelem, ngl, nq)

  allocate(iside(nside,4),jeside(nelem,4),stat=array_stat)
  if(array_stat .ne. 0) then
     write (*,*) 'Failure to allocate isde and jeside'
  endif

write (*,*) 'Computing side data'
  ! create side data
  call create_side_2d(iside,jeside,intma,bsido,npoin,nelem, &
       nboun, nside, ngl)

  allocate(psideh(nside,4),imapl(4,2,ngl),imapr(4,2,ngl),&
       stat=array_stat)
  if(array_stat .ne. 0) then
     write (*,*) 'Failure to allocate psideh,imapl,imapr'
  endif

  ! create side data for DG
  call create_side_2d_dg(psideh,imapl,imapr,iside,intma,nside, &
       nelem,ngl)

!!$write(*,*) 'after create_side2d_dg, psideh = '
!!$write(*,'(1X,4I3)') ((psideh(ie,i),i=1,4),ie=1,nside)








  allocate(nx(nside,nq),ny(nside,nq),jac_side(nside,nq),&
       stat=array_stat)
  if(array_stat .ne. 0) then
     write (*,*) 'Failure to allocate nx,ny,jac_side'
  endif

  ! compute normals for flux calculation
  call compute_normals_2d(nx,ny,jac_side,psideh,intma,coord, &
       nside,ngl,nq,wnq,psi,dpsi,nelem,npoin)

write (*,*) 'Set-up for periodic boundary conditions'
  ! periodicity
  call create_periodicity_2d(psideh,iside,coord,nside,nboun,npoin)

!!$write(*,*) 'after create_periodicity_2d, psideh = '
!!$write(*,'(1X,4I3)') ((psideh(ie,i),i=1,4),ie=1,nside)

  npts = ngl*ngl
  nqs = nq*nq
  allocate(psi2d(npts,nqs),psi2d_ksi(npts,nqs),psi2d_eta(npts,nqs), &
       ksi2d_x(nelem,nqs),ksi2d_y(nelem,nqs),eta2d_x(nelem,nqs), &
       eta2d_y(nelem,nqs),jac2d(nelem,nqs),intma2d(nelem,npts),& 
       stat=array_stat)
  if(array_stat .ne. 0) then
     write (*,*) 'Failure to allocate non-tensor-product arrays'
  endif

write (*,*) 'Re-shaping arrays to Non-Tensor Product form'
  ! for solvers that want the non-tensor-product form 
  call reshape_1d_to_2d(psi2d,psi2d_ksi, psi2d_eta, &
       ksi2d_x,ksi2d_y,eta2d_x,eta2d_y, &
       jac2d,intma2d,npts,nqs,psi,dpsi, &
       ksi_x,ksi_y,eta_x,eta_y,jac,intma,&
       nelem,ngl,nq)

  allocate(qa(npoin),ua(npoin),va(npoin), &
       stat=array_stat)
  if(array_stat .ne. 0) then
     write (*,*) 'Failure to allocate solution/parameter arrays'
  endif

  time = 0.
  call exact_solution(qa,ua,va,coord,npoin,time,icase)






  allocate(Mmatrix(nelem,npts,npts),stat=array_stat)
  if(array_stat .ne. 0) then
     write (*,*) 'Failure to allocate Mass Matrix'
  endif
write(*,*) 'Computing the Mass Matrix'
  ! generate the mass matrix
  call create_Mmatrix2d_dg_NonTensorProduct(Mmatrix,jac2d, &
       psi2d,nelem,npts,nqs)

  allocate(Mmatrix_inv(nelem,npts,npts),Mtemp(npts,npts), &
       Mtemp_inv(npts,npts),rtemp(npts),stat=array_stat)

  ! invert the mass matrix
  Mmatrix_inv = 0.
  Mtemp = 0.
  Mtemp_inv = 0.
  rtemp = 0.
write (*,*) 'Inverting the Mass Matrix'
! have suspicion for the matrix invert routine, but it 
! seems to be doing the right thing
  do e = 1,nelem
     Mtemp(:,:) = Mmatrix(e,:,:)
     call diag_matrix_invert(Mtemp,Mtemp_inv,npts)
     !call matrix_invert(Mtemp,Mtemp_inv,npts)
     Mmatrix_inv(e,:,:) = Mtemp_inv
  end do

  allocate(qe(nelem,npts),ue(nelem,npts),ve(nelem,npts),&
       q0(nelem,npts),stat=array_stat)
  q0 = 0.; qe = 0.; ue = 0.; ve = 0.;
  ! re-shape arrays here...
  do e = 1,nelem
     do i = 1,npts
        ip = intma2d(e,i)
        qe(e,i) = qa(ip)
        ue(e,i) = ua(ip)
        ve(e,i) = va(ip)
     end do
  end do

!!$write(fmt_string,*) '(1X,',npts,'ES12.3)'
!!$write (*,*) 'qe:'
!!$write(*,fmt_string) ((qe(i,ie),ie=1,npts),i=1,nelem)



  ! ----- Ready to send data to the solver ---------------------
  
  q0 = qe; 

  write (*,*) 'Set-up complete, commencing time-stepping'

  call cpu_time(time_start)

  select case (solver)

  case (1)

     call dg_advection_2d_non_tensor_product(ntime,&
          dt,kstages,q0,ue,ve,psi,ksi2d_x,ksi2d_y,&
          eta2d_x,eta2d_y,jac2d,psi2d,psi2d_ksi,&
          psi2d_eta,intma2d,psideh,jac_side,imapl,imapr,&
          nx,ny,nelem,npts,nqs,ngl,nq,Mmatrix_inv,nside)

  case (2)
     call dg_advection_2d_ntp_c(ntime,dt,kstages,q0,ue,ve,psi, &
          ksi2d_x,ksi2d_y,eta2d_x,eta2d_y,jac2d,psi2d,psi2d_ksi, &
          psi2d_eta,intma2d,psideh,jac_side,imapl,imapr, &
          nx,ny,nelem,npts,nqs,ngl,nq,Mmatrix_inv,nside)

  case (3)
     call dg_advection_2d_ntp_cuda(ntime,dt,kstages,q0,ue,ve,psi, &
          ksi2d_x,ksi2d_y,eta2d_x,eta2d_y,jac2d,psi2d,psi2d_ksi, &
          psi2d_eta,intma2d,psideh,jac_side,imapl,imapr, &
          nx,ny,nelem,npts,nqs,ngl,nq,Mmatrix_inv,nside)


  end select

  call cpu_time(time_end)
!!$write(fmt_string,*) '(1X,',npts,'ES12.3)'
!!$write(*,*) 'q0 after return = '
!!$write(*,fmt_string) ((q0(i,ie),ie=1,npts),i=1,nelem)

  ! ------ Solver returned, upack the data and check the result ---
  time = dt*ntime ! update the time

  call exact_solution(qa,ua,va,coord,npoin,time,icase)

  ! -- reshape arrays ---
  do e = 1,nelem
     do i = 1,npts
        ip = intma2d(e,i)
        qe(e,i)=qa(ip)
        ue(e,i)=ua(ip)
        ve(e,i)=va(ip)
     end do!i
  end do !e

  ! --- compute norm
  top = 0.; bot = 0.;
  do e = 1,nelem
     do i = 1,npts
        top = top + (q0(e,i)-qe(e,i))**2
        bot = bot + qe(e,i)**2
     end do !i
  end do !e

  l2_norm = sqrt(top/bot)

time_per_timestep = (time_end - time_start)/ntime
write(*,*) 'Time per timestep = ', time_per_timestep

write (*,*) 'l2_norm = ', l2_norm

  ! --- all done


  ! deallocate dynamic arrays
  deallocate(xgl,wgl,xnq,wnq,psi,dpsi,coord,intma, & 
       bsido, iperiodic,ksi_x,ksi_y, eta_x,eta_y, &
       jac,iside,jeside, psideh, imapl, imapr, &
       nx,ny,jac_side,psi2d,psi2d_ksi,psi2d_eta, &
       ksi2d_x,ksi2d_y,eta2d_x,eta2d_y, &
       jac2d,intma2d,qa,ua,va,q0,qe,ue,ve, &
       Mmatrix,Mmatrix_inv,Mtemp,Mtemp_inv, &
       rtemp,stat=array_stat)
  if(array_stat .ne. 0) then
     write (*,*) 'Failure to deallocate dynamic arrays'
  endif


end program stuDG
