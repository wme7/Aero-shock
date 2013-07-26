!   the common block
module data_module
  implicit none
  integer nx,md,nv

  !real*8, allocatable :: ph_q_v(:), ph_q_w(:)
  !real*8 pv_begin,pv_end

  real*8 nh
  parameter (nx=2560,md=5,nh=nx/2,nv=20)

  real*8 u(0:md,-md:nx+md,1:nv),ueq(0:md,-md:nx+md,1:nv)
  real*8 v(0:md,0:nx+1,1:nv),hg(0:md,0:nx+1,1:nv)
  real*8 rki(0:md,-md:nx+md),uki(0:md,-md:nx+md),pki(0:md,-md:nx+md),tki(0:md,-md:nx+md)
  real*8 eki(0:md,-md:nx+md)

  real x(-md:nx+md),dx(-md:nx+md),bl(0:10),br(0:10),bi(0:md)
  real*8 ai(0:10,0:10),am1(0:md),xx(200),uu(200),alim(3,3)
  real*8 aintt(0:10,0:10,0:10)
  real coef(10,10,10),rco(10,2,10),gauss(10,2),coef9(5,5,6)
  real*8 rco9(5,5,6)
  real sigma(2,6),gau(6,2),aii(0:10,0:10),sr
  real*8 gh(1:nv),w(1:nv),cc(1:nv),vv(1:nv)

  real*8 pi,cfld,dxmin,xm2
  integer mp,kcmax,indexmax,indexmin
  integer indexnum,indx1(0:nx), index
  real*8 t,dt,mm

  !! Input Param !!
  integer mo     !ooa in space
  integer mt     !ooa in time
  integer kflux  !common flux '1=Roe, 2=LF, 3=LLF'
  real*8 cflc      !CFL number, 0.3 for P1, 0.18 for P2 and 0.1 for P3
  integer ierror !0=initial runs; 1=restart
  real*8 tprint    !input the terminal time
  integer n      !n cells
  real*8 xmmm      !TVB constant M
  real*8 tau
  integer phase_quadrature_rule
  integer init_value    !init values 0: 0.5+sin(pi*x)  1:sod problem

  namelist/proj_list/mo,mt,kflux,cflc,tau,phase_quadrature_rule,init_value,ierror,tprint,n,xmmm
end module data_module
