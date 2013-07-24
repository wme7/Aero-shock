!   the common block
module data_module
  implicit none
  integer nv,nx,md
  real nh
  parameter (nv=20,nx=2560,md=5,nh=nx/2)

  real u(0:md,-md:nx+md,1:nv),v(0:md,0:nx+1,1:nv),hg(0:md,0:nx+1,1:nv)
  real x(-md:nx+md),dx(-md:nx+md),bl(0:10),br(0:10),bi(0:md)
  real ai(0:10,0:10),am1(0:md),xx(200),uu(200),alim(3,3)
  real aintt(0:10,0:10,0:10)
  real coef(10,10,10),rco(10,2,10),gauss(10,2),coef9(5,5,6)
  real rco9(5,5,6)
  real sigma(2,6),gau(6,2),aii(0:10,0:10),sr
  real gh(1:nv),w(1:nv),vv(1:nv),cc(1:nv)
  real rki(0:md,-md:nx+md),uki(0:md,-md:nx+md),eki(0:md,-md:nx+md),tki(0:md,-md:nx+md),
        pki(0:md,-md:nx+md),ueq(0:md,-md:nx+md,1:nv)
  real pi,cfld,dxmin,xm2
  integer mp,kcmax,indexmax,indexmin
  integer indexnum,indx1(0:nx), index
  real t,dt,mm

  !! Input Param !!
  integer mo     !ooa in space
  integer mt     !ooa in time
  integer kflux  !common flux '1=Roe, 2=LF, 3=LLF'
  real cflc      !CFL number, 0.3 for P1, 0.18 for P2 and 0.1 for P3
  integer ierror !0=initial runs; 1=restart
  real tprint    !input the terminal time
  integer nv     !nv discrete velocity points
  integer n      !n cells
  real xmmm      !TVB constant M

  namelist/proj_list/mo,mt,kflux,cflc,ierror,tprint,n,xmmm
end module data_module
