      module mod_wave
      integer ::ne,nlevel
      integer ::nproblem,nexist,nbc
      integer ::nsch,imp,npost
      integer ::nx,k,kt,nt,ks,nf
      double precision::t,t0
      double precision::dx,dt,cfl,dt0
      double precision::cc,a,rl
      double precision::theta
      double precision,dimension(100)::tf
      double precision,dimension(:,:),allocatable ::alpha
      double precision,dimension(:),allocatable ::x,w,ddx
      double precision,dimension(:),allocatable ::xt,wt
      double precision,dimension(:,:),allocatable ::pl,plt
      double precision,dimension(:,:),allocatable ::cbin,dl
      double precision,dimension(:,:,:),allocatable ::u,v
      double precision,dimension(:),allocatable ::fluxu,fluxv 
      double precision,dimension(:),allocatable ::c11,c12,c22,c
      double precision,dimension(:),allocatable ::errl1,ordl1,constl1
      double precision,dimension(:),allocatable ::errli,ordli,constli
      double precision,dimension(:),allocatable ::errl2,ordl2,constl2
      double precision,dimension(:),allocatable ::pordli,pordl1,pordl2
      double precision,dimension(:),allocatable ::serrli
      double precision,dimension(:),allocatable ::perrli,perrl1,perrl2 
      double precision,dimension(:),allocatable ::AA
      double precision,dimension(:),allocatable ::f
!      double precision,dimension(:),allocatable ::ffunc
!      double precision,dimension(:,:),allocatable ::phi

      end module mod_wave