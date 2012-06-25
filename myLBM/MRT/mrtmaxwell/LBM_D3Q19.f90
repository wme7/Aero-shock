
      PROGRAM LBM3D
	  use parameters
      implicit none  
	  real*8, external :: Fermi_Bose_func
	  !real(wp), parameter :: piy=3.14159265358979323846264338328_wp
	  !integer, parameter :: Arraysize=81
!                             XMAX,YMAX.ZMAX.grid size in x y dimension   
!-----integer  XMAX,YMAX,ZMAX,NSTART,NUMB,NPLUS,NUMAX
!                                 density........fluid density per link
      integer  XMAX,YMAX,NSTART,NUMB,NPLUS,NUMAX
	  
	  real*8  density
!                                 omega............collision2 parameter
      real*8  omega 
!                                 TMAX.....maximum number of iterations
      integer TMAX
!                                 time................iteration counter
      integer time
!                                 obst(XMAX,YMAX,ZMAX)........obstacle array
!-----logical, dimension(:,:,:), allocatable :: obst !obst(XMAX,YMAX,ZMAX)
!                                 g(0:18,XMAX,YMAX,ZMAX)......fluid densities
!-----real*8, dimension(:,:,:,:), allocatable :: g, gprop  !g(0:18,XMAX,YMAX,ZMAX)
!                                 gprop...rray of temporareYMAX storage
                                 !real*8  gprop(0:18,XMAX,YMAX,ZMAX)
!                                 vel............mean velocity computed
      
	  
	  logical, dimension(:,:), allocatable :: obst !obst(XMAX,YMAX)
	  real*8, dimension(:,:,:), allocatable :: g, gprop  !g(0:8,XMAX,YMAX)

	  real*8  vel,Re

      real*8  UW
	  real*8  original_nu
!=========================================================================================
!================================ MRT ====================================================
!=========================================================================================
real*8, dimension(:,:) , allocatable :: tminv,tm,stmiv,ev,sumcc

real a1
integer e,f,h
!=========================================================================================
!================================ MRT ====================================================
!=========================================================================================

!    ==================================================================

!    |                  Begin initialisation                          |
!    ==================================================================

!
!                       Input calculation parameter

      OPEN(30,file='input0.dat')
      read(30,*)
      read(30,*) TMAX
      read(30,*)
      read(30,*) XMAX
      read(30,*)
      read(30,*) YMAX 
!-----read(30,*)
!-----read(30,*) ZMAX 
      read(30,*)
      read(30,*) NUMAX 
      read(30,*)
      read(30,*) UW
	  read(30,*)
	  read(30,*) Re
	  read(30,*)
	  read(30,*) IF_FERMI
	  read(30,*)
	  read(30,*) T_inf
	  read(30,*)
	  read(30,*) z_inf

!---allocate(obst(XMAX,YMAX,ZMAX))
!---allocate(g(0:18,XMAX,YMAX,ZMAX))
!---allocate(gprop(0:18,XMAX,YMAX,ZMAX))

    allocate(obst(XMAX,YMAX))
    allocate(g(0:8,XMAX,YMAX))
    allocate(gprop(0:8,XMAX,YMAX))

	allocate(tminv(0:8,0:8))
    allocate(tm(0:8,0:8))
    allocate(stmiv(0:8,0:8))
	allocate(ev(0:8,0:8))
    allocate(sumcc(0:8,0:8))
 

!===========================================================!
!!                           Martix M                      !!
!===========================================================!

tm(0,:)=(/1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0/)
tm(1,:)=(/-4.,-1.,-1.,-1.,-1.,2.0,2.0,2.0,2.0/)
tm(2,:)=(/4.0,-2.,-2.,-2.,-2.,1.0,1.0,1.0,1.0/)
tm(3,:)=(/0.0,1.0,0.0,-1.,0.0,1.0,-1.,-1.,1.0/)
tm(4,:)=(/0.0,-2.,0.0,2.0,0.0,1.0,-1.,-1.,1.0/)
tm(5,:)=(/0.0,0.0,1.0,0.0,-1.,1.0,1.0,-1.,-1./)
tm(6,:)=(/0.0,0.0,-2.,0.0,2.0,1.0,1.0,-1.,-1./)
tm(7,:)=(/0.0,1.0,-1.,1.0,-1.,0.0,0.0,0.0,0.0/)
tm(8,:)=(/0.0,0.0,0.0,0.0,0.0,1.0,-1.,1.0,-1./)

!===========================================================!
!!                           Martix M-1                    !!
!===========================================================!

a1=1./36.

tminv(0,:)=(/4.*a1,-4.*a1,4.*a1,0.0,0.0,0.0,0.0,0.0,0.0/)
tminv(1,:)=(/4.*a1,-a1,-2.*a1,6.*a1,-6.*a1,0.0,0.0,9.*a1,0.0/)
tminv(2,:)=(/4.*a1,-a1,-2.*a1,0.0,0.0,6.*a1,-6.*a1,-9.*a1,0.0/)
tminv(3,:)=(/4.*a1,-a1,-2.*a1,-6.*a1,6.*a1,0.0,0.0,9.*a1,0.0/)
tminv(4,:)=(/4.*a1,-a1,-2.*a1,0.0,0.0,-6.*a1,6.*a1,-9.*a1,0.0/)
tminv(5,:)=(/4.*a1,2.*a1,a1,6.*a1,3.*a1,6.*a1,3.*a1,0.0,9.*a1/)
tminv(6,:)=(/4.*a1,2.*a1,a1,-6.*a1,-3.*a1,6.*a1,3.*a1,0.0,-9.*a1/)
tminv(7,:)=(/4.*a1,2.*a1,a1,-6.*a1,-3.*a1,-6.*a1,-3.*a1,0.0,9.*a1/)
tminv(8,:)=(/4.*a1,2.*a1,a1,6.*a1,3.*a1,-6.*a1,-3.*a1,0.0,-9.*a1/)

!===========================================================!
!!                           Martix M-1M                   !!
!===========================================================!

do e=0,8

  do f=0,8

    sumcc(e,f)=0.0

     do h=0,8

       sumcc(e,f)=sumcc(e,f)+tminv(e,h)*tm(h,f)

	 end do

	ev(e,f) = sumcc(e,f)

   end do

end do


	do e=0,8

	   

	 end do 
!===========================================================!
!!                           Martix M-1M                   !!
!===========================================================!	
	
!---density=T_inf**1.5*Fermi_Bose_func(1.5d0,z_inf)
    density=T_inf**1.0d0*Fermi_Bose_func(1.0d0,z_inf)            

!---call wall_coordinate(obst,XMAX,YMAX,ZMAX)     !read obstacle
	call wall_coordinate(obst,XMAX,YMAX)     !read obstacle
!---call init_density(XMAX,YMAX,ZMAX,density,g)
    call init_density(XMAX,YMAX,density,g)

	!call velocity_bc(density,UW,XMAX,YMAX,ZMAX,g) ! FOR TEST

    !omega = 1.0 / ( (XMAX * UW / Re) / T_inf * Fermi_Bose_func(1.5d0,z_inf) / Fermi_Bose_func(2.5d0,z_inf)  +  1.0/2)
	original_nu = XMAX * UW * 3.d0/ Re
!-----------------------------------------------------

print*,'XMAX=',XMAX       !!!!!
print*,'YMAX=',YMAX       !!!!!
print*,'Re=',Re           !!!!!
print*,'UW=',UW           !!!!!
print*,'T_inf=',T_inf     !!!!!
print*,'z_inf=',z_inf     !!!!!
print*,'original_nu=',original_nu !!!!!
!-----------------------------------------------------
	PRINT*,'omega is variable'
!                                   mean sqaure flow is  stored in file
!
!    ==================================================================
!    |                  End initialisation                            |
!    ==================================================================



!    ==================================================================
!    |                  Begin iterations                              |
!    ==================================================================
!      
open(12,file='error.txt')

      NSTART=1
      NUMB=1      
     do  time = NSTART, TMAX
	     IF (NUMB.LE.NUMAX) THEN
	         NUMB=NUMB+1
           END IF
!
		write(*,*) time

!-------call velocity_bc(density,UW,XMAX,YMAX,ZMAX,g)
        call velocity_bc(density,UW,XMAX,YMAX,g)
!                                                    density collision2  
        if (NUMB.GT.NUMAX) then !output results
!--------call write_results(XMAX,YMAX,ZMAX,obst,g,density,time,NUMAX,UW)
		 call write_results(XMAX,YMAX,obst,g,density,time,NUMAX,UW)
         NUMB=1
        end if   
!                                                   density propagation
!-------call streaming(XMAX,YMAX,ZMAX,g,gprop)
		call streaming(XMAX,YMAX,g,gprop)
!                                             bounc back from obstacles
       ! call bounceback(XMAX,YMAX,obst,g,gprop)                           !!!!!!改改改改
			
		!call bounceback(XMAX,YMAX,ZMAX,obst,g,gprop)
!                                           velocity boundary condition 

!-------call collision(density,original_nu,XMAX,YMAX,ZMAX,g,gprop,obst)
		call collision(density,original_nu,XMAX,YMAX,g,gprop,obst,tminv,tm)
!                                                each TMAX/10 iteration
     
	 end do

!    ==================================================================
!    |                    End iterations                              |
!    ==================================================================
!-------call write_results(XMAX,YMAX,ZMAX,obst,g,density,time,NUMAX,UW)
		call write_results(XMAX,YMAX,obst,g,density,time,NUMAX,UW)
      close(10)
	  deallocate(obst,g,gprop)
	  
	  deallocate(tminv,tm,stmiv,ev,sumcc)

      end


!     ****************************************************************
!     *                                                              *
!     *       No obstacles coordinates of a driven cavity flow       *
!     *                                                              *
!     ****************************************************************
!
!-----subroutine wall_coordinate(obst,XMAX,YMAX,ZMAX)
	  subroutine wall_coordinate(obst,XMAX,YMAX)
      implicit none  
!-----integer  XMAX,YMAX,ZMAX
!-----logical  obst(XMAX,YMAX,ZMAX)
	  integer  XMAX,YMAX
      logical  obst(XMAX,YMAX)
!                                                       local variables
!設定六邊為 boundaries, 內部為 fluids, 注意上方頂蓋部分為流体。
!-----integer  x,y,z
	  integer  x,y

!      z=1
!      y=1
        do  y = 1, YMAX
          do  x = 1, XMAX
!-----------do  z = 1, ZMAX

!-------------obst(x,y,z) = .false. !false 為流體
			  obst(x,y) = .false. !false 為流體

!-----------end do
          end do
        end do


      DO x=1, XMAX
!--------DO z=1, ZMAX	         

!-----------obst(x,1   ,z)= .true.
!-----------obst(x,YMAX,z)= .true.
			obst(x,1)= .true.     !下
            
!--------END DO
      END DO


!-----DO x=1, XMAX
         DO y=1, YMAX	         

!----------obst(x,y,1)= .true. !只有下方為固體，上方為流體。
         
		 obst(1,y)= .true.       !左
         obst(XMAX,y)= .true.    !右

         END DO
!-----END DO


!-----DO y=1, YMAX
!-----        DO z=1, ZMAX 	        

!-----            obst(1   ,y,z)= .true.
!-----            obst(XMAX,y,z)= .true.

!-----        END DO
!-----      END DO

      END

!
!
!     ****************************************************************
!     *                                                              *
!     *  Initialize density distribution function n with equilibrium *
!     *  for zero velocity                                           *
!     *                                                              *
!     ****************************************************************
!
!-----      subroutine init_density(XMAX,YMAX,ZMAX,density,g)
      subroutine init_density(XMAX,YMAX,density,g)
	  use parameters
      implicit none
!-----integer  XMAX,YMAX,ZMAX
      integer  XMAX,YMAX
!-----real*8  density,g(0:18,XMAX,YMAX,ZMAX)
	  real*8  density,g(0:8,XMAX,YMAX)
!                                                       local variables
!-----integer  x,y,z
	  integer  x,y
      real*8  t0,t1,t2,g2515
	  real*8, external :: Fermi_Bose_func
!             compute weighting factors (depending on lattice geometry)
!-----t0 = density / 3.0d0
!-----t1 = density / 18.d0
!-----t2 = density / 36.d0
      
      t0 = 4.0d0*density / 9.0d0
      t1 = density / 9.d0
      t2 = density / 36.d0
	   
!      z=1
!      y=1
!                                        loop over computational domain
!-----g2515 = Fermi_Bose_func(2.5d0,z_inf)/Fermi_Bose_func(1.5d0,z_inf)
      g2515 = Fermi_Bose_func(2.0d0,z_inf)/Fermi_Bose_func(1.0d0,z_inf)
      


print*,'XMAX=',XMAX       !!!!!
print*,'YMAX=',YMAX       !!!!!
      do x = 1, XMAX
        do y = 1, YMAX
!---------DO z = 1, ZMAX
!----write(*,*) z
  write(*,*)

!----           g(0,x,y,z) =t0 *(1.d0 + 0.5*(T_inf*g2515-1)*(-3))
!----           g(1,x,y,z) =t1
!----           g(2,x,y,z) =t1
!----           g(3,x,y,z) =t1
!----           g(4,x,y,z) =t1
!----           g(5,x,y,z) =t1
!----           g(6,x,y,z) =t1
!----           g(7,x,y,z) =t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(8,x,y,z) =t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(9,x,y,z) =t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(10,x,y,z)=t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(11,x,y,z)=t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(12,x,y,z)=t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(13,x,y,z)=t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(14,x,y,z)=t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(15,x,y,z)=t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(16,x,y,z)=t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(17,x,y,z)=t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
!----           g(18,x,y,z)=t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*3.d0)
          
          g2515 = 1.0d0
		  T_inf = 1.0d0

          g(0,x,y) =t0 *(1.d0 + 0.5*(T_inf*g2515-1)*(-2))
          
		  g(1,x,y) =t1 *(1.d0 + 0.5*(T_inf*g2515-1.d0))
          g(2,x,y) =t1 *(1.d0 + 0.5*(T_inf*g2515-1.d0))
          g(3,x,y) =t1 *(1.d0 + 0.5*(T_inf*g2515-1.d0))
          g(4,x,y) =t1 *(1.d0 + 0.5*(T_inf*g2515-1.d0))
          
          g(7,x,y) =t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*4.d0)
          g(8,x,y) =t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*4.d0)
          g(5,x,y) =t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*4.d0)
          g(6,x,y) =t2 *(1.d0 + 0.5*(T_inf*g2515-1.d0)*4.d0)
        

!---write(*,*) g(0,1,1,z)

	write(*,*)
!------END DO
        end do
	end do

!test!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
write(*,*) density

write(*,*)
      end


!     ****************************************************************
!     *                                                              *
!     *                Velocity bounday condition                    *
!     *                                                              *
!     ****************************************************************
!
!-----subroutine velocity_bc(density,UW,XMAX,YMAX,ZMAX,g)
	  subroutine velocity_bc(density,UW,XMAX,YMAX,g)
	  use parameters
      implicit none

!-----integer  XMAX,YMAX,ZMAX
	  integer  XMAX,YMAX
!-----real*8   density,g(0:18,XMAX,YMAX,ZMAX),UW
	  real*8   density,g(0:8,XMAX,YMAX),UW

!-----integer  x,y,z,i
	  integer  x,y,i
!-----real*8  t0,t1,t2,U,V,W,un(18),usqu,d_loc,temperature_loc,z_loc
	  real*8  t0,t1,t2,U,V,un(8),usqu,d_loc,temperature_loc,z_loc
      real*8  dt0,dt1,dt2,g2515  !g2515 = Fermi_Bose_func(2.5d0,z_loc)/Fermi_Bose_func(1.5d0,z_loc)
	  real*8, external :: Fermi_Bose_func

!
!-----t0 = 1.d0 /  3.d0
!-----t1 = 1.d0 / 18.d0
!-----t2 = 1.d0 / 36.d0

      t0 = 4.0d0 /  9.0d0
      t1 = 1.0d0 /  9.0d0
      t2 = 1.0d0 /  36.0d0

!                                                 square speed of sound


      d_loc = density

      dt0=t0*d_loc
      dt1=t1*d_loc
      dt2=t2*d_loc
	  temperature_loc=T_inf
	  
	  
	  z_loc=z_inf

!-----g2515 = Fermi_Bose_func(2.5d0,z_loc)/Fermi_Bose_func(1.5d0,z_loc)
      g2515 = Fermi_Bose_func(2.0d0,z_loc)/Fermi_Bose_func(1.0d0,z_loc)
			U =  UW
            V =  0.0d0
!-----------W =  0.0d0
!
!----            usqu = (U * U + V * V + W * W)
!----            un(1) =   U * dsqrt(3.d0)
!----            un(2) =  -U * dsqrt(3.d0)
!----            un(3) =   V * dsqrt(3.d0)
!----            un(4) = - V * dsqrt(3.d0)
!----            un(5) =   W * dsqrt(3.d0)
!----            un(6) = - W * dsqrt(3.d0)
!----            un(7) =   (U + V) * dsqrt(3.d0)
!----            un(8) =   (U - V) * dsqrt(3.d0)      
!----	         un(9) = (- U + V) * dsqrt(3.d0)
!----            un(10)= (- U - V) * dsqrt(3.d0)
!----            un(11)= (  U + W) * dsqrt(3.d0)
!----            un(12)= (- U + W) * dsqrt(3.d0)
!----            un(13)= (  U - W) * dsqrt(3.d0)
!----            un(14)= (- U - W) * dsqrt(3.d0)
!----            un(15)= (  V + W) * dsqrt(3.d0)
!----            un(16)= (  V - W) * dsqrt(3.d0)
!----            un(17)= (- V + W) * dsqrt(3.d0)
!----            un(18)= (- V - W) * dsqrt(3.d0)

            usqu = (U * U + V * V)
            un(1) =   U * dsqrt(3.d0)
            un(3) =  -U * dsqrt(3.d0)
            un(2) =   V * dsqrt(3.d0)
            un(4) =  -V * dsqrt(3.d0)

            un(5) =   (U + V) * dsqrt(3.d0)
            un(8) =   (U - V) * dsqrt(3.d0)      
	        un(6) = (- U + V) * dsqrt(3.d0)
            un(7) = (- U - V) * dsqrt(3.d0)
            

!----z = ZMAX
     y = YMAX 
	  do x = 1, XMAX
!-------do y = 1, YMAX
!----       g(0, x,y,z)= dt0*(1.d0 - 0.5d0 * usqu + 0.5 * (temperature_loc*g2515-1.d0)*(-3.d0) )
!----       g(1, x,y,z)= dt1*(1.d0+un(1) +0.5*un(1) **2.d0-0.5*usqu) 
!----       g(2, x,y,z)= dt1*(1.d0+un(2) +0.5*un(2) **2.d0-0.5*usqu) 
!----       g(3, x,y,z)= dt1*(1.d0+un(3) +0.5*un(3) **2.d0-0.5*usqu) 
!----       g(4, x,y,z)= dt1*(1.d0+un(4) +0.5*un(4) **2.d0-0.5*usqu) 
!----       g(5, x,y,z)= dt1*(1.d0+un(5) +0.5*un(5) **2.d0-0.5*usqu)
!----       g(6, x,y,z)= dt1*(1.d0+un(6) +0.5*un(6) **2.d0-0.5*usqu)
!----       g(7, x,y,z)= dt2*(1.d0+un(7) +0.5*un(7) **2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(8, x,y,z)= dt2*(1.d0+un(8) +0.5*un(8) **2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(9, x,y,z)= dt2*(1.d0+un(9) +0.5*un(9) **2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(10, x,y,z)=dt2*(1.d0+un(10)+0.5*un(10)**2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(11, x,y,z)=dt2*(1.d0+un(11)+0.5*un(11)**2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(12, x,y,z)=dt2*(1.d0+un(12)+0.5*un(12)**2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(13, x,y,z)=dt2*(1.d0+un(13)+0.5*un(13)**2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(14, x,y,z)=dt2*(1.d0+un(14)+0.5*un(14)**2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(15, x,y,z)=dt2*(1.d0+un(15)+0.5*un(15)**2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(16, x,y,z)=dt2*(1.d0+un(16)+0.5*un(16)**2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(17, x,y,z)=dt2*(1.d0+un(17)+0.5*un(17)**2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
!----       g(18, x,y,z)=dt2*(1.d0+un(18)+0.5*un(18)**2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*3.d0 )
      g2515 = 1.0d0
	  temperature_loc = 1.0d0

      g(0, x,y)= dt0*(1.d0 - 0.5d0 * usqu + 0.5 * (temperature_loc*g2515-1.d0)*(-2.d0) )
      g(1, x,y)= dt1*(1.d0+un(1) +0.5*un(1) **2.d0-0.5*usqu+ 0.5 * (temperature_loc*g2515-1.d0)) 
      g(3, x,y)= dt1*(1.d0+un(3) +0.5*un(3) **2.d0-0.5*usqu+ 0.5 * (temperature_loc*g2515-1.d0)) 
      g(2, x,y)= dt1*(1.d0+un(2) +0.5*un(2) **2.d0-0.5*usqu+ 0.5 * (temperature_loc*g2515-1.d0)) 
      g(4, x,y)= dt1*(1.d0+un(4) +0.5*un(4) **2.d0-0.5*usqu+ 0.5 * (temperature_loc*g2515-1.d0)) 

      g(5, x,y)= dt2*(1.d0+un(5) +0.5*un(5) **2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*4.d0 )
      g(8, x,y)= dt2*(1.d0+un(8) +0.5*un(8) **2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*4.d0 )
      g(6, x,y)= dt2*(1.d0+un(6) +0.5*un(6) **2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*4.d0 )
      g(7, x,y)= dt2*(1.d0+un(7) +0.5*un(7) **2.d0-0.5*usqu + 0.5 * (temperature_loc*g2515-1.d0)*4.d0 )


!--------end do
       end do
!--g(0, XMAX,1:YMAX,1:ZMAX)=dt0 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*(-3.d0))
!--g(1:6, XMAX,1:YMAX,1:ZMAX)=dt1
!--g(7:18, XMAX,1:YMAX,1:ZMAX)=dt2 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*3.d0)

!--g(0, 1,1:YMAX,1:ZMAX)=dt0 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*(-3.d0))
!--g(1:6, 1,1:YMAX,1:ZMAX)=dt1 
!--g(7:18, 1,1:YMAX,1:ZMAX)=dt2 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*3.d0)

!--g(0, 1:XMAX,1:YMAX,1)=dt0 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*(-3.d0))
!--g(1:6, 1:XMAX,1:YMAX,1)=dt1 
!--g(7:18, 1:XMAX,1:YMAX,1)=dt2 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*3.d0)

!--g(0, 1:XMAX,1,1:ZMAX)=dt0 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*(-3.d0))
!--g(1:6, 1:XMAX,1,1:ZMAX)=dt1
!--g(7:18, 1:XMAX,1,1:ZMAX)=dt2 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*3.d0)

!--g(0, 1:XMAX,YMAX,1:ZMAX)=dt0 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*(-3.d0))
!--g(1:6, 1:XMAX,YMAX,1:ZMAX)=dt1 
!--g(7:18, 1:XMAX,YMAX,1:ZMAX)=dt2 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*3.d0)
      g2515 = 1.0d0
	  temperature_loc = 1.0d0

g(0, XMAX,1:YMAX)=dt0 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*(-2.d0))
g(1:4, XMAX,1:YMAX)=dt1 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0))
g(5:8, XMAX,1:YMAX)=dt2 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*4.d0)

g(0, 1,1:YMAX)=dt0 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*(-2.d0))
g(1:4, 1,1:YMAX)=dt1 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)) 
g(5:8, 1,1:YMAX)=dt2 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*4.d0)

g(0, 1:XMAX,1)=dt0 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*(-2.d0))
g(1:4, 1:XMAX,1)=dt1 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0))
g(5:8, 1:XMAX,1)=dt2 * (1.d0 + 0.5 * (temperature_loc*g2515-1.d0)*4.d0)

      end

!
!
!     ****************************************************************
!     *                                                              *
!     *   streaming2 fluid densities to their next neighbour nodes    *
!     *                                                              *
!     ****************************************************************
!
!-----subroutine streaming(XMAX,YMAX,ZMAX,g,gprop)
      subroutine streaming(XMAX,YMAX,g,gprop)
      implicit none
!-----integer  XMAX,YMAX,ZMAX
      integer  XMAX,YMAX
!-----real*8   g(0:18,XMAX,YMAX,ZMAX),gprop(0:18,XMAX,YMAX,ZMAX)
	  real*8   g(0:8,XMAX,YMAX),gprop(0:8,XMAX,YMAX)

!-----integer  x,y,z,x_e,x_w,y_n,y_s,z_t,z_b
	  integer  x,y,x_e,x_w,y_n,y_s
!      z=1
!      y=1

      do x = 1, XMAX
        do y = 1, YMAX
!---------do z = 1, ZMAX

!--          x_e = mod(x,XMAX) + 1
!--          y_n = mod(y,YMAX) + 1
!--          z_t = mod(z,ZMAX) + 1
!--          x_w = XMAX - mod(XMAX + 1 - x, XMAX)
!--          y_s = YMAX - mod(YMAX + 1 - y, YMAX)
!--          z_b = ZMAX - mod(ZMAX + 1 - z, ZMAX)

!--          gprop(0, x  ,y  ,z  ) = g(0, x,y,z)
!--          gprop(1, x_e,y  ,z  ) = g(1, x,y,z)
!--          gprop(2, x_w,y  ,z  ) = g(2, x,y,z)
!--          gprop(3, x  ,y_n,z  ) = g(3, x,y,z)
!--          gprop(4, x  ,y_s,z  ) = g(4, x,y,z)
!--		     gprop(5, x  ,y,  z_t) = g(5, x,y,z)
!--          gprop(6, x  ,y,  z_b) = g(6, x,y,z)
!--          gprop(7, x_e,y_n,z  ) = g(7, x,y,z)
!--          gprop(8, x_e,y_s,z  ) = g(8, x,y,z)
!--          gprop(9, x_w,y_n,z  ) = g(9, x,y,z)
!--          gprop(10,x_w,y_s,z  ) = g(10,x,y,z)
!--          gprop(11,x_e,y,  z_t) = g(11,x,y,z)
!--          gprop(12,x_w,y,  z_t) = g(12,x,y,z)
!--          gprop(13,x_e,y,  z_b) = g(13,x,y,z)
!--          gprop(14,x_w,y,  z_b) = g(14,x,y,z)
!--          gprop(15,x,  y_n,z_t) = g(15,x,y,z)
!--          gprop(16,x,  y_n,z_b) = g(16,x,y,z)
!--          gprop(17,x,  y_s,z_t) = g(17,x,y,z)
!--          gprop(18,x,  y_s,z_b) = g(18,x,y,z)

		  x_e = mod(x,XMAX) + 1     !右一
          y_n = mod(y,YMAX) + 1     !上一
          
          x_w = XMAX - mod(XMAX + 1 - x, XMAX)   !左一
          y_s = YMAX - mod(YMAX + 1 - y, YMAX)   !下一
          

          gprop(0, x  ,y    ) = g(0, x,y)
          gprop(1, x_e,y    ) = g(1, x,y)
          gprop(3, x_w,y    ) = g(3, x,y)
          gprop(2, x  ,y_n  ) = g(2, x,y)
          gprop(4, x  ,y_s  ) = g(4, x,y)

          gprop(5, x_e,y_n  ) = g(5, x,y)
          gprop(8, x_e,y_s  ) = g(8, x,y)
          gprop(6, x_w,y_n  ) = g(6, x,y)
          gprop(7, x_w,y_s  ) = g(7, x,y)
          

!----------end do
        end do
      end do
      end

!
!
!     ****************************************************************
!     *                                                              *
!     *  Fluid densities are rotated. By the next propagation step,  *
!     *    this results in a bounce back from obstacle nodes.        *
!     *                                                              *
!     ****************************************************************
!
!-----subroutine bounceback(XMAX,YMAX,ZMAX,obst,g,gprop)
	  subroutine bounceback(XMAX,YMAX,obst,g,gprop)
      implicit none
!-----integer  XMAX,YMAX,ZMAX
	  integer  XMAX,YMAX
!-----logical  obst(XMAX,YMAX,ZMAX)
	  logical  obst(XMAX,YMAX)
!-----real*8   g(0:18,XMAX,YMAX,ZMAX),gprop(0:18,XMAX,YMAX,ZMAX)
	  real*8   g(0:8,XMAX,YMAX),gprop(0:8,XMAX,YMAX)
!                                                       local variables
!-----integer  x,y,z
	  integer  x,y

!      z=1
!      y=1
!流體和固體都作 bounce back

      do x = 1, XMAX
        do y = 1, YMAX
!----------do z= 1, ZMAX
!-----------if (obst(x,y,z)) then
			if (obst(x,y)) then

!--            g(1, x,y,z) = gprop(2, x,y,z)
!--            g(2, x,y,z) = gprop(1, x,y,z)
!--            g(3, x,y,z) = gprop(4, x,y,z)
!--            g(4, x,y,z) = gprop(3, x,y,z)
!--            g(5, x,y,z) = gprop(6, x,y,z)
!--            g(6, x,y,z) = gprop(5, x,y,z)
!--            g(7, x,y,z) = gprop(10,x,y,z)
!--            g(8, x,y,z) = gprop(9, x,y,z)
!--            g(9, x,y,z) = gprop(8, x,y,z)
!--            g(10,x,y,z) = gprop(7, x,y,z)
!--            g(11,x,y,z) = gprop(14,x,y,z)
!--            g(12,x,y,z) = gprop(13,x,y,z)
!--            g(13,x,y,z) = gprop(12,x,y,z)
!--            g(14,x,y,z) = gprop(11,x,y,z)
!--            g(15,x,y,z) = gprop(18,x,y,z)
!--            g(16,x,y,z) = gprop(17,x,y,z)
!--            g(17,x,y,z) = gprop(16,x,y,z)
!--            g(18,x,y,z) = gprop(15,x,y,z)


            g(1, x,y) = gprop(3, x,y)
            g(3, x,y) = gprop(1, x,y)
            g(2, x,y) = gprop(4, x,y)
            g(4, x,y) = gprop(2, x,y)

            g(5, x,y) = gprop(7, x,y)
            g(8, x,y) = gprop(6, x,y)
            g(6, x,y) = gprop(8, x,y)
            g(7, x,y) = gprop(5, x,y)
            

            end if     
!---------end do
        end do
      end do

      end
!
!
!     ****************************************************************
!     *                                                              *
!     *          One-step density collision2 process                 *
!     *                                                              *
!     ****************************************************************
!
!-----subroutine collision(density,original_nu,XMAX,YMAX,ZMAX,g,gprop,obst)
	  subroutine collision(density,original_nu,XMAX,YMAX,g,gprop,obst,tminv,tm)
      implicit none
	  real*8, external :: Fermi_Bose_func
	  real*8, external :: f
	  real*8, external :: df
!-----integer  XMAX,YMAX,ZMAX
!-----logical  obst(XMAX,YMAX,ZMAX)
!-----real*8   density,g(0:18,XMAX,YMAX,ZMAX),gprop(0:18,XMAX,YMAX,ZMAX)
	  integer  XMAX,YMAX
      logical  obst(XMAX,YMAX)
      real*8   density
      real*8   g(0:8,XMAX,YMAX)
	  real*8   gprop(0:8,XMAX,YMAX)
!
!                                                       local variables
!-----integer  x,y,z,i
	  integer  x,y,i
!-----real*8  t0,t1,t2,U,V,W,un(18),gq(0:18),usqu,d_loc,energy_loc,z_loc,temperature_loc
!-----real*8  dt0,dt1,dt2,omega,omega1,g2515,original_nu
	  real*8  t0,t1,t2,U,V,un(8),usqu,d_loc,energy_loc,z_loc,temperature_loc
      real*8  dt0,dt1,dt2,omega,g2515,original_nu
!====================================================================================================
real*8  sm(0:8) , stmiv(0:8,0:8),tminv(0:8,0:8),tm(0:8,0:8),tu
integer a,b,c,d,l
!real*8  fmom(0:8,XMAX,YMAX),suma(0:8,XMAX,YMAX),sumb(0:8,XMAX,YMAX)
real*8  gq(0:8,XMAX,YMAX)

real*8 , dimension(:,:,:) , allocatable :: fmom,sumb
!====================================================================================================
!--      t0 = 1.d0 /  3.d0
!--      t1 = 1.d0 / 18.d0
!--      t2 = 1.d0 / 36.d0

	  t0 = 4.d0 /  9.d0
      t1 = 1.d0 /  9.d0
      t2 = 1.d0 / 36.d0

	  !omega = 1.d0 / ( original_nu +  1.0/2)
      !t0 = omega*t0
      !t1 = omega*t1
      !t2 = omega*t2
	  
      !omega1=1.0d0-omega

!       z=1
!       y=1
!$acc region
      do x = 1, XMAX
        do y = 1, YMAX
!-----     do z = 1, ZMAX
!                                   only free nodes are considered here
!-----    if (.not. obst(x,y,z)) then
		  if (.not. obst(x,y)) then
!                                         calculate local density d_loc
!-      d_loc = gprop(0, x,y,z) &
!-     &      + gprop(1, x,y,z) + gprop(2, x,y,z)&
!-     &      + gprop(3, x,y,z) + gprop(4, x,y,z)&
!-     &      + gprop(5, x,y,z) + gprop(6, x,y,z)&
!-     &      + gprop(7, x,y,z) + gprop(8, x,y,z)&
!-     &      + gprop(9, x,y,z) + gprop(10,x,y,z)&
!-     &      + gprop(11,x,y,z) + gprop(12,x,y,z)&
!-     &      + gprop(13,x,y,z) + gprop(14,x,y,z)&
!-     &      + gprop(15,x,y,z) + gprop(16,x,y,z)&
!-     &      + gprop(17,x,y,z) + gprop(18,x,y,z)


	 d_loc = gprop(0, x,y) &
     &      + gprop(1, x,y) + gprop(2, x,y)&
     &      + gprop(3, x,y) + gprop(4, x,y)&
     &      + gprop(5, x,y) + gprop(6, x,y)&
     &      + gprop(7, x,y) + gprop(8, x,y)


!                                         calculate velocity components


		U =  gprop(1, x,y) - gprop(3, x,y)& 
     &     + gprop(5, x,y) + gprop(8, x,y)&
     &     - gprop(6, x,y) - gprop(7,x,y)

        U = dsqrt(3.d0)*U  / d_loc

        V =  gprop(2, x,y) - gprop(4, x,y)& 
     &     + gprop(5, x,y) - gprop(8, x,y)&
     &     + gprop(6, x,y) - gprop(7, x,y)

        V = dsqrt(3.d0)*V  / d_loc

            usqu = (U * U + V * V)
            un(1) =   U * dsqrt(3.d0)
            un(3) = - U * dsqrt(3.d0)
            un(2) =   V * dsqrt(3.d0)
            un(4) = - V * dsqrt(3.d0)
            
            un(5) =   (U + V) * dsqrt(3.d0)
            un(8) =   (U - V) * dsqrt(3.d0)      
	        un(6) = (- U + V) * dsqrt(3.d0)
            un(7) = (- U - V) * dsqrt(3.d0)


      energy_loc = -4.0d0*gprop(0, x,y) -1.0d0 *(gprop(1, x,y) + gprop(2, x,y)+ gprop(3, x,y) + gprop(4, x,y))+2.d0*gprop(5, x,y) + 2.d0*gprop(6, x,y)+ 2.d0*gprop(7, x,y) + 2.d0*gprop(8, x,y) 
      

!---call newton(z_loc,d_loc,U,V,W,energy_loc)
	call newton(z_loc,d_loc,U,V,energy_loc)
!---g2515 = Fermi_Bose_func(2.5d0,z_loc)/Fermi_Bose_func(1.5d0,z_loc)
    
	g2515 = Fermi_Bose_func(2.0d0,z_loc)/Fermi_Bose_func(1.0d0,z_loc)
!---temperature_loc = (energy_loc/d_loc - 0.5*usqu) / g2515 * 2 / 3
	
    g2515 = 1.0d0
	

	
	temperature_loc = (energy_loc/(2.0d0*d_loc) +2.0d0 - 0.5*usqu) / g2515
    
	temperature_loc = 1.0d0
!omega = 1.0 / ( (XMAX * UW / Re) / T_inf * Fermi_Bose_func(1.5d0,z_loc) / Fermi_Bose_func(2.5d0,z_loc)  +  1.0/2)
omega = 1.d0 / ( original_nu / temperature_loc / g2515  +  1.0/2)
!write(*,*) omega
call matrix1(omega,tminv,stmiv)

!========================================================================================
!================================     M-1S ==============================================
!========================================================================================
!========================================================================================
!tu=1.0*omega

!sm(0:8)=(/1.0d0,1.4d0,1.4d0,1.0d0,1.2d0,1.0d0,1.2d0 , tu , tu /)

!do a=0,8

 !   do b=0,8

!	stmiv(a,b)=tminv(a,b)*sm(b)

!	end do

!end do 

!do a=0,8
  
 ! print *, (stmiv(a,b),b=0,8)

!end do
!========================================================================================
!========================================================================================
!========================================================================================
!========================================================================================

      g2515 = 1.0d0
	  temperature_loc = 1.0d0


   gq(0,x,y)=d_loc

   gq(1,x,y)=d_loc*(-4.0+1.0*(U*U+V*V)+2.0*(temperature_loc*g2515))

   gq(2,x,y)=d_loc*(3.0-1.0*(U*U+V*V)-2.0*(temperature_loc*g2515))

   gq(3,x,y)=d_loc*U/dsqrt(3.d0)

   gq(4,x,y)=-d_loc*U/dsqrt(3.d0)

   gq(5,x,y)=d_loc*V/dsqrt(3.d0)

   gq(6,x,y)=-d_loc*V/dsqrt(3.d0)

   gq(7,x,y)=d_loc*(U*U-V*V)/3.0d0

   gq(8,x,y)=d_loc*U*V/3.0d0

	allocate(fmom(0:8,XMAX,YMAX))

    allocate(sumb(0:8,XMAX,YMAX))

!do c=0,8

!	suma(c,x,y)=0.0

!	   do l=0,8

!	   suma(c,x,y)=suma(c,x,y)+tm(c,l)*gprop(l,x,y)

!	   end do
    
!	fmom(c,x,y)=suma(c,x,y)

!fmom(c,x,y)=tm(c,0)*gprop(0,x,y)+tm(c,1)*gprop(1,x,y)+tm(c,2)*gprop(2,x,y)+tm(c,3)*gprop(3,x,y)  &
!&  +tm(c,4)*gprop(4,x,y)+tm(c,5)*gprop(5,x,y)+tm(c,6)*gprop(6,x,y)+tm(c,7)*gprop(7,x,y)+tm(c,8)*gprop(8,x,y)
!end do

fmom(1,x,y)=tm(1,0)*gprop(0,x,y)+tm(1,1)*gprop(1,x,y)+tm(1,2)*gprop(2,x,y)+tm(1,3)*gprop(3,x,y)  &
&  +tm(1,4)*gprop(4,x,y)+tm(1,5)*gprop(5,x,y)+tm(1,6)*gprop(6,x,y)+tm(1,7)*gprop(7,x,y)+tm(1,8)*gprop(8,x,y)

fmom(0,x,y)=tm(0,0)*gprop(0,x,y)+tm(0,1)*gprop(1,x,y)+tm(0,2)*gprop(2,x,y)+tm(0,3)*gprop(3,x,y)  &
&  +tm(0,4)*gprop(4,x,y)+tm(0,5)*gprop(5,x,y)+tm(0,6)*gprop(6,x,y)+tm(0,7)*gprop(7,x,y)+tm(0,8)*gprop(8,x,y)

fmom(2,x,y)=tm(2,0)*gprop(0,x,y)+tm(2,1)*gprop(1,x,y)+tm(2,2)*gprop(2,x,y)+tm(2,3)*gprop(3,x,y)  &
&  +tm(2,4)*gprop(4,x,y)+tm(2,5)*gprop(5,x,y)+tm(2,6)*gprop(6,x,y)+tm(2,7)*gprop(7,x,y)+tm(2,8)*gprop(8,x,y)

fmom(3,x,y)=tm(3,0)*gprop(0,x,y)+tm(3,1)*gprop(1,x,y)+tm(3,2)*gprop(2,x,y)+tm(3,3)*gprop(3,x,y)  &
&  +tm(3,4)*gprop(4,x,y)+tm(3,5)*gprop(5,x,y)+tm(3,6)*gprop(6,x,y)+tm(3,7)*gprop(7,x,y)+tm(3,8)*gprop(8,x,y)

fmom(4,x,y)=tm(4,0)*gprop(0,x,y)+tm(4,1)*gprop(1,x,y)+tm(4,2)*gprop(2,x,y)+tm(4,3)*gprop(3,x,y)  &
&  +tm(4,4)*gprop(4,x,y)+tm(4,5)*gprop(5,x,y)+tm(4,6)*gprop(6,x,y)+tm(4,7)*gprop(7,x,y)+tm(4,8)*gprop(8,x,y)

fmom(5,x,y)=tm(5,0)*gprop(0,x,y)+tm(5,1)*gprop(1,x,y)+tm(5,2)*gprop(2,x,y)+tm(5,3)*gprop(3,x,y)  &
&  +tm(5,4)*gprop(4,x,y)+tm(5,5)*gprop(5,x,y)+tm(5,6)*gprop(6,x,y)+tm(5,7)*gprop(7,x,y)+tm(5,8)*gprop(8,x,y)

fmom(6,x,y)=tm(6,0)*gprop(0,x,y)+tm(6,1)*gprop(1,x,y)+tm(6,2)*gprop(2,x,y)+tm(6,3)*gprop(3,x,y)  &
&  +tm(6,4)*gprop(4,x,y)+tm(6,5)*gprop(5,x,y)+tm(6,6)*gprop(6,x,y)+tm(6,7)*gprop(7,x,y)+tm(6,8)*gprop(8,x,y)

fmom(7,x,y)=tm(7,0)*gprop(0,x,y)+tm(7,1)*gprop(1,x,y)+tm(7,2)*gprop(2,x,y)+tm(7,3)*gprop(3,x,y)  &
&  +tm(7,4)*gprop(4,x,y)+tm(7,5)*gprop(5,x,y)+tm(7,6)*gprop(6,x,y)+tm(7,7)*gprop(7,x,y)+tm(7,8)*gprop(8,x,y)

fmom(8,x,y)=tm(8,0)*gprop(0,x,y)+tm(8,1)*gprop(1,x,y)+tm(8,2)*gprop(2,x,y)+tm(8,3)*gprop(3,x,y)  &
&  +tm(8,4)*gprop(4,x,y)+tm(8,5)*gprop(5,x,y)+tm(8,6)*gprop(6,x,y)+tm(8,7)*gprop(7,x,y)+tm(8,8)*gprop(8,x,y)


!	  sumb(c,x,y)=0.00

!	     do l=0,8

!		 sumb(c,x,y)=sumb(c,x,y)+stmiv(c,l)*(fmom(l,x,y)-gq(l,x,y))

!		 end do

sumb(0,x,y) = stmiv(0,1)*(fmom(1,x,y)-gq(1,x,y))+stmiv(0,2)*(fmom(2,x,y)-gq(2,x,y))  &
&          +stmiv(0,3)*(fmom(3,x,y)-gq(3,x,y))+stmiv(0,4)*(fmom(4,x,y)-gq(4,x,y))  &
&          +stmiv(0,5)*(fmom(5,x,y)-gq(5,x,y))+stmiv(0,6)*(fmom(6,x,y)-gq(6,x,y))  &
&          +stmiv(0,7)*(fmom(7,x,y)-gq(7,x,y))+stmiv(0,8)*(fmom(8,x,y)-gq(8,x,y))  &
&          +stmiv(0,0)*(fmom(0,x,y)-gq(0,x,y))  

sumb(1,x,y) = stmiv(1,1)*(fmom(1,x,y)-gq(1,x,y))+stmiv(1,2)*(fmom(2,x,y)-gq(2,x,y))  &
&          +stmiv(1,3)*(fmom(3,x,y)-gq(3,x,y))+stmiv(1,4)*(fmom(4,x,y)-gq(4,x,y))  &
&          +stmiv(1,5)*(fmom(5,x,y)-gq(5,x,y))+stmiv(1,6)*(fmom(6,x,y)-gq(6,x,y))  &
&          +stmiv(1,7)*(fmom(7,x,y)-gq(7,x,y))+stmiv(1,8)*(fmom(8,x,y)-gq(8,x,y))  &
&          +stmiv(1,0)*(fmom(0,x,y)-gq(0,x,y))  

sumb(2,x,y) = stmiv(2,1)*(fmom(1,x,y)-gq(1,x,y))+stmiv(2,2)*(fmom(2,x,y)-gq(2,x,y))  &
&          +stmiv(2,3)*(fmom(3,x,y)-gq(3,x,y))+stmiv(2,4)*(fmom(4,x,y)-gq(4,x,y))  &
&          +stmiv(2,5)*(fmom(5,x,y)-gq(5,x,y))+stmiv(2,6)*(fmom(6,x,y)-gq(6,x,y))  &
&          +stmiv(2,7)*(fmom(7,x,y)-gq(7,x,y))+stmiv(2,8)*(fmom(8,x,y)-gq(8,x,y))  &
&          +stmiv(2,0)*(fmom(0,x,y)-gq(0,x,y))  

sumb(3,x,y) = stmiv(3,1)*(fmom(1,x,y)-gq(1,x,y))+stmiv(3,2)*(fmom(2,x,y)-gq(2,x,y))  &
&          +stmiv(3,3)*(fmom(3,x,y)-gq(3,x,y))+stmiv(3,4)*(fmom(4,x,y)-gq(4,x,y))  &
&          +stmiv(3,5)*(fmom(5,x,y)-gq(5,x,y))+stmiv(3,6)*(fmom(6,x,y)-gq(6,x,y))  &
&          +stmiv(3,7)*(fmom(7,x,y)-gq(7,x,y))+stmiv(3,8)*(fmom(8,x,y)-gq(8,x,y))  &
&          +stmiv(3,0)*(fmom(0,x,y)-gq(0,x,y))  

sumb(4,x,y) = stmiv(4,1)*(fmom(1,x,y)-gq(1,x,y))+stmiv(4,2)*(fmom(2,x,y)-gq(2,x,y))  &
&          +stmiv(4,3)*(fmom(3,x,y)-gq(3,x,y))+stmiv(4,4)*(fmom(4,x,y)-gq(4,x,y))  &
&          +stmiv(4,5)*(fmom(5,x,y)-gq(5,x,y))+stmiv(4,6)*(fmom(6,x,y)-gq(6,x,y))  &
&          +stmiv(4,7)*(fmom(7,x,y)-gq(7,x,y))+stmiv(4,8)*(fmom(8,x,y)-gq(8,x,y))  &
&          +stmiv(4,0)*(fmom(0,x,y)-gq(0,x,y))  

sumb(5,x,y) = stmiv(5,1)*(fmom(1,x,y)-gq(1,x,y))+stmiv(5,2)*(fmom(2,x,y)-gq(2,x,y))  &
&          +stmiv(5,3)*(fmom(3,x,y)-gq(3,x,y))+stmiv(5,4)*(fmom(4,x,y)-gq(4,x,y))  &
&          +stmiv(5,5)*(fmom(5,x,y)-gq(5,x,y))+stmiv(5,6)*(fmom(6,x,y)-gq(6,x,y))  &
&          +stmiv(5,7)*(fmom(7,x,y)-gq(7,x,y))+stmiv(5,8)*(fmom(8,x,y)-gq(8,x,y))  &
&          +stmiv(5,0)*(fmom(0,x,y)-gq(0,x,y))  

sumb(6,x,y) = stmiv(6,1)*(fmom(1,x,y)-gq(1,x,y))+stmiv(6,2)*(fmom(2,x,y)-gq(2,x,y))  &
&          +stmiv(6,3)*(fmom(3,x,y)-gq(3,x,y))+stmiv(6,4)*(fmom(4,x,y)-gq(4,x,y))  &
&          +stmiv(6,5)*(fmom(5,x,y)-gq(5,x,y))+stmiv(6,6)*(fmom(6,x,y)-gq(6,x,y))  &
&          +stmiv(6,7)*(fmom(7,x,y)-gq(7,x,y))+stmiv(6,8)*(fmom(8,x,y)-gq(8,x,y))  &
&          +stmiv(6,0)*(fmom(0,x,y)-gq(0,x,y))  

sumb(7,x,y) = stmiv(7,1)*(fmom(1,x,y)-gq(1,x,y))+stmiv(7,2)*(fmom(2,x,y)-gq(2,x,y))  &
&          +stmiv(7,3)*(fmom(3,x,y)-gq(3,x,y))+stmiv(7,4)*(fmom(4,x,y)-gq(4,x,y))  &
&          +stmiv(7,5)*(fmom(5,x,y)-gq(5,x,y))+stmiv(7,6)*(fmom(6,x,y)-gq(6,x,y))  &
&          +stmiv(7,7)*(fmom(7,x,y)-gq(7,x,y))+stmiv(7,8)*(fmom(8,x,y)-gq(8,x,y))  &
&          +stmiv(7,0)*(fmom(0,x,y)-gq(0,x,y))  

sumb(8,x,y) = stmiv(8,1)*(fmom(1,x,y)-gq(1,x,y))+stmiv(8,2)*(fmom(2,x,y)-gq(2,x,y))  &
            &          +stmiv(8,3)*(fmom(3,x,y)-gq(3,x,y))+stmiv(8,4)*(fmom(4,x,y)-gq(4,x,y))  &
            &          +stmiv(8,5)*(fmom(5,x,y)-gq(5,x,y))+stmiv(8,6)*(fmom(6,x,y)-gq(6,x,y))  &
            &          +stmiv(8,7)*(fmom(7,x,y)-gq(7,x,y))+stmiv(8,8)*(fmom(8,x,y)-gq(8,x,y))  &
            &          +stmiv(8,0)*(fmom(0,x,y)-gq(0,x,y))  	  
	 
	 
	  g(0,x,y)=gprop(0,x,y)-sumb(0,x,y)
	  g(1,x,y)=gprop(1,x,y)-sumb(1,x,y)
	  g(2,x,y)=gprop(2,x,y)-sumb(2,x,y)
	  g(3,x,y)=gprop(3,x,y)-sumb(3,x,y)
	  g(4,x,y)=gprop(4,x,y)-sumb(4,x,y)
	  g(5,x,y)=gprop(5,x,y)-sumb(5,x,y)
	  g(6,x,y)=gprop(6,x,y)-sumb(6,x,y)
	  g(7,x,y)=gprop(7,x,y)-sumb(7,x,y)
	  g(8,x,y)=gprop(8,x,y)-sumb(8,x,y)
     
deallocate(fmom,sumb)

          end if
!------   end do
        end do
      end do
	  !write(*,*) omega, 1.0 / (original_nu +  1.0/2)
!$acc end region

      end
!


!     ****************************************************************
!     *                                                              *
!     *           Output of rsults to file 'result.dat'              *
!     *                                                              *
!     ****************************************************************
!
!-----subroutine write_results(XMAX,YMAX,ZMAX,obst,g,density,time,NUMAX,UW)
	  subroutine write_results(XMAX,YMAX,obst,g,density,time,NUMAX,UW)
	  use parameters
      implicit none  
	  real*8, external :: Fermi_Bose_func
!-----integer  XMAX,YMAX,ZMAX,time,NUMAX
!-----real*8  g(0:18,XMAX,YMAX,ZMAX),density,UW
!-----logical  obst(XMAX,YMAX,ZMAX)
	  integer  XMAX,YMAX,time,NUMAX
      real*8  g(0:8,XMAX,YMAX),density,UW
      logical  obst(XMAX,YMAX)
!                                                       local variables
!-----integer  x,y,z,i
!-----real*8  U,V,W,d_loc,press,c_squ,energy_loc,z_loc,temperature_loc,usqu,g2515
	  integer  x,y,i
	  integer  ax,by
      real*8  U,V,d_loc,press_xx,press_xy,c_squ,energy_loc,z_loc,temperature_loc,usqu,g2515
	  real*8  U1,U2,U3,V1,V2,V3
      real*8  UT1,UT2,UT3
	  real*8  UA1,UA2,UA3,VA1,VA2,VA3,error
      CHARACTER *6 rr

!                                              open results output file
!      rr=character(time)

      print*,'time=',time

         WRITE (rr,200) 100000+time !/NUMAX
200        FORMAT (I6) 

!--open(9,file='centerlineX00'//TRIM(rr)//'.plt')
!--open(10,file='centerline00Z'//TRIM(rr)//'.plt')
open(11,file='Comp'//TRIM(rr)//'.plt')

!
!=====================write header for postprocessing with TECPLOT================
!--      write(9,*) 'VARIABLES = X, W' 
!--     write(9,*) 'ZONE J=',YMAX,', F=POINT'
!=================================================================================
!=====================write header for postprocessing with TECPLOT================
!--       write(10,*) 'VARIABLES = Z, U' 
!--       write(10,*) 'ZONE I=',XMAX,', F=POINT'
!=================================================================================
!=====================write header for postprocessing with TECPLOT================
!--      write(11,*) 'VARIABLES = X, Y, Z, U, V, W, Fugasity, Temperature, PRESS' 
!--      write(11,*) 'ZONE I=',XMAX,', J=', YMAX,', K=',ZMAX,', F=POINT'

	  write(11,*) 'VARIABLES = X, Y, U, V, Fugasity, Temperature, PRESS_XX , PRESS_XY ' 
      write(11,*) 'ZONE I=',XMAX,', J=', YMAX,', F=POINT'



!      y=20


!---   do z = 1, ZMAX
         do y = 1, YMAX
            do x = 1, XMAX
!                                         calculate local density d_loc
!---      d_loc = g(0, x,y,z)&
!---     &      + g(1, x,y,z) + g(2, x,y,z)&
!---     &      + g(3, x,y,z) + g(4, x,y,z)&
!---     &      + g(5, x,y,z) + g(6, x,y,z)&
!---     &      + g(7, x,y,z) + g(8, x,y,z)&
!---     &      + g(9, x,y,z) + g(10,x,y,z)&
!---     &      + g(11,x,y,z) + g(12,x,y,z)&
!---     &      + g(13,x,y,z) + g(14,x,y,z)&
!---     &      + g(15,x,y,z) + g(16,x,y,z)&
!---     &      + g(17,x,y,z) + g(18,x,y,z)

	       d_loc = g(0, x,y)&
     &      + g(1, x,y) + g(2, x,y)&
     &      + g(3, x,y) + g(4, x,y)&
     &      + g(5, x,y) + g(6, x,y)&
     &      + g(7, x,y) + g(8, x,y)
    


!
!--        U =  g(1, x,y,z) - g(2, x,y,z)& 
!--     &     + g(7, x,y,z) + g(8, x,y,z)&
!--     &     - g(9, x,y,z) - g(10,x,y,z)&
!--     &     + g(11,x,y,z) - g(12,x,y,z)&
!--     &     + g(13,x,y,z) - g(14,x,y,z)
!--        U = dsqrt(3.d0)*U/ d_loc

        U =  g(1, x,y) - g(3, x,y)& 
     &     + g(5, x,y) + g(8, x,y)&
     &     - g(6, x,y) - g(7, x,y)

        U = dsqrt(3.d0)*U/ d_loc


!
!--        V =  g(3, x,y,z) - g(4, x,y,z)& 
!--     &     + g(7, x,y,z) - g(8, x,y,z)&
!--     &     + g(9, x,y,z) - g(10,x,y,z)&
!--     &     + g(15,x,y,z) + g(16,x,y,z)&
!--     &     - g(17,x,y,z) - g(18,x,y,z)    
!--        V = dsqrt(3.d0)*V/ d_loc 

	    V =  g(2, x,y) - g(4, x,y)& 
     &     + g(5, x,y) - g(8, x,y)&
     &     + g(6, x,y) - g(7, x,y)
    
        V = dsqrt(3.d0)*V/ d_loc


!--        W =  g(5, x,y,z) - g(6, x,y,z)&
!--     &     + g(11,x,y,z) + g(12,x,y,z)&
!--     &     - g(13,x,y,z) - g(14,x,y,z)&
!--     &     + g(15,x,y,z) - g(16,x,y,z)&
!--     &     + g(17,x,y,z) - g(18,x,y,z)      
!--        W = dsqrt(3.d0)*W/ d_loc 
!
!--    usqu = (U * U + V * V + W * W)
	   usqu = (U * U + V * V)


!--	   energy_loc =  1.5 *( &
!--     &      + g(1, x,y,z) + g(2, x,y,z)&
!--     &      + g(3, x,y,z) + g(4, x,y,z)&
!--     &      + g(5, x,y,z) + g(6, x,y,z) )+3.0*(&
!--     &      + g(7, x,y,z) + g(8, x,y,z)&
!--     &      + g(9, x,y,z) + g(10,x,y,z)&
!--     &      + g(11,x,y,z) + g(12,x,y,z)&
!--     &      + g(13,x,y,z) + g(14,x,y,z)&
!--     &      + g(15,x,y,z) + g(16,x,y,z)&
!--     &      + g(17,x,y,z) + g(18,x,y,z) )


  energy_loc = -4.0d0*g(0, x,y) -1.0d0 *(g(1, x,y) + g(2, x,y)+ g(3, x,y) + g(4, x,y))+2.d0*g(5, x,y) + 2.d0*g(6, x,y)+ 2.d0*g(7, x,y) + 2.d0*g(8, x,y)
	
!---call newton(z_loc,d_loc,U,V,W,energy_loc)
	call newton(z_loc,d_loc,U,V,energy_loc)
!---g2515 = Fermi_Bose_func(2.5d0,z_loc)/Fermi_Bose_func(1.5d0,z_loc)
    g2515 = Fermi_Bose_func(2.0d0,z_loc)/Fermi_Bose_func(1.0d0,z_loc)
!---temperature_loc = (energy_loc/d_loc - 0.5*usqu) / g2515 * 2 / 3
    
	temperature_loc = (energy_loc/(2.0d0*d_loc) +2.0d0 - 0.5*usqu) / g2515

print*,'temperature_loc =',temperature_loc  

      g2515 = 1.0d0
	  temperature_loc = 1.0d0

print*,'temperature_loc =',temperature_loc  

	press_xx = 1.0d0 *(g(1, x,y) - g(2, x,y)+ g(3, x,y) - g(4, x,y))
    press_xy = 1.0d0 *(g(5, x,y) - g(6, x,y)+ g(7, x,y) - g(8, x,y))
!                                                 write results to file
!----write(11,100) x,y,z,U,V,W,z_loc,temperature_loc,press
	 write(11,100) x,y,U,V,z_loc,temperature_loc,press_xx,press_xy

!--	 if(z==ZMAX/2 .and. y==YMAX/2) then
!--	  write(9,99) -0.5d0 + 1.d0/(XMAX-1)*(x-1), W / UW /2
!--	 endif

!--	 if (x==XMAX/2 .and. y==YMAX/2) then
!--	  write(10,99) -0.5d0 + 1.d0/(ZMAX-1)*(z-1), U / UW /2
!--	 endif

99        format(2E15.6)
100       format(2I8,6E15.6)

          end do
        end do
!--   end do

!                                                close file 'result.dat'
!---close(9)
!---close(10)
close(11)

!---------------------------收斂後停止條件----------------------------
!---------------------------------------------------------------------
!---------------------------------------------------------------------
   ax=((XMAX-1)/4)
   by=((YMAX-1)/4)
   
   
   
        U1 =  (g(1, ax,by) - g(3, ax,by)+ g(5, ax,by) + g(8, ax,by)- g(6, ax,by) - g(7, ax,by))  &
		&     /(g(0,  ax,by)+ g(1, ax,by) + g(2, ax,by)+ g(3, ax,by) + g(4, ax,by)+ g(5, ax,by) + g(6, ax,by)+ g(7, ax,by) + g(8, ax,by))


	    V1 =  (g(2,ax,by) - g(4, ax,by)+ g(5,ax,by) - g(8,ax,by)+ g(6,ax,by) - g(7,ax,by))&
		&     /(g(0,  ax,by)+ g(1, ax,by) + g(2, ax,by)+ g(3, ax,by) + g(4, ax,by)+ g(5, ax,by) + g(6, ax,by)+ g(7, ax,by) + g(8, ax,by))
  
		U3 =  (g(1, 3.0*ax,3.0*by)-g(3,3.0*ax,3.0*by)+ g(5, 3.0*ax,3.0*by)+g(8,3.0*ax,3.0*by)- g(6, 3.0*ax,3.0*by)-g(7,3.0*ax,3.0*by))   &  
		&     /(g(0,3.0*ax,3.0*by)+ g(1,3.0*ax,3.0*by) + g(2,3.0*ax,3.0*by)+g(3,3.0*ax,3.0*by) + g(4,3.0*ax,3.0*by)&
		&     + g(5,3.0*ax,3.0*by)+g(6,3.0*ax,3.0*by)+ g(7,3.0*ax,3.0*by)+ g(8,3.0*ax,3.0*by))


	    V3 =  (g(2,3.0*ax,3.0*by) - g(4,3.0*ax,3.0*by)+ g(5,3.0*ax,3.0*by)-g(8,3.0*ax,3.0*by)+ g(6,3.0*ax,3.0*by) - g(7,3.0*ax,3.0*by))&
		&     /(g(0,3.0*ax,3.0*by)+ g(1,3.0*ax,3.0*by) + g(2,3.0*ax,3.0*by)+ g(3,3.0*ax,3.0*by) + g(4,3.0*ax,3.0*by)&
		&     + g(5,3.0*ax,3.0*by)+ g(6,3.0*ax,3.0*by)+ g(7,3.0*ax,3.0*by) + g(8,3.0*ax,3.0*by))
    
        U2 =  (g(1,2.0*ax,2.0*by) - g(3,2.0*ax,2.0*by)+ g(5,2.0*ax,2.0*by)      +g(8,2.0*ax,2.0*by)- g(6,2.0*ax,2.0*by) - g(7,2.0*ax,2.0*by))      /(g(0,2.0*ax,2.0*by)+ g(1,2.0*ax,2.0*by) + g(2,2.0*ax,2.0*by)     + g(3,2.0*ax,2.0*by) + g(4,2.0*ax,2.0*by)+ g(5,2.0*ax,2.0*by)      + g(6,2.0*ax,2.0*by)+ g(7,2.0*ax,2.0*by) + g(8,2.0*ax,2.0*by))


	    V2 =  (g(2,2.0*ax,2.0*by) - g(4,2.0*ax,2.0*by)+ g(5,2.0*ax,2.0*by)      -g(8,2.0*ax,2.0*by)+ g(6,2.0*ax,2.0*by) - g(7,2.0*ax,2.0*by))    /(g(0,2.0*ax,2.0*by)+ g(1,2.0*ax,2.0*by) + g(2,2.0*ax,2.0*by)     +g(3,2.0*ax,2.0*by) + g(4,2.0*ax,2.0*by)+ g(5,2.0*ax,2.0*by)      +g(6,2.0*ax,2.0*by)+ g(7,2.0*ax,2.0*by) + g(8,2.0*ax,2.0*by))

    Ut1=U1**2+V1**2

	Ut2=U2**2+V2**2

	UT3=U3**2+V3**2
    error = (U1-UA1)**2+(U2-UA2)**2+(U3-UA3)**2+(V1-VA1)**2+(V2-VA2)**2+(V3-VA3)**2/(Ut1+Ut2+Ut3)

101    format(I10,E15.6)
	    
		write(12,101) time , error
	
	if ( time >= 10000 .and. ( error <= 1e-10  )) then

       close (12)
	   stop
    endif
	
	UA1=U1
	UA2=U2
	UA3=U3
	VA1=V1
	VA2=V2
	VA3=V3
	

      end


!================================ Tools for Semiclassical===============================
!=======================================================================================
FUNCTION Fermi_Bose_func(v,z) !calculate Fermi or Bose function value
 use parameters
 implicit none
 real*8, INTENT(IN):: v ,z
 real*8 :: Fermi_Bose_func ,term_order
 integer :: l ,steps

 Fermi_Bose_func = 0
 term_order = 1.d0
 l = 1
 steps = 0
if(IF_Fermi) then !Fermi

  do while(abs(term_order) >= 1e-6)
  term_order = ((-z)**l) / (dble(l)**v)
  Fermi_Bose_func = Fermi_Bose_func + term_order
  l = l+1
  steps = steps + 1
  enddo
  Fermi_Bose_func = -Fermi_Bose_func

else !Bose

  do while(abs(term_order) > 1e-6)
   term_order = (z**l) / (dble(l)**v)
   Fermi_Bose_func = Fermi_Bose_func + term_order
   l = l+1
   steps = steps + 1
  enddo

endif
END FUNCTION Fermi_Bose_func

!---subroutine newton(z,rho,u0,u1,u2,E)
subroutine newton(z,rho,u0,u1,E) 
use parameters   
implicit none                                        
!---real*8, intent (in) :: rho,u0,u1,u2,E 
real*8, intent (in) :: rho,u0,u1,E                                         
real*8, intent (out) :: z
real*8 ,external :: f ,df

z = z_inf
!---do while( abs(f(z,rho,u0,u1,u2,E)) > 1e-6 )
do while( abs(f(z,rho,u0,u1,E)) > 1e-6 )
   
!---z = z - f(z,rho,u0,u1,u2,E)/df(z,rho,u0,u1,u2,E)
   z = z - f(z,rho,u0,u1,E)/df(z,rho,u0,u1,E)
   z = dabs(z)
print*,'z=',z 
   if(z > 1.d0 ) then
	 print*,'z=',z
	 z = 0.9d0
     
	 write(*,*) 'Fugacity is too big !!'
   endif

enddo


end subroutine newton


!---function f(z,rho,u0,u1,u2,E)
function f(z,rho,u0,u1,E)
implicit none                                        
real*8, intent (in) :: z,rho,u0,u1,E
!---real*8, intent (in) :: z,rho,u0,u1,u2,E                                         
real*8 :: f
real*8, external :: Fermi_Bose_func
!---f = 2*E - 3*(rho/Fermi_Bose_func(1.5d0,z))**(5.d0/3)*Fermi_Bose_func(2.5d0,z) - rho*u0*u0 - rho*u1*u1 - rho*u2*u2
!---f = 2*E - 3*(rho/Fermi_Bose_func(1.0d0,z))**(5.d0/3)*Fermi_Bose_func(2.0d0,z) - rho*u0*u0 - rho*u1*u1
f = E - 2*(rho/Fermi_Bose_func(1.0d0,z))**(2.0d0)*Fermi_Bose_func(2.0d0,z) - 3.0d0*rho*u0*u0 - 3.0d0*rho*u1*u1 + 4.0d0*rho

end function f

!---function df(z,rho,u0,u1,u2,E)
function df(z,rho,u0,u1,E)
implicit none
!---real*8, intent (in) :: z,rho,u0,u1,u2,E        
real*8, intent (in) :: z,rho,u0,u1,E                                   
real*8 :: df
real*8, external :: Fermi_Bose_func
!---df = (5 * rho**(5.d0/3)) / z * Fermi_Bose_func(2.0d0,z) * Fermi_Bose_func(0.0d0,z) * Fermi_Bose_func(1.0d0,z) ** (-8.d0/3) - (3 * rho**(5.0/3)) / z * Fermi_Bose_func(1.0d0,z) ** (-2.0/3)
df = (4 * rho**(2.0d0)) / z * Fermi_Bose_func(2.0d0,z) * Fermi_Bose_func(0.0d0,z) * Fermi_Bose_func(1.0d0,z) ** (-3.0d0) - (2 * rho**(2.0d0)) / z * Fermi_Bose_func(1.0d0,z) ** (-1.0d0)
  
end function df



subroutine matrix1(omega,tminv,stmiv) 
use parameters
implicit none
real*8 sm(0:8),tu
real*8, intent (in) :: omega,tminv(0:8,0:8)                                         
real*8, intent (out) :: stmiv(0:8,0:8)

integer a,b


tu=1.0*omega

sm(0:8)=(/1.0d0,1.4d0,1.4d0,1.0d0,1.2d0,1.0d0,1.2d0 , tu , tu /)

do a=0,8

    do b=0,8

	stmiv(a,b)=tminv(a,b)*sm(b)

	end do

end do 

do a=0,8
  
  

end do
end subroutine matrix1

