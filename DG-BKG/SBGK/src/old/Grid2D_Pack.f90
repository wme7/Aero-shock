subroutine Init_Physical_Grid_2D
  use MD2D_Grid
  use Legendre
  implicit none

  ! Declare local arguments
  ! Boundary grid for setting physical space grid by transfinite blending
  ! Function : Need to deallocate before leaving this subroutine
  
  real(kind=8), allocatable :: Edge_x(:,:), Edge_y(:,:)
  integer ierr
  integer N, i, j, k
  
  !Linear mapping parameters
  real(kind=8) :: x_start, x_end
  real(kind=8) :: y_start, y_end
  
  integer:: ND_Max
  real(kind=8) :: xi1, xi2 
  
  
  !  Start subroutine
  !  Find the max of degree in all directions
  ND1=PolyDegN_Max(1); ND2=PolyDegN_Max(2);
  ND_Max=max(ND1,ND2,ND3)
  
  allocate( Edge_x(0:ND_Max,1:12), &
       Edge_y(0:ND_Max,1:12), stat= ierr)
  
  if ( ierr .ne. 0 ) then
     write(*,*) 'Grid2D_Pack.f90:'
     write(*,*) 'Can not allocate Edge_x Edge_y and in Grid2D_Pack'
     write(*,*) 'Abort!'
     stop
  endif
  
! Construct Computational Grids by Transfinite Blending Mapping
! Step 1: Construct Points on Egdes (Totally 4 Edges)
! Step 2: Construct Points on Surfaces (Total 6 Surfaces) based on the Edges
!
! Start looping over all domain
  do DDK=1,TotNum_DM
     
     ND1=PolyDegN_DM(1,DDK);
     ND2=PolyDegN_DM(2,DDK);
     
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     !
     ! Set up grid points on Edges 
     !
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     
     ! Edge 1: v1 (-1,-1) to v2 ( 1,-1)
     ! ^^^^^^ 
     x_start = DM_Vertex(1,1,DDK)
     y_start = DM_vertex(1,2,DDK) 
     
     x_end  = DM_Vertex(2,1,DDK)
     y_end  = DM_vertex(2,2,DDK) 
     
     ND=PolyDegN_DM(1,DDK)
     
     Edge_x(0:ND,1) = x_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (x_end-x_start)
     
     Edge_y(0:ND,1) = y_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (y_end-y_start)
     
     !-------------------------------------------------------------!
     
     ! Edge 2: v2 ( 1,-1) to v3 ( 1, 1)
     ! ^^^^^^ 
     x_start = DM_Vertex(2,1,DDK)
     y_start = DM_vertex(2,2,DDK) 
     
     x_end  = DM_Vertex(3,1,DDK)
     y_end  = DM_vertex(3,2,DDK) 
     
     ND=PolyDegN_DM(2,DDK)
     
     Edge_x(0:ND,2) = x_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (x_end-x_start)
     
     Edge_y(0:ND,2) = y_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (y_end-y_start)
     
     !-------------------------------------------------------------!
     
     ! Edge 3: v4 (-1, 1) to v3 ( 1, 1) 
     ! ^^^^^^ 
     x_start = DM_Vertex(4,1,DDK)
     y_start = DM_vertex(4,2,DDK) 
     
     x_end  = DM_Vertex(3,1,DDK)
     y_end  = DM_vertex(3,2,DDK) 
     
     
     ND=PolyDegN_DM(1,DDK)
     
     Edge_x(0:ND,3) = x_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (x_end-x_start)
     
     Edge_y(0:ND,3) = y_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (y_end-y_start)
     
     !-------------------------------------------------------------!
     
     ! Edge 4: v1 (-1,-1) to v4 (-1, 1)
     ! ^^^^^^ 
     x_start = DM_Vertex(1,1,DDK)
     y_start = DM_vertex(1,2,DDK) 
     
     x_end  = DM_Vertex(4,1,DDK)
     y_end  = DM_vertex(4,2,DDK) 
     
     ND=PolyDegN_DM(2,DDK)
     
     Edge_x(0:ND,4) = x_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (x_end-x_start)
     
     Edge_y(0:ND,4) = y_start + (LGLCoord(0:ND,ND)+1.d0)/2.d0 * &
          (y_end-y_start)
     
     !-------------------------------------------------------------!
     
     !-------------------------------------------------------------!
     ! This completes the points on Edges
     !-------------------------------------------------------------!
     
     !=============================================================!
     
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     !
     ! Set up grid points on Surfaces 
     !
     !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
     
     do j=0,ND2
        
        xi2 = LGLCoord(j,ND2)
        
        do i=0,ND1
           
           xi1 = LGLCoord(i,ND1)
           
           ! Original Two D version
           x1(i,j,DDK) = &
                 (1.d0-xi2)/2.d0 * Edge_x(i,1) + &
                 (1.d0+xi2)/2.d0 * Edge_x(i,3) + &
                 (1.d0-xi1)/2.d0 * Edge_x(j,4) + &
                 (1.d0+xi1)/2.d0 * Edge_x(j,2) - &
                 (1.d0-xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(1,1,DDK)- &
                 (1.d0+xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(2,1,DDK)- &
                 (1.d0+xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(3,1,DDK)- &
                 (1.d0-xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(4,1,DDK)

           x2(i,j,DDK) = &
                 (1.d0-xi2)/2.d0 * Edge_y(i,1) + &
                 (1.d0+xi2)/2.d0 * Edge_y(i,3) + &
                 (1.d0-xi1)/2.d0 * Edge_y(j,4) + &
                 (1.d0+xi1)/2.d0 * Edge_y(j,2) - &
                 (1.d0-xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(1,2,DDK)- &
                 (1.d0+xi1)*(1.d0-xi2)/4.d0 * DM_Vertex(2,2,DDK)- &
                 (1.d0+xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(3,2,DDK)- &
                 (1.d0-xi1)*(1.d0+xi2)/4.d0 * DM_Vertex(4,2,DDK)

        enddo ! enddo i
        
     enddo ! enddo j
     
     
     !-------------------------------------------------------------!
     ! This completes the points on Quadalateron
     !-------------------------------------------------------------!
     
     !=============================================================!
     
  enddo ! enddo DDK
  
  
!!$   open(77,file='Surface_Grid.dat')
!!$   write(77,*) 'VARIABLES = "X", "Y", "Z"'
!!$
!!$   do surf_num=5,6
!!$      ND1=PolyDegN_DM(1,1); ND2=PolyDegN_DM(2,1)
!!$      
!!$      write(77,*) 'ZONE I=', ND1+1, ', J=', ND2+1, ', F=POINT' 
!!$
!!$      do j=0,ND2
!!$         do i=0,ND1
!!$            write(77, 10000) Surf_x(i,j,surf_num), &
!!$                             Surf_y(i,j,surf_num), &
!!$                             Surf_z(i,j,surf_num)
!!$         enddo
!!$      enddo
!!$   enddo
!!$   close(77)
!!$
!!$
!!$   open(78,file='LGL_Grid3D.dat')
!!$   write(78,*) 'VARIABLES = "X", "Y", "Z", "V"'
!!$   
!!$   do DDK=1,TotNum_DM
!!$      ND1=PolyDegN_DM(1,DDK)
!!$      ND2=PolyDegN_DM(2,DDK)
!!$      ND3=PolyDegN_DM(3,DDK)
!!$      
!!$      write(78,*) 'ZONE I=', ND1+1, ', J=', ND2+1, ', K=', ND3+1,  'F=POINT'
!!$
!!$      do k=0,ND3
!!$!      write(78,*) 'ZONE I=', ND1+1, ', J=', ND2+1, ', F=POINT'
!!$         do j=0,ND2
!!$            do i=0,ND1
!!$               write(78, 10001) x1(i,j,k,DDK), x2(i,j,k,DDK), x3(i,j,k,DDK), &
!!$                                x3(i,j,k,DDK)
!!$            enddo
!!$         enddo
!!$      enddo
!!$    enddo ! DDK
!!$    close(78) 
!!$   
!!$
!!$
!!$
!!$
!!$
!!$10000       format(3e23.15)
!!$10001       format(4e23.15)


!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!  Deallocate variable Edge_x and Edge_y
  if ( allocated( Edge_x ) ) Deallocate( Edge_x, stat=ierr )
  if ( ierr .ne. 0 ) then
     write(*,*) 'Can not deallocate Bnd_x and in Grid'
     write(*,*) 'Abort'
     stop
  endif
  
  if ( allocated( Edge_y ) ) Deallocate( Edge_y, stat=ierr )
  if ( ierr .ne. 0 ) then
     write(*,*) 'Can not deallocate Bnd_y and in Grid'
     write(*,*) 'Abort'
     stop
  endif
  
  return
  
  
end subroutine Init_Physical_Grid_2D


